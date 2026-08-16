# Internal: the Track B differentially private synthesiser for linked
# multi-table data. Not exported.
#
# The privacy unit is the ROOT entity (e.g. a patient). Removing one entity
# removes its single root row and, cascading down the key hierarchy, all of its
# descendant rows in every child table. The release is therefore accounted at
# the root-entity grain across the whole database at once.
#
# Contribution bounding is hierarchical. `max_rows_per_person` gives a *branching
# cap* per child table (the maximum children kept per parent unit); the root cap
# is always 1 (one root row per entity). Capping top-down bounds how much of each
# table one entity can touch: the per-entity row cap for a table T is the product
# of the branching caps along the path root -> T (`path_cap`), and the number of
# T's parent units one entity owns is `path_cap[parent(T)]`.
#
# Under those caps the whole release is one concatenation of histograms:
#   * root table  : variable marginals (one-way, + pairwise for a Chow-Liu tree),
#                   person-sensitivity path_cap[root] = 1 each;
#   * child table : a children-per-parent count histogram (person-sensitivity
#                   path_cap[parent], one bin per possible count 0..cap), plus its
#                   own variable marginals at person-sensitivity path_cap[T];
#   * numeric bin edges are DP-estimated per table at sensitivity path_cap[T]
#     (domain = "dp"), sharing the domain_frac slice.
# Because histograms concatenate, the budget composes exactly (see
# `dp_make_noise()`): total L1 sensitivity and summed squared L2 sensitivity add
# over every histogram, and uniform per-cell noise calibrated to those totals is
# (eps, delta)-DP end to end.
#
# Generation is parent-first: draw the root entities from the root's noisy
# marginals, then for each synthetic parent draw a child count from the count
# model, copy the parent's (surrogate) key down as the foreign key so referential
# integrity holds by construction, regenerate the own index as a within-parent
# position, and fill the child's variables from its own noisy marginals. Drawing
# `m` datasets from the fitted models is post-processing (no extra budget).
#
# Cross-table conditioning is opt-in via dp_control(cross_table = TRUE): a child
# table's variables are then conditioned on the synthetic parent's attributes by
# measuring parent-by-child joint marginals (at the child's path-cap sensitivity)
# and folding the parent variables into the child's Chow-Liu structure as fixed
# context nodes (see dp-linked-cross.R). With cross_table = FALSE (default) child
# variables are modelled by their own within-table marginals - the synthetic child
# links to a synthetic parent but is statistically independent of it.
#
# Within-table longitudinal structure is opt-in via dp_control(longitudinal = ...):
# a child table's repeated rows are then modelled as a DP Markov trajectory over
# the child's own key index (an initial-state model + per-variable transition
# matrices), with the children-per-parent count model doubling as the length
# model, so within-unit autocorrelation survives (see dp-linked-longitudinal.R).
# When cross_table = TRUE is set together with a longitudinal model on the same
# table the two combine: the table's INITIAL-STATE model is cross-conditioned on
# the synthetic parent (parent-by-child joints at the first-row sensitivity), and
# the transition chain then carries that parent dependence across the trajectory.
# The transitions themselves stay parent-free, so the extra cost is exactly the
# initial-state joints.
#
# A longitudinal child also carries the flat DP Markov engine's two transition
# controls, applied per table: dp_control(baseline = ...) names subject-invariant
# columns that are modelled once in the initial state and then held constant within
# a unit (so they contribute no transition histogram, sharpening the rest at the
# same budget), and dp_control(transition_order / transition_cross) deepen each
# time-varying column's transition to its own last `order` values plus the lag-1
# values of its `cross` most associated companions. Extra conditioning is
# budget-neutral (a tuple lands in one cell); a higher order lowers the transition
# sensitivity to path_cap[parent] * (branching_cap - order).
#
# Remaining scope (documented in the DP vignette): only the IMMEDIATE parent
# conditions a child (deeper ancestors reach it through the parent's synthesised
# values). Constraints and `unit = "row"` are refused.

# Resolve the branching cap for one table from dp$max_rows_per_person, which is
# either a single integer (applied to every child table) or a named integer
# vector keyed by table name. Returns NULL when unset for this table.
dp_branch_cap <- function(dp, table) {
  m <- dp$max_rows_per_person
  if (is.null(m)) return(NULL)
  if (!is.null(names(m))) {
    v <- m[[table]]
    if (is.null(v)) NULL else as.integer(v)
  } else {
    as.integer(m)
  }
}

# Per-table local branching caps and per-entity path caps (product root -> T).
dp_linked_caps <- function(hierarchy, dp) {
  tn <- hierarchy$names
  local <- stats::setNames(rep(NA_integer_, length(tn)), tn)
  for (t in tn) {
    if (is.na(hierarchy$parent[[t]])) {
      local[[t]] <- 1L                       # one root row per entity
    } else {
      cap <- dp_branch_cap(dp, t)
      if (is.null(cap) || cap < 1L) {
        stop(sprintf(paste0(
          "linked DP: no branching cap for child table '%s'. Set ",
          "`max_rows_per_person` in dp_control() to the public maximum children ",
          "per parent - a single integer for every child table, or a named ",
          "list like list(%s = 5) per table."), t, t), call. = FALSE)
      }
      local[[t]] <- cap
    }
  }
  path <- stats::setNames(rep(NA_integer_, length(tn)), tn)
  for (t in hierarchy$order) {
    p <- hierarchy$parent[[t]]
    path[[t]] <- if (is.na(p)) local[[t]] else path[[p]] * local[[t]]
  }
  list(local = local, path = path)
}

# The modelled variables of a table: everything except its key columns and any
# structural (nesting) index. Keys are copied / regenerated, not modelled.
dp_table_vars <- function(hierarchy, t, tables) {
  structural <- union(hierarchy$keys[[t]], hierarchy$structure[[t]]$nested)
  setdiff(names(tables[[t]]), structural)
}

# Hierarchical, top-down contribution bounding. For each child table (parents
# first) drop rows orphaned by an already-capped parent, then keep at most
# `cap_t` children per surviving parent unit. Over-cap units are trimmed at
# random, except for longitudinal tables (`longi[[t]]` TRUE) where the
# temporally-first `cap_t` rows (by the own key index) are kept as a contiguous
# prefix so consecutive-pair transitions are not corrupted. Guarantees each entity
# owns <= path_cap[T] rows of T and <= path_cap[parent] units of parent(T).
# Returns the capped tables and per-table dropped counts.
dp_cap_hierarchy <- function(tables, hierarchy, caps, longi = NULL) {
  capped  <- tables
  dropped <- stats::setNames(integer(length(hierarchy$names)), hierarchy$names)
  for (t in hierarchy$order) {
    p <- hierarchy$parent[[t]]
    if (is.na(p)) next                       # root: key unique -> 1 row / entity
    fk    <- hierarchy$fk[[t]]
    cap_t <- caps$local[[t]]
    is_long <- isTRUE(longi[[t]])
    own <- hierarchy$own[[t]]
    ck <- key_string(capped[[t]], fk)
    pk <- key_string(capped[[p]], fk)        # surviving parent units
    idx <- which(ck %in% pk)                 # drop orphans from earlier capping
    groups <- split(idx, ck[idx])
    keep <- unlist(lapply(groups, function(rows) {
      if (length(rows) <= cap_t) rows
      else if (is_long) {                    # keep the temporal prefix
        o <- order(capped[[t]][[own]][rows])
        sort(rows[o[seq_len(cap_t)]])
      } else sort(rows[sample.int(length(rows), cap_t)])
    }), use.names = FALSE)
    keep <- sort(keep)
    dropped[[t]] <- nrow(capped[[t]]) - length(keep)
    capped[[t]] <- capped[[t]][keep, , drop = FALSE]
    rownames(capped[[t]]) <- NULL
  }
  list(tables = capped, dropped = dropped)
}

# Number of variable marginals a table contributes: one-way for every variable,
# plus all pairwise when a Chow-Liu tree is requested and there are >= 2 vars.
dp_n_var_marg <- function(nvars, dp) {
  if (dp$dependence == "tree" && nvars > 1L) nvars + nvars * (nvars - 1L) / 2L
  else nvars
}

# Enforce the per-table domain contract (identical to the flat engine, but
# reported with the table name) before any budget is spent.
dp_linked_domain_check <- function(cdata, vars, dp, table) {
  is_char <- vapply(vars, function(v) is.character(cdata[[v]]), logical(1))
  if (dp$domain %in% c("dp", "public") && any(is_char)) {
    stop(sprintf(paste0(
      "linked DP: table '%s' has character column(s) %s but domain = \"%s\" ",
      "needs a public category set. Convert them to factors (with full ",
      "`levels`), or use domain = \"data\" (not accounted)."),
      table, paste(vars[is_char], collapse = ", "), dp$domain), call. = FALSE)
  }
  numeric_vars <- vars[vapply(vars, function(v) is.numeric(cdata[[v]]), logical(1))]
  has_bound <- vapply(numeric_vars, function(v) !is.null(dp$bounds[[v]]), logical(1))
  missing_bound <- numeric_vars[!has_bound]
  if (dp$domain == "public" && length(missing_bound)) {
    stop(sprintf(paste0(
      "linked DP: domain = \"public\" needs `bounds` for every numeric variable; ",
      "table '%s' is missing %s. Supply public ranges or use domain = \"dp\"."),
      table, paste(missing_bound, collapse = ", ")), call. = FALSE)
  }
  missing_bound
}

# Draw a children-per-parent count vector (each in 0..cap) for `n_units` parents
# from the noisy, normalised count model `prob` (index i = P(count = i - 1)).
dp_draw_counts <- function(prob, n_units) {
  if (!n_units) return(integer(0))
  dp_sample_cat(n_units, prob) - 1L
}

# Top-level linked DP synthesiser. `hierarchy` is the parsed link hierarchy (see
# link_hierarchy()); `dp` a dp_control(); `m` the number of synthetic databases.
synth_linked_dp <- function(tables, hierarchy, dp, tuning, m, seed) {
  if (!is.null(seed)) set.seed(seed)
  if (identical(dp$select, "adaptive")) {
    stop(paste0("select = \"adaptive\" (AIM-style) is currently flat-table only ",
                "and is not supported for a linked DP release. Use the default ",
                "select = \"fixed\"."), call. = FALSE)
  }
  if (dp$degree > 1L) {
    stop(paste0("degree > 1 (a Bayesian network) is currently flat-table only ",
                "and is not supported for a linked DP release. Use degree = 1."),
         call. = FALSE)
  }
  if (dp$unit == "row") {
    stop(paste0("linked DP: unit = \"row\" is not supported; the privacy unit ",
                "of a linked release is the root entity. Use unit = \"person\"."),
         call. = FALSE)
  }
  order  <- hierarchy$order
  roots  <- order[vapply(order, function(t) is.na(hierarchy$parent[[t]]), logical(1))]
  if (!length(roots)) {
    stop("linked DP: no root table found in the hierarchy.", call. = FALSE)
  }

  caps <- dp_linked_caps(hierarchy, dp)

  # Which child tables are modelled as within-unit DP Markov trajectories (opt-in
  # via dp_control(longitudinal = ...)). Resolved before capping so those tables
  # are prefix-truncated in temporal order, not subsampled at random.
  use_long <- dp_resolve_longitudinal(dp, hierarchy, caps, tables)

  # Structural guard: a table may only nest along its key path plus (for a
  # longitudinal table) its own temporal index. A nesting index that is neither
  # part of the key nor the own index means a deeper within-table structure that
  # is not modelled under linked DP - the own-key-index Markov model covers the
  # common repeated-measures case.
  for (t in order) {
    extra <- setdiff(hierarchy$structure[[t]]$nested, hierarchy$keys[[t]])
    if (length(extra)) {
      stop(sprintf(paste0(
        "linked DP: table '%s' has a nesting index (%s) outside its key (%s). ",
        "A deeper within-table structure is not modelled under linked DP; drop ",
        "the extra index, or for repeated measures key the table by its own ",
        "temporal index and set dp_control(longitudinal = ...)."),
        t, paste(extra, collapse = "/"),
        paste(hierarchy$keys[[t]], collapse = ", ")), call. = FALSE)
    }
  }

  # Cap the whole hierarchy top-down, then work on the capped tables.
  capped_res <- dp_cap_hierarchy(tables, hierarchy, caps, longi = use_long)
  cdata <- capped_res$tables
  dropped <- capped_res$dropped

  # Per-table variable sets + domain contract; collect DP-estimated numeric vars.
  vars      <- stats::setNames(vector("list", length(order)), order)
  est_vars  <- stats::setNames(vector("list", length(order)), order)
  for (t in order) {
    vt <- dp_table_vars(hierarchy, t, cdata)
    if (is.na(hierarchy$parent[[t]]) && !length(vt)) {
      stop(sprintf(paste0(
        "linked DP: root table '%s' has no modellable variable (only keys). ",
        "The DP engine needs at least one non-key column to synthesise from."),
        t), call. = FALSE)
    }
    vars[[t]] <- vt
    mb <- dp_linked_domain_check(cdata[[t]], vt, dp, t)
    est_vars[[t]] <- if (dp$domain == "dp") mb else character(0)
  }

  # Which child tables condition on their parent (cross-table DP): opt-in, needs a
  # modellable immediate parent and >= 1 own variable. A non-longitudinal child
  # cross-conditions ALL its rows (use_cross); a longitudinal child instead
  # cross-conditions only its INITIAL-STATE model (use_long_cross), letting the
  # transition chain carry that parent dependence forward. Computed once so
  # composition and the fitted models agree exactly.
  cross_ok <- function(t) {
    p <- hierarchy$parent[[t]]
    isTRUE(dp$cross_table) && !is.na(p) &&
      length(vars[[t]]) > 0L && length(vars[[p]]) > 0L
  }
  use_cross      <- stats::setNames(logical(length(order)), order)
  use_long_cross <- stats::setNames(logical(length(order)), order)
  for (t in order) {
    ok <- cross_ok(t)
    use_cross[[t]]      <- ok && !use_long[[t]]
    use_long_cross[[t]] <- ok &&  use_long[[t]]
  }

  # Longitudinal children carry the flat engine's two within-unit transition
  # controls, applied per table: dp_control(baseline = ...) names subject-invariant
  # columns held constant within a unit (no transition histogram), and
  # dp_control(transition_order / transition_cross) deepen each time-varying
  # column's transition (own lags + lag-1 companions). `held_flag[[t]]` is a
  # logical over vars[[t]]; `ord` / `cross` are release-wide scalars. A higher
  # order needs a branching cap of at least ord + 1 so an order-deep within-unit
  # transition can be measured, so it is validated against each longitudinal
  # child's cap before any budget is spent.
  ord   <- if (is.null(dp$transition_order)) 1L else as.integer(dp$transition_order)
  cross <- if (is.null(dp$transition_cross)) 0L else as.integer(dp$transition_cross)
  tpar  <- if (is.null(dp$transition_parent)) 0L else as.integer(dp$transition_parent)
  base_cols <- if (is.null(dp$baseline)) character(0) else dp$baseline
  held_flag <- stats::setNames(vector("list", length(order)), order)
  for (t in order) {
    held_flag[[t]] <- vars[[t]] %in% base_cols
    if (use_long[[t]] && ord > caps$local[[t]] - 1L) {
      stop(sprintf(paste0(
        "linked DP longitudinal: transition_order (%d) must be <= branching cap - ",
        "1 (%d) for child table '%s' so order-%d within-unit transitions can be ",
        "measured under its cap. Raise its `max_rows_per_person` or lower ",
        "transition_order."), ord, caps$local[[t]] - 1L, t, ord), call. = FALSE)
    }
    # Parent-conditioned transitions reuse the parent-by-child joints that the
    # cross-conditioned initial state already measures, so they are budget-neutral
    # but need that model to exist: require cross_table on this longitudinal child.
    if (use_long[[t]] && tpar > 0L && !use_long_cross[[t]]) {
      stop(sprintf(paste0(
        "linked DP longitudinal: transition_parent (%d) for child table '%s' needs ",
        "its initial state cross-conditioned on a modellable parent - set ",
        "cross_table = TRUE (and ensure table '%s' has a modellable parent), or ",
        "transition_parent = 0. The parent-conditioned transitions reuse the ",
        "parent-by-child joints the cross-conditioned initial state measures, so ",
        "they add no budget."), tpar, t, t), call. = FALSE)
    }
  }

  # Per-table VARIABLE-model release: total L1, summed squared L2, and histogram
  # count. A longitudinal child releases initial-state marginals (at the parent
  # path cap) plus one transition matrix per variable (at parent path cap *
  # (branching cap - 1)); when it is also cross-conditioned the initial-state
  # count picks up the nC * nP parent-by-child joints (still at the parent path
  # cap). Other tables release plain / cross-conditioned variable marginals at
  # their own path cap. The children-per-parent count histogram is added
  # separately below (it doubles as the length model for longitudinal tables, so
  # no extra length histogram is charged).
  var_release <- function(t) {
    nC <- length(vars[[t]])
    if (use_long[[t]]) {
      pcp <- as.numeric(caps$path[[hierarchy$parent[[t]]]])
      lc  <- as.numeric(caps$local[[t]])
      # Initial-state marginals: plain n_init, or (combined cross-table) the
      # cross-conditioned count that adds nC * nP parent-by-child joints - all at
      # the same first-row (parent path-cap) sensitivity as the init marginals.
      # The init model spans ALL variables (baseline included), so ni is unchanged
      # by baseline; only the time-varying columns (nT) contribute a transition
      # histogram, at sensitivity pcp * (cap - order) - a higher order or any
      # cross-parent leaves ni untouched but lowers the transition sensitivity.
      ni  <- if (use_long_cross[[t]])
        dp_child_nvarmarg(nC, length(vars[[hierarchy$parent[[t]]]]), dp)
      else dp_longi_n_init(nC, dp)
      nT  <- sum(!held_flag[[t]])                          # time-varying columns
      ts  <- pcp * (lc - ord)                              # transition sensitivity
      list(l1 = ni * pcp + nT * ts,
           sq = ni * pcp^2 + nT * ts^2,
           nhist = as.integer(ni + nT))
    } else {
      pc  <- as.numeric(caps$path[[t]])
      nvm <- if (use_cross[[t]])
        dp_child_nvarmarg(nC, length(vars[[hierarchy$parent[[t]]]]), dp)
      else dp_n_var_marg(nC, dp)
      list(l1 = nvm * pc, sq = nvm * pc^2, nhist = as.integer(nvm))
    }
  }

  # ---- Exact composition over the whole release ---------------------------
  # total L1 (Laplace) and summed squared L2 (Gaussian zCDP) across every
  # histogram: root/child variable marginals (incl. parent-by-child joints when
  # cross_table, or initial-state + transitions when longitudinal) + child count
  # histograms.
  total_l1 <- 0
  sum_sq   <- 0
  n_hist   <- 0L
  for (t in order) {
    vr <- var_release(t)
    total_l1 <- total_l1 + vr$l1
    sum_sq   <- sum_sq   + vr$sq
    n_hist   <- n_hist + vr$nhist
    if (!is.na(hierarchy$parent[[t]])) {
      cs <- as.numeric(caps$path[[hierarchy$parent[[t]]]])   # count sensitivity
      total_l1 <- total_l1 + cs
      sum_sq   <- sum_sq   + cs^2
      n_hist   <- n_hist + 1L
    }
  }

  # ---- DP domain estimation slice (shared across all tables) --------------
  n_dom <- 2L * sum(vapply(est_vars, length, integer(1)))
  est_bounds  <- stats::setNames(vector("list", length(order)), order)
  marg_frac   <- 1
  domain_info <- list(mode = dp$domain, vars = character(0),
                      eps_per_query = NA_real_, frac = 0)
  if (n_dom > 0L) {
    eps_q <- dp_quantile_eps(dp, n_dom, dp$domain_frac)
    all_est <- character(0)
    for (t in order) {
      if (length(est_vars[[t]])) {
        pc <- caps$path[[t]]
        est_bounds[[t]] <- stats::setNames(
          lapply(est_vars[[t]],
                 function(v) dp_estimate_bounds(cdata[[t]][[v]], eps_q, pc)),
          est_vars[[t]])
        all_est <- c(all_est, paste0(t, ".", est_vars[[t]]))
      }
    }
    marg_frac <- 1 - dp$domain_frac
    domain_info <- list(mode = dp$domain, vars = all_est,
                        eps_per_query = eps_q, frac = dp$domain_frac)
  }

  calib <- dp_make_noise(dp, total_l1, sum_sq, budget_frac = marg_frac)

  # ---- Fit per-table models on the capped data ----------------------------
  # Topological order guarantees each parent's fit (domain + variable set) is
  # available before its children, which cross-table conditioning needs.
  fit <- stats::setNames(vector("list", length(order)), order)
  for (t in order) {
    vt  <- vars[[t]]
    dom <- if (length(vt)) dp_build_domain(cdata[[t]], vt, dp, est_bounds[[t]]) else list()
    nbins <- if (length(vt)) vapply(vt, function(v) dom[[v]]$nbin, integer(1)) else integer(0)
    var_model  <- NULL
    cross_model <- NULL
    long_model  <- NULL
    if (length(vt)) {
      if (use_long[[t]]) {
        parent_ctx <- NULL; parent_nbins <- NULL
        if (use_long_cross[[t]]) {                # cross-condition the init state
          p <- hierarchy$parent[[t]]
          pvars <- fit[[p]]$vars
          parent_dom <- fit[[p]]$dom
          parent_ctx <- dp_parent_ctx_codes(cdata, hierarchy, t, p, parent_dom,
                                             pvars)
          parent_nbins <- stats::setNames(
            vapply(pvars, function(u) parent_dom[[u]]$nbin, integer(1)), pvars)
        }
        long_model <- dp_fit_child_longitudinal(
          cdata[[t]], hierarchy$fk[[t]], hierarchy$own[[t]], dom, vt, nbins,
          dp, calib, parent_ctx, parent_nbins,
          held = held_flag[[t]], ord = ord, cross = cross, tran_parent = tpar)
      } else {
        codes <- stats::setNames(
          lapply(vt, function(v) dp_encode(dom[[v]], cdata[[t]][[v]])), vt)
        if (use_cross[[t]]) {
          p <- hierarchy$parent[[t]]
          pvars <- fit[[p]]$vars
          parent_dom <- fit[[p]]$dom
          pctx <- dp_parent_ctx_codes(cdata, hierarchy, t, p, parent_dom, pvars)
          parent_nbins <- stats::setNames(
            vapply(pvars, function(u) parent_dom[[u]]$nbin, integer(1)), pvars)
          cross_model <- dp_fit_child_cross(codes, nbins, pctx, parent_nbins,
                                            dp, calib)
        } else {
          var_model <- dp_fit_model(codes, nbins, dp, calib)
        }
      }
    }
    count_prob <- NULL
    if (!is.na(hierarchy$parent[[t]])) {
      p     <- hierarchy$parent[[t]]
      fk    <- hierarchy$fk[[t]]
      cap_t <- caps$local[[t]]
      pk <- key_string(cdata[[p]], fk)
      ck <- key_string(cdata[[t]], fk)
      cnt <- as.integer(table(factor(ck, levels = pk)))   # 0-inflated, aligned
      cnt <- pmin(cnt, cap_t)
      hist <- tabulate(cnt + 1L, cap_t + 1L)              # bin i = count (i - 1)
      count_prob <- dp_normalise(pmax(calib$add_noise(hist), 0))
    }
    fit[[t]] <- list(vars = vt, dom = dom, var_model = var_model,
                     cross_model = cross_model, long_model = long_model,
                     count_prob = count_prob)
  }

  # ---- Generation ---------------------------------------------------------
  n_root_target <- if (!is.null(tuning$k)) max(1L, as.integer(round(tuning$k)))
                   else NULL

  make_one <- function() {
    syn <- stats::setNames(vector("list", length(order)), order)
    # Cell codes behind each synthetic table (cols = that table's variables),
    # kept so a cross-table child can condition on its synthetic parent's codes.
    codes <- stats::setNames(vector("list", length(order)), order)
    for (t in order) {
      p  <- hierarchy$parent[[t]]
      f  <- fit[[t]]
      if (is.na(p)) {
        # Root: draw entities, assign a surrogate key 1..n.
        n_root <- n_root_target %||% max(1L, f$var_model$n_est)
        cmat <- dp_sample_codes(f$var_model, n_root)
        codes[[t]] <- cmat
        frame <- dp_decode_frame(cmat, f$dom)
        frame[[hierarchy$keys[[t]][1L]]] <- seq_len(n_root)
        syn[[t]] <- frame[names(tables[[t]])]
      } else {
        fk      <- hierarchy$fk[[t]]
        own_col <- hierarchy$own[[t]]
        parent_syn <- syn[[p]]
        counts <- dp_draw_counts(f$count_prob, nrow(parent_syn))
        keep   <- counts > 0L
        if (!any(keep)) {
          out <- tables[[t]][0, names(tables[[t]]), drop = FALSE]
          rownames(out) <- NULL
          syn[[t]] <- out
          codes[[t]] <- matrix(NA_integer_, 0L, length(f$vars),
                               dimnames = list(NULL, f$vars))
          next
        }
        rep_rows <- rep(which(keep), counts[keep])
        ppos <- sequence(counts[keep])          # 1..count within each unit
        skel <- parent_syn[rep_rows, fk, drop = FALSE]
        skel[[own_col]] <- dp_index_as_position(tables[[t]][[own_col]], ppos)
        if (length(f$vars)) {
          if (!is.null(f$long_model)) {
            lm <- f$long_model
            # When the initial state is cross-conditioned, seed each unit's first
            # row from the parent-conditioned model using the synthetic parent's
            # codes at that first row; otherwise the Markov engine draws it itself.
            first_codes <- NULL
            pctx_gen <- NULL
            if (isTRUE(lm$init_cross)) {
              im     <- lm$init_model
              frows  <- which(ppos == 1L)
              pctx_f <- stats::setNames(
                lapply(im$pvars, function(u) codes[[p]][rep_rows[frows], u]),
                im$pvars)
              first_codes <- dp_sample_child_codes(im, pctx_f, length(frows))
              # When the transitions also condition on parent attributes, carry the
              # synthetic parent's codes down to every child row (aligned to ppos).
              if (isTRUE(lm$tran_parent > 0L))
                pctx_gen <- stats::setNames(
                  lapply(im$pvars, function(u) codes[[p]][rep_rows, u]), im$pvars)
            }
            cmat <- if (isTRUE(lm$use_tensor))
              dp_markov_codes_tensor(lm$init_model, lm$tensors, ppos, f$vars,
                                     lm$held, lm$tv_vars, lm$order,
                                     first_codes = first_codes,
                                     parent_ctx = pctx_gen)
            else
              dp_markov_codes(lm$init_model, lm$tran, ppos, f$vars, lm$held,
                              first_codes = first_codes)
          } else if (!is.null(f$cross_model)) {
            pvars <- f$cross_model$pvars
            pctx <- stats::setNames(
              lapply(pvars, function(u) codes[[p]][rep_rows, u]), pvars)
            cmat <- dp_sample_child_codes(f$cross_model, pctx, length(rep_rows))
          } else {
            cmat <- dp_sample_codes(f$var_model, length(rep_rows))
          }
          codes[[t]] <- cmat
          vframe <- dp_decode_frame(cmat, f$dom)
          # Baseline (held) codes are equal within a unit, but numeric decoding
          # draws a fresh within-bin value per row; broadcast each unit's first-row
          # decoded value so a baseline column is byte-for-byte constant per unit.
          lm <- f$long_model
          if (!is.null(lm) && any(lm$held)) {
            g <- cumsum(ppos == 1L)                # unit id (rows contiguous)
            first_idx <- match(g, g)               # first-row index of each unit
            for (v in f$vars[lm$held]) vframe[[v]] <- vframe[[v]][first_idx]
          }
          for (v in f$vars) skel[[v]] <- vframe[[v]]
        } else {
          codes[[t]] <- matrix(NA_integer_, length(rep_rows), 0L)
        }
        rownames(skel) <- NULL
        syn[[t]] <- skel[names(tables[[t]])]
      }
    }
    syn
  }

  syns <- lapply(seq_len(m), function(i) make_one())
  syn <- if (m == 1L) syns[[1L]] else syns

  # ---- Accounting ---------------------------------------------------------
  tinfo <- lapply(order, function(t) {
    p  <- hierarchy$parent[[t]]
    lm <- fit[[t]]$long_model
    list(name = t,
         role = if (is.na(p)) "root" else "child",
         parent = if (is.na(p)) NA_character_ else p,
         local_cap = caps$local[[t]],
         path_cap = caps$path[[t]],
         n_var_marg = var_release(t)$nhist,
         cross = use_cross[[t]],
         cross_init = use_long_cross[[t]],
         longitudinal = use_long[[t]],
         baseline = if (!is.null(lm)) vars[[t]][lm$held] else character(0),
         tran_order = if (!is.null(lm)) lm$order else 1L,
         tran_cross = if (!is.null(lm)) lm$cross else 0L,
         tran_cross_parents = if (!is.null(lm)) lm$cross_parents else NULL,
         tran_parent = if (!is.null(lm)) lm$tran_parent else 0L,
         tran_parent_parents = if (!is.null(lm)) lm$parent_parents else NULL,
         count_sensitivity = if (is.na(p)) NA_integer_ else caps$path[[p]],
         rows_dropped = dropped[[t]])
  })
  linked <- list(tables = tinfo, n_tables = length(order))
  acct <- new_dp_accounting(
    dp, calib, cap = caps$local, n_marginals = n_hist,
    variables = unlist(lapply(order, function(t)
      if (length(vars[[t]])) paste0(t, ".", vars[[t]]) else character(0)),
      use.names = FALSE),
    dropped = sum(dropped), domain = domain_info, linked = linked)

  new_synth_linked_result(
    syn = syn, m = as.integer(m),
    n = vapply(tables, nrow, integer(1)), hierarchy = hierarchy,
    method = stats::setNames(rep(paste0("dp-", dp$mechanism), length(order)), order),
    privacy = acct, seed = seed, call = NULL
  )
}
