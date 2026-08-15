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
# A longitudinally-modelled table is not simultaneously cross-conditioned.
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
  # modellable immediate parent and >= 1 own variable, and is not superseded by a
  # longitudinal model on the same table (which conditions on its own past, not
  # the parent). Computed once so composition and the fitted models agree exactly.
  use_cross <- stats::setNames(logical(length(order)), order)
  for (t in order) {
    p <- hierarchy$parent[[t]]
    use_cross[[t]] <- isTRUE(dp$cross_table) && !use_long[[t]] && !is.na(p) &&
      length(vars[[t]]) > 0L && length(vars[[p]]) > 0L
  }

  # Per-table VARIABLE-model release: total L1, summed squared L2, and histogram
  # count. A longitudinal child releases initial-state marginals (at the parent
  # path cap) plus one transition matrix per variable (at parent path cap *
  # (branching cap - 1)); other tables release plain / cross-conditioned variable
  # marginals at their own path cap. The children-per-parent count histogram is
  # added separately below (it doubles as the length model for longitudinal
  # tables, so no extra length histogram is charged).
  var_release <- function(t) {
    nC <- length(vars[[t]])
    if (use_long[[t]]) {
      pcp <- as.numeric(caps$path[[hierarchy$parent[[t]]]])
      lc  <- as.numeric(caps$local[[t]])
      ni  <- dp_longi_n_init(nC, dp)
      ts  <- pcp * (lc - 1)                                # transition sensitivity
      list(l1 = ni * pcp + nC * ts,
           sq = ni * pcp^2 + nC * ts^2,
           nhist = as.integer(ni + nC))
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
        long_model <- dp_fit_child_longitudinal(
          cdata[[t]], hierarchy$fk[[t]], hierarchy$own[[t]], dom, vt, nbins,
          dp, calib)
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
            cmat <- dp_markov_codes(lm$init_model, lm$tran, ppos, f$vars)
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
    p <- hierarchy$parent[[t]]
    list(name = t,
         role = if (is.na(p)) "root" else "child",
         parent = if (is.na(p)) NA_character_ else p,
         local_cap = caps$local[[t]],
         path_cap = caps$path[[t]],
         n_var_marg = var_release(t)$nhist,
         cross = use_cross[[t]],
         longitudinal = use_long[[t]],
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
