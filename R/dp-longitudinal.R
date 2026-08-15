# Internal: the Track B DP Markov synthesiser for nested / longitudinal tables.
# Not exported.
#
# When the `structure` formula declares a nesting index (e.g. ~ id / visit) the
# flat marginal release throws away within-unit temporal structure. This engine
# preserves it under differential privacy with a first-order Markov model, the
# DP analogue of Track A's initial-state + lag-1 transition model:
#
#   * length histogram   L(len)         — how many rows a person contributes;
#   * initial-state model (first rows)  — one-way (+ pairwise, for a Chow-Liu
#                                          tree) marginals over the t = 1 row;
#   * transition matrices P(v_t | v_{t-1}) per time-varying variable — the joint
#                                          of consecutive within-person values.
#
# Sensitivity (person level, at a public row cap `cap`): adding / removing one
# person moves the length histogram by 1 and each initial marginal by 1 (one
# first-row per person), and moves each transition histogram by at most cap - 1
# (a person of length <= cap contributes at most cap - 1 consecutive pairs). The
# whole release is one concatenated set of histograms, so the budget composes
# exactly (see `dp_make_noise()`): total L1 = 1 + n_init_marg + |V| * (cap - 1),
# and summed squared L2 = 1 + n_init_marg + |V| * (cap - 1)^2.
#
# Generation draws a length, an initial row, then steps each variable's
# transition matrix. Producing `m` datasets from one fitted model is
# post-processing and costs no extra budget.

# Prefix-truncate each person to their first `cap` rows in temporal order. Unlike
# the flat `dp_contribution_bound()` (which subsamples at random), the prefix is
# kept contiguous so consecutive-pair transitions are not corrupted. `data` must
# already be sorted so each unit's rows are contiguous and in temporal order.
dp_truncate_prefix <- function(data, id, cap) {
  gid <- data[[id]]
  pos <- stats::ave(seq_along(gid), gid, FUN = seq_along)
  keep <- pos <= cap
  list(data = data[keep, , drop = FALSE], dropped = sum(!keep))
}

# Fit the per-variable transition matrices P(v_t | v_{t-1}). `prev` / `cur` are
# named lists (per variable) of aligned code vectors over every consecutive
# within-person pair. Each matrix is row = previous cell, column = current cell,
# noised and row-normalised to a conditional distribution.
dp_fit_transitions <- function(prev, cur, nbins, add_noise) {
  vars <- names(nbins)
  out <- lapply(vars, function(v) {
    nb <- nbins[[v]]
    tt <- table(factor(prev[[v]], levels = seq_len(nb)),
                factor(cur[[v]],  levels = seq_len(nb)))
    tab <- matrix(as.integer(tt), nb, nb)
    noisy <- pmax(matrix(add_noise(as.vector(tab)), nb, nb), 0)
    t(apply(noisy, 1L, dp_normalise))        # rows -> P(current | previous)
  })
  stats::setNames(out, vars)
}

# Draw synthetic trajectory lengths (each in 1..cap). With a row budget `k` the
# loop accumulates whole persons until the budget is met; otherwise it draws
# `n_persons` lengths directly. `Lprob` is the (noisy, normalised) length model.
dp_draw_lengths <- function(Lprob, n_persons, k) {
  cap <- length(Lprob)
  draw1 <- function(m) sample.int(cap, m, replace = TRUE, prob = Lprob)
  if (is.null(k)) return(draw1(n_persons))
  lens <- integer(0); total <- 0L
  while (total < k) {
    l <- draw1(1L); lens <- c(lens, l); total <- total + l
  }
  lens
}

# Build a code matrix for trajectories laid out contiguously per unit, whose
# within-unit positions are `ppos` (each unit's rows are consecutive and its
# positions run 1..len in temporal order). The first row of each unit is drawn
# from `init_model`; every later row steps each variable through its transition
# matrix `tran[[vi]]` conditioned on the same unit's previous-position code. Used
# by both the flat DP longitudinal engine and a longitudinally-modelled linked
# child table. `vars` fixes the column order. `tran` is looked up by variable
# name, so it may hold only the time-varying variables. `held` is a logical
# vector over `vars` (default all FALSE): a held (subject-invariant / baseline)
# variable is drawn once in the initial row and then carried forward unchanged, so
# it stays exactly constant within the unit instead of stepping a transition.
dp_markov_codes <- function(init_model, tran, ppos, vars, held = NULL) {
  nV <- length(vars)
  if (is.null(held)) held <- logical(nV)
  N  <- length(ppos)
  cmat <- matrix(NA_integer_, N, nV, dimnames = list(NULL, vars))
  if (!N) return(cmat)
  frows <- which(ppos == 1L)
  cmat[frows, ] <- dp_sample_codes(init_model, length(frows))
  tmax <- max(ppos)
  for (tt in seq_len(tmax)[-1L]) {           # t = 2, 3, ...
    rows_t <- which(ppos == tt)
    if (!length(rows_t)) next
    prv <- rows_t - 1L                        # previous position, same unit
    for (vi in seq_len(nV)) {
      if (held[vi]) {                         # baseline: carry the value forward
        cmat[rows_t, vi] <- cmat[prv, vi]
        next
      }
      condv <- tran[[vars[vi]]]
      pv <- cmat[prv, vi]
      child <- integer(length(rows_t))
      for (a in unique(pv)) {
        sel <- which(pv == a)
        child[sel] <- dp_sample_cat(length(sel), condv[a, ])
      }
      cmat[rows_t, vi] <- child
    }
  }
  cmat
}

# --- Higher-order / cross-variable transition tensors -----------------------
# When `transition_order > 1` or `transition_cross > 0` the flat DP Markov engine
# conditions each time-varying variable on more than its own single previous
# value. The extra conditioning does not change the (eps, delta) budget (a
# transition tuple still lands in exactly one cell of one histogram); it trades
# budget-free structure for cell sparsity. See dp_control() for the semantics.

# Column-major linear index (1-based) into a set of aligned parent-code vectors.
# `codes` holds one integer vector per parent (values 1..nb[[i]]); `nb` the
# matching bin counts. The first parent varies fastest, matching the row order of
# a table built as matrix(<array with parent dims first, current last>, nparent, K).
dp_combo_index <- function(codes, nb) {
  if (!length(codes)) return(integer(0))
  mult <- cumprod(c(1, nb[-length(nb)]))
  idx <- 1L
  for (i in seq_along(codes)) idx <- idx + (codes[[i]] - 1L) * mult[i]
  as.integer(idx)
}

# For each time-varying variable, pick its `cross` most strongly associated other
# variables from the (noisy) pairwise mutual-information matrix `W` (dimnamed by
# variable name). Selection reads only already-released marginals, so it spends no
# budget and adds no leak. Returns a named list (per tv var) of parent var names.
dp_select_cross_parents <- function(W, tv_vars, all_vars, cross) {
  out <- stats::setNames(vector("list", length(tv_vars)), tv_vars)
  if (cross <= 0L || is.null(W)) {
    for (v in tv_vars) out[[v]] <- character(0)
    return(out)
  }
  for (v in tv_vars) {
    others <- setdiff(all_vars, v)
    if (!length(others)) { out[[v]] <- character(0); next }
    w <- W[v, others]
    out[[v]] <- others[utils::head(order(w, decreasing = TRUE), cross)]
  }
  out
}

# Fit one conditional transition tensor per time-varying variable. For variable v
# the tensor counts the tuples (v_t, v_{t-1..t-order}, u_{t-1} for each cross
# parent u) over the within-unit rows that have a full order-deep history
# (position > order), then adds calibrated noise. From that single noisy count
# array it precomputes, for every context depth d = 1..order, the conditional
# table P(v_t | v_{t-1..t-d}, cross) by marginalising out the deeper own-lags
# (pure post-processing), so early rows (position <= order) can still be generated
# without needing their own separate histograms. `codes` is the named list of full
# code vectors over all variables; `pos` the within-unit position of each row
# (rows already contiguous and temporal per unit).
dp_fit_transition_tensors <- function(codes, nbins, pos, order, cross_parents,
                                      tv_vars, add_noise) {
  cur <- which(pos > order)                  # rows with a full order-deep history
  fit_one <- function(v) {
    K  <- nbins[[v]]
    cp <- cross_parents[[v]]; if (is.null(cp)) cp <- character(0)
    cross_nb <- if (length(cp)) as.integer(nbins[cp]) else integer(0)
    dims <- c(K, rep(K, order), cross_nb)     # [current, own1..own_order, cross...]
    ncell <- prod(dims)
    if (length(cur)) {
      comps <- vector("list", length(dims))
      comps[[1L]] <- codes[[v]][cur]                       # current value
      for (j in seq_len(order)) comps[[1L + j]] <- codes[[v]][cur - j]
      for (k in seq_along(cp))
        comps[[1L + order + k]] <- codes[[cp[k]]][cur - 1L]
      counts <- tabulate(dp_combo_index(comps, dims), nbins = ncell)
    } else {
      counts <- integer(ncell)
    }
    A <- array(pmax(add_noise(counts), 0), dim = dims)
    cross_dims <- if (length(cp)) (order + 2L):(order + 1L + length(cp))
                  else integer(0)
    depth <- lapply(seq_len(order), function(d) {
      keep <- c(2L:(d + 1L), cross_dims, 1L)   # own1..d, cross, then current LAST
      red  <- if (length(keep) == length(dims)) aperm(A, keep)
              else apply(A, keep, sum)
      nparent <- prod(dims[keep[-length(keep)]])
      tab <- matrix(as.vector(red), nparent, K)
      t(apply(tab, 1L, dp_normalise))          # row = parent combo -> P(current|.)
    })
    list(depth_tables = depth, own = v, cross = cp,
         cross_nbins = cross_nb, K = as.integer(K), order = as.integer(order))
  }
  stats::setNames(lapply(tv_vars, fit_one), tv_vars)
}

# Higher-order / cross-variable analogue of dp_markov_codes(). Time-varying
# variables are stepped through their transition tensors: at position t the
# context is min(t - 1, order) own lags plus each variable's cross parents at lag
# 1, looked up in the precomputed depth-d conditional table. Held (baseline)
# variables are carried forward unchanged, exactly as in dp_markov_codes().
dp_markov_codes_tensor <- function(init_model, tensors, ppos, vars, held,
                                   tv_vars, order) {
  nV <- length(vars)
  N  <- length(ppos)
  cmat <- matrix(NA_integer_, N, nV, dimnames = list(NULL, vars))
  if (!N) return(cmat)
  if (is.null(held)) held <- logical(nV)
  vcol <- stats::setNames(seq_len(nV), vars)
  frows <- which(ppos == 1L)
  cmat[frows, ] <- dp_sample_codes(init_model, length(frows))
  tmax <- max(ppos)
  for (tt in seq_len(tmax)[-1L]) {
    rows_t <- which(ppos == tt)
    if (!length(rows_t)) next
    prv <- rows_t - 1L
    d <- min(tt - 1L, order)
    for (vi in seq_len(nV)) {
      if (held[vi]) { cmat[rows_t, vi] <- cmat[prv, vi]; next }
      te  <- tensors[[vars[vi]]]
      tab <- te$depth_tables[[d]]
      comps <- vector("list", d + length(te$cross))
      for (j in seq_len(d)) comps[[j]] <- cmat[rows_t - j, vi]
      for (k in seq_along(te$cross))
        comps[[d + k]] <- cmat[rows_t - 1L, vcol[[te$cross[k]]]]
      combo <- dp_combo_index(comps, c(rep(te$K, d), te$cross_nbins))
      child <- integer(length(rows_t))
      for (a in unique(combo)) {
        sel <- which(combo == a)
        child[sel] <- dp_sample_cat(length(sel), tab[a, ])
      }
      cmat[rows_t, vi] <- child
    }
  }
  cmat
}

# Regenerate a structural index column as the within-unit position, matching the
# original column's type (integer / numeric / factor). Factor levels are indexed
# by position and clamped to the available levels.
dp_index_as_position <- function(template, pos) {
  if (is.factor(template)) {
    lv <- levels(template)
    factor(lv[pmin(pos, length(lv))], levels = lv)
  } else if (is.integer(template)) {
    as.integer(pos)
  } else if (is.numeric(template)) {
    as.numeric(pos)
  } else {
    pos
  }
}

# Top-level DP Markov synthesiser. `st$nested` is the (regenerated) nesting
# index; the unit id becomes a fresh surrogate key. Every other column is a
# time-varying variable modelled by the initial-state + transition model.
synth_dp_longitudinal <- function(data, st, structure, dp, tuning, m, seed) {
  if (!is.null(seed)) set.seed(seed)
  id <- st$id
  vars <- setdiff(names(data), c(id, st$nested))
  if (!length(vars)) {
    stop(paste0("DP longitudinal synthesis needs at least one non-identifier, ",
                "non-structural column to model over time."), call. = FALSE)
  }

  cap <- dp_scalar_cap(dp)
  if (dp$unit == "row") cap <- 1L
  if (cap < 2L) {
    stop(sprintf(paste0(
      "DP longitudinal synthesis (structure = ~ %s) needs max_rows_per_person ",
      ">= 2 so within-person transitions can be measured; got %d. Set it from ",
      "public domain knowledge (the maximum visits per person), or drop the ",
      "nesting index for a flat DP release."),
      paste(st$vars, collapse = " / "), cap), call. = FALSE)
  }

  # Order so each unit's rows are contiguous and temporal, then prefix-truncate.
  rdat <- data[order_rows(data, st), , drop = FALSE]
  bounded <- dp_truncate_prefix(rdat, id, cap)
  cdata <- bounded$data

  # Domain contract (identical to the flat engine): rigorous modes require public
  # factor levels for categoricals and, in "public" mode, public numeric bounds.
  is_num  <- vapply(vars, function(v) is.numeric(cdata[[v]]), logical(1))
  is_char <- vapply(vars, function(v) is.character(cdata[[v]]), logical(1))
  numeric_vars <- vars[is_num]
  if (dp$domain %in% c("dp", "public") && any(is_char)) {
    stop(sprintf(paste0(
      "DP synthesis with domain = \"%s\" needs a public category set for %s: ",
      "convert the character column(s) to factors (with their full `levels`), ",
      "or use domain = \"data\" to read levels from the data (not accounted)."),
      dp$domain, paste(vars[is_char], collapse = ", ")), call. = FALSE)
  }
  has_bound <- vapply(numeric_vars,
                      function(v) !is.null(dp$bounds[[v]]), logical(1))
  missing_bound <- numeric_vars[!has_bound]
  if (dp$domain == "public" && length(missing_bound)) {
    stop(sprintf(paste0(
      "domain = \"public\" requires `bounds` for every numeric variable; ",
      "missing: %s. Supply public ranges, or use the default domain = \"dp\" ",
      "to estimate them privately."),
      paste(missing_bound, collapse = ", ")), call. = FALSE)
  }

  # Declared baseline (subject-invariant) columns are modelled once in the
  # initial state and then held constant within a unit, so they contribute NO
  # transition histogram. Only the time-varying variables get transitions.
  held_flag <- vars %in% dp$baseline
  held_vars <- vars[held_flag]
  tv_vars   <- vars[!held_flag]
  nT        <- length(tv_vars)                 # number of transition histograms

  # Transition model order (own lags) and cross-parent count. Higher order needs a
  # deeper history within the row cap; a person then contributes cap - order tuples
  # per variable, so the transition sensitivity drops with order.
  ord   <- if (is.null(dp$transition_order)) 1L else as.integer(dp$transition_order)
  cross <- if (is.null(dp$transition_cross)) 0L else as.integer(dp$transition_cross)
  if (ord > cap - 1L) {
    stop(sprintf(paste0(
      "transition_order (%d) must be <= max_rows_per_person - 1 (%d) so ",
      "order-%d within-person transitions can be measured under the row cap."),
      ord, cap - 1L, ord), call. = FALSE)
  }
  use_tensor <- (ord > 1L || cross > 0L) && nT > 0L

  # Release composition: 1 length histogram, the initial-state marginals (over
  # ALL variables, baseline included, so baseline correlations are kept), and one
  # transition histogram per time-varying variable. A transition histogram has
  # person-sensitivity cap - order regardless of how many columns condition it (a
  # tuple lands in one cell), so cross-conditioning is budget-free; only the order
  # moves the sensitivity. We hand the exact totals to the shared calibrator.
  nV <- length(vars)
  tree <- dp$dependence == "tree" && nV > 1L
  n_init_marg <- nV + if (tree) nV * (nV - 1L) / 2L else 0L
  capf <- as.numeric(cap)
  total_l1 <- 1 + n_init_marg + nT * (capf - ord)
  sum_sq   <- 1 + n_init_marg + nT * (capf - ord)^2

  # Under domain = "dp", privately estimate numeric bin edges (an accounted
  # `domain_frac` slice); the remaining budget pays for the histograms above.
  est_vars <- if (dp$domain == "dp") missing_bound else character(0)
  est_bounds <- NULL
  marg_frac <- 1
  domain_info <- list(mode = dp$domain, vars = character(0),
                      eps_per_query = NA_real_, frac = 0)
  if (length(est_vars) > 0L) {
    n_dom <- 2L * length(est_vars)
    eps_q <- dp_quantile_eps(dp, n_dom, dp$domain_frac)
    est_bounds <- stats::setNames(
      lapply(est_vars, function(v) dp_estimate_bounds(cdata[[v]], eps_q, cap)),
      est_vars)
    marg_frac <- 1 - dp$domain_frac
    domain_info <- list(mode = dp$domain, vars = est_vars,
                        eps_per_query = eps_q, frac = dp$domain_frac)
  }

  calib <- dp_make_noise(dp, total_l1, sum_sq, budget_frac = marg_frac)

  # Discretise every variable to the shared grid.
  dom <- dp_build_domain(cdata, vars, dp, est_bounds)
  nbins <- vapply(vars, function(v) dom[[v]]$nbin, integer(1))
  codes <- stats::setNames(
    lapply(vars, function(v) dp_encode(dom[[v]], cdata[[v]])), vars)

  # Split rows into first (initial state) and consecutive pairs (transitions).
  ids <- cdata[[id]]
  pos <- stats::ave(seq_along(ids), ids, FUN = seq_along)
  first_rows <- which(pos == 1L)
  cur_rows   <- which(pos > 1L)
  prev_rows  <- cur_rows - 1L                 # contiguous & ordered -> same person

  # Length model: one bin per length in 1..cap, person-sensitivity 1.
  plen <- as.integer(tabulate(pos[!duplicated(ids, fromLast = TRUE)], cap))
  # (tabulate over each person's final position = that person's trajectory length)
  Lnoisy <- pmax(calib$add_noise(plen), 0)
  Lprob <- dp_normalise(Lnoisy)

  # Initial-state model over first rows (reuses the flat marginal fitter).
  init_codes <- stats::setNames(lapply(vars, function(v) codes[[v]][first_rows]),
                                vars)
  init_model <- dp_fit_model(init_codes, nbins, dp, calib)

  # Transition model over consecutive within-person tuples, for the time-varying
  # variables only (baseline columns are held constant, so they need no transition
  # and drew no budget above). A first-order, own-lag-only model uses the simple
  # per-variable matrices; higher-order / cross-variable models use conditional
  # tensors. Both draw from `calib$add_noise` at the exact totals set above.
  tran <- NULL; tensors <- NULL; cross_parents <- NULL
  if (use_tensor) {
    cross_parents <- dp_select_cross_parents(init_model$pairwise_mi, tv_vars,
                                             vars, cross)
    tensors <- dp_fit_transition_tensors(codes, nbins, pos, ord, cross_parents,
                                         tv_vars, calib$add_noise)
  } else if (nT && length(cur_rows)) {
    prev_codes <- stats::setNames(
      lapply(tv_vars, function(v) codes[[v]][prev_rows]), tv_vars)
    cur_codes  <- stats::setNames(
      lapply(tv_vars, function(v) codes[[v]][cur_rows]), tv_vars)
    tran <- dp_fit_transitions(prev_codes, cur_codes, nbins[tv_vars],
                               calib$add_noise)
  } else {
    # No repeated visits survived (or nothing is time-varying): transitions are
    # unidentified but also unused. Fall back to an identity-ish uniform over the
    # time-varying variables.
    tran <- stats::setNames(
      lapply(tv_vars, function(v) matrix(1 / nbins[[v]], nbins[[v]], nbins[[v]])),
      tv_vars)
  }

  n_persons_est <- max(1L, round(sum(Lnoisy)))

  make_one <- function() {
    lens <- dp_draw_lengths(Lprob, n_persons_est, tuning$k)
    np <- length(lens)
    person <- rep(seq_len(np), lens)
    ppos <- sequence(lens)                    # 1..len within each person

    cmat <- if (use_tensor)
      dp_markov_codes_tensor(init_model, tensors, ppos, vars, held_flag,
                             tv_vars, ord)
    else
      dp_markov_codes(init_model, tran, ppos, vars, held_flag)

    out <- dp_decode_frame(cmat, dom)
    # Held codes are equal within a unit, but numeric decoding draws a fresh
    # within-bin value per row; broadcast each unit's first-row decoded value so a
    # baseline column is byte-for-byte constant within the unit.
    if (length(held_vars)) {
      first_idx <- match(person, person)       # first row index of each unit
      for (v in held_vars) out[[v]] <- out[[v]][first_idx]
    }
    out[[id]] <- person
    for (nk in st$nested) out[[nk]] <- dp_index_as_position(data[[nk]], ppos)
    out[names(data)]
  }

  syns <- lapply(seq_len(m), function(i) make_one())
  syn <- if (m == 1L) syns[[1L]] else syns

  n_hist <- 1L + as.integer(n_init_marg) + nT
  longitudinal <- list(n_init_marg = as.integer(n_init_marg),
                       n_transitions = nT, length_bins = as.integer(cap),
                       baseline = held_vars, order = ord, cross = cross,
                       cross_parents = if (use_tensor && cross > 0L) cross_parents
                                       else NULL)
  acct <- new_dp_accounting(dp, calib, cap, n_hist, vars, bounded$dropped,
                            domain_info, longitudinal = longitudinal)
  new_synth_result(
    syn = syn, m = as.integer(m), n = nrow(data),
    method = stats::setNames(rep(paste0("dp-", dp$mechanism), length(vars)), vars),
    structure = structure, visit_sequence = vars, fixed = st$nested,
    subject = character(0), privacy = acct, seed = seed, call = NULL
  )
}
