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

  # Release composition: 1 length histogram, the initial-state marginals, and one
  # transition matrix per variable. Sensitivities differ (transitions scale with
  # cap - 1), so we hand the exact totals to the shared calibrator.
  nV <- length(vars)
  tree <- dp$dependence == "tree" && nV > 1L
  n_init_marg <- nV + if (tree) nV * (nV - 1L) / 2L else 0L
  capf <- as.numeric(cap)
  total_l1 <- 1 + n_init_marg + nV * (capf - 1)
  sum_sq   <- 1 + n_init_marg + nV * (capf - 1)^2

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

  # Transition matrices over consecutive within-person pairs.
  if (length(cur_rows)) {
    prev_codes <- stats::setNames(lapply(vars, function(v) codes[[v]][prev_rows]),
                                  vars)
    cur_codes  <- stats::setNames(lapply(vars, function(v) codes[[v]][cur_rows]),
                                  vars)
    tran <- dp_fit_transitions(prev_codes, cur_codes, nbins, calib$add_noise)
  } else {
    # No repeated visits survived: transitions are unidentified but also unused
    # (every drawn length collapses to 1). Fall back to an identity-ish uniform.
    tran <- stats::setNames(
      lapply(vars, function(v) matrix(1 / nbins[[v]], nbins[[v]], nbins[[v]])),
      vars)
  }

  n_persons_est <- max(1L, round(sum(Lnoisy)))

  make_one <- function() {
    lens <- dp_draw_lengths(Lprob, n_persons_est, tuning$k)
    np <- length(lens)
    person <- rep(seq_len(np), lens)
    ppos <- sequence(lens)                    # 1..len within each person
    N <- length(person)

    cmat <- matrix(NA_integer_, N, nV, dimnames = list(NULL, vars))
    frows <- which(ppos == 1L)
    cmat[frows, ] <- dp_sample_codes(init_model, length(frows))
    tmax <- if (N) max(ppos) else 0L
    for (t in seq_len(tmax)[-1L]) {           # t = 2, 3, ...
      rows_t <- which(ppos == t)
      if (!length(rows_t)) next
      prv <- rows_t - 1L                      # previous position, same person
      for (vi in seq_len(nV)) {
        condv <- tran[[vi]]
        pv <- cmat[prv, vi]
        child <- integer(length(rows_t))
        for (a in unique(pv)) {
          sel <- which(pv == a)
          child[sel] <- dp_sample_cat(length(sel), condv[a, ])
        }
        cmat[rows_t, vi] <- child
      }
    }

    out <- dp_decode_frame(cmat, dom)
    out[[id]] <- person
    for (nk in st$nested) out[[nk]] <- dp_index_as_position(data[[nk]], ppos)
    out[names(data)]
  }

  syns <- lapply(seq_len(m), function(i) make_one())
  syn <- if (m == 1L) syns[[1L]] else syns

  n_hist <- 1L + as.integer(n_init_marg) + nV
  longitudinal <- list(n_init_marg = as.integer(n_init_marg),
                       n_transitions = nV, length_bins = as.integer(cap))
  acct <- new_dp_accounting(dp, calib, cap, n_hist, vars, bounded$dropped,
                            domain_info, longitudinal = longitudinal)
  new_synth_result(
    syn = syn, m = as.integer(m), n = nrow(data),
    method = stats::setNames(rep(paste0("dp-", dp$mechanism), length(vars)), vars),
    structure = structure, visit_sequence = vars, fixed = st$nested,
    subject = character(0), privacy = acct, seed = seed, call = NULL
  )
}
