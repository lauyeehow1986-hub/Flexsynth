# Internal: the Track B differentially private marginal synthesiser.
# Not exported.
#
# Pipeline: bound each person's contribution -> discretise every variable to a
# public grid -> measure low-order marginals under the calibrated noise
# mechanism (with correct budget composition) -> fit a generative model over the
# discrete cells (independent one-way, or a Chow-Liu tree learned from the same
# noisy marginals) -> draw synthetic cells and decode to typed values. Drawing
# the synthetic records is post-processing of the noisy marginals, so producing
# `m` datasets from one fitted model costs no extra privacy budget.

# Normalise a nonneg weight vector to a probability vector; all-zero -> uniform.
dp_normalise <- function(x) {
  x <- pmax(x, 0)
  s <- sum(x)
  if (s <= 0) rep(1 / length(x), length(x)) else x / s
}

# Sample n categorical draws (indices 1..length(prob)) from a weight vector.
dp_sample_cat <- function(n, prob) {
  prob <- dp_normalise(prob)
  sample.int(length(prob), n, replace = TRUE, prob = prob)
}

# Pairwise mutual information from a noisy joint count matrix.
dp_mutual_information <- function(joint) {
  p <- joint / sum(joint)
  pr <- rowSums(p); pc <- colSums(p)
  mi <- 0
  nz <- which(p > 0, arr.ind = TRUE)
  for (k in seq_len(nrow(nz))) {
    a <- nz[k, 1L]; b <- nz[k, 2L]
    denom <- pr[a] * pc[b]
    if (denom > 0) mi <- mi + p[a, b] * log(p[a, b] / denom)
  }
  max(mi, 0)
}

# Maximum-weight spanning tree by Prim's algorithm, rooted at node 1. Returns a
# list of directed c(parent, child) edges in an order where each edge's parent is
# already reachable from the root (a valid sampling order).
dp_max_spanning_tree <- function(W) {
  d <- nrow(W)
  if (d <= 1L) return(list())
  in_tree <- logical(d); in_tree[1L] <- TRUE
  best <- W[1L, ]; best_from <- rep(1L, d)
  edges <- vector("list", d - 1L)
  for (step in seq_len(d - 1L)) {
    cand <- which(!in_tree)
    j <- cand[which.max(best[cand])]
    in_tree[j] <- TRUE
    edges[[step]] <- c(best_from[j], j)
    upd <- !in_tree & (W[j, ] > best)
    best[upd] <- W[j, upd]
    best_from[upd] <- j
  }
  edges
}

# Measure the noisy marginals and fit the generative model.
dp_fit_model <- function(codes, nbins, dp, calib) {
  d <- length(codes)
  vars <- names(codes)

  # One-way noisy counts (clipped nonneg).
  c1 <- lapply(seq_len(d), function(i)
    pmax(calib$add_noise(tabulate(codes[[i]], nbins[i])), 0))
  names(c1) <- vars

  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))

  if (dp$dependence == "independent" || d <= 1L) {
    return(list(kind = "independent", vars = vars, nbins = nbins,
                marginals = lapply(c1, dp_normalise), n_est = n_est))
  }

  # All pairwise noisy joints (i < j), reused for both structure and parameters.
  joints <- matrix(list(), d, d)
  W <- matrix(0, d, d)
  for (i in seq_len(d - 1L)) for (j in (i + 1L):d) {
    tab <- matrix(0L, nbins[i], nbins[j])
    ok <- !is.na(codes[[i]]) & !is.na(codes[[j]])
    idx <- cbind(codes[[i]][ok], codes[[j]][ok])
    tt <- table(factor(idx[, 1L], levels = seq_len(nbins[i])),
                factor(idx[, 2L], levels = seq_len(nbins[j])))
    tab[] <- as.integer(tt)
    noisy <- pmax(matrix(calib$add_noise(as.vector(tab)), nbins[i], nbins[j]), 0)
    joints[[i, j]] <- noisy
    W[i, j] <- W[j, i] <- dp_mutual_information(noisy)
  }

  edges <- dp_max_spanning_tree(W)
  # Root marginal = the one-way of node 1 (the Prim root).
  cond <- lapply(edges, function(e) {
    p <- e[1L]; ch <- e[2L]
    joint <- if (p < ch) joints[[p, ch]] else t(joints[[ch, p]])  # rows=parent
    # Row-normalise to P(child | parent); empty rows -> uniform.
    t(apply(joint, 1L, dp_normalise))
  })

  list(kind = "tree", vars = vars, nbins = nbins, edges = edges,
       cond = cond, root = 1L, root_marginal = dp_normalise(c1[[1L]]),
       marginals = lapply(c1, dp_normalise), n_est = n_est)
}

# Draw one synthetic dataset of `n` rows (matrix of cell codes, n x d).
dp_sample_codes <- function(model, n) {
  d <- length(model$vars)
  out <- matrix(NA_integer_, n, d, dimnames = list(NULL, model$vars))
  if (model$kind == "independent" || d <= 1L) {
    for (i in seq_len(d))
      out[, i] <- dp_sample_cat(n, model$marginals[[i]])
    return(out)
  }
  out[, model$root] <- dp_sample_cat(n, model$root_marginal)
  for (k in seq_along(model$edges)) {
    e <- model$edges[[k]]; p <- e[1L]; ch <- e[2L]
    cond <- model$cond[[k]]
    pc <- out[, p]
    child <- integer(n)
    for (a in unique(pc)) {
      rows <- which(pc == a)
      child[rows] <- dp_sample_cat(length(rows), cond[a, ])
    }
    out[, ch] <- child
  }
  out
}

# Decode a code matrix into a typed data.frame in the domain's variable order.
dp_decode_frame <- function(codes, dom) {
  vars <- colnames(codes)
  cols <- lapply(vars, function(v) dp_decode(dom[[v]], codes[, v]))
  as.data.frame(stats::setNames(cols, vars), stringsAsFactors = FALSE)
}

# Top-level DP synthesiser. `st` is the parsed structure; the unit id becomes a
# surrogate key (1..n) and every other column is modelled as a discrete variable.
# When the structure declares a nesting index (`~ id / visit`) the DP Markov
# engine preserves within-unit temporal structure (see dp-longitudinal.R);
# otherwise this produces a flat marginal release.
synth_dp <- function(data, st, structure, dp, tuning, m, seed) {
  if (length(st$nested) > 0L) {
    return(synth_dp_longitudinal(data, st, structure, dp, tuning, m, seed))
  }
  if (!is.null(seed)) set.seed(seed)
  id <- st$id
  vars <- setdiff(names(data), id)
  if (!length(vars)) {
    stop("DP synthesis needs at least one non-identifier column.", call. = FALSE)
  }

  cap <- dp_scalar_cap(dp)
  if (dp$unit == "row") cap <- 1L
  bounded <- dp_contribution_bound(data, id, cap)
  cdata <- bounded$data

  # Classify columns and enforce the domain contract before spending budget.
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

  n_marginals <- length(vars) +
    if (dp$dependence == "tree" && length(vars) > 1L)
      length(vars) * (length(vars) - 1L) / 2L else 0L

  # Under domain = "dp", privately estimate bin edges for numeric variables that
  # have no public bounds, spending an accounted `domain_frac` slice of budget.
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

  calib <- dp_calibrate(dp, n_marginals, cap, budget_frac = marg_frac)

  dom <- dp_build_domain(cdata, vars, dp, est_bounds)
  nbins <- vapply(vars, function(v) dom[[v]]$nbin, integer(1))
  codes <- stats::setNames(
    lapply(vars, function(v) dp_encode(dom[[v]], cdata[[v]])), vars)

  model <- dp_fit_model(codes, nbins, dp, calib)
  n_syn <- if (!is.null(tuning$k)) max(1L, as.integer(round(tuning$k)))
           else model$n_est

  make_one <- function() {
    out <- dp_decode_frame(dp_sample_codes(model, n_syn), dom)
    out[[id]] <- seq_len(n_syn)
    out[names(data)]
  }
  syns <- lapply(seq_len(m), function(i) make_one())
  syn <- if (m == 1L) syns[[1L]] else syns

  acct <- new_dp_accounting(dp, calib, cap, n_marginals, vars, bounded$dropped,
                            domain_info)
  new_synth_result(
    syn = syn, m = as.integer(m), n = nrow(data),
    method = stats::setNames(rep(paste0("dp-", dp$mechanism), length(vars)), vars),
    structure = structure, visit_sequence = vars, fixed = character(0),
    subject = character(0), privacy = acct, seed = seed, call = NULL
  )
}
