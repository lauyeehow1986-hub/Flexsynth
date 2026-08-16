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

# Integer joint-count matrix of two code vectors, `a` as rows (1..na), `b` as
# columns (1..nb); NA-paired rows dropped. Shared by the naive and the
# budget-efficient tree fitters.
dp_joint_counts <- function(a, b, na, nb) {
  ok <- !is.na(a) & !is.na(b)
  tt <- table(factor(a[ok], levels = seq_len(na)),
              factor(b[ok], levels = seq_len(nb)))
  matrix(as.integer(tt), na, nb)
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
#
# `calib` calibrates the noise for the released parameters (one-way marginals and,
# for a tree, its conditional edges). `calib_struct` is optional: when supplied
# (budget-efficient structure learning) the Chow-Liu tree is *selected* from a
# separate, cheaper pass of all pairwise joints measured under `calib_struct`,
# and then only the chosen tree's edges are re-measured under `calib` — so the
# bulk of the budget lands on the parameters that survive into the model instead
# of on all C(d, 2) pairs. When `calib_struct` is NULL (the default) the pairwise
# joints are measured once under `calib` and reused for both structure and
# parameters, exactly as before.
dp_fit_model <- function(codes, nbins, dp, calib, calib_struct = NULL) {
  d <- length(codes)
  vars <- names(codes)

  # One-way noisy counts (clipped nonneg), measured under the parameter budget.
  c1 <- lapply(seq_len(d), function(i)
    pmax(calib$add_noise(tabulate(codes[[i]], nbins[i])), 0))
  names(c1) <- vars

  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))

  if (dp$dependence == "independent" || d <= 1L) {
    return(list(kind = "independent", vars = vars, nbins = nbins,
                marginals = lapply(c1, dp_normalise), n_est = n_est))
  }

  if (is.null(calib_struct)) {
    # Naive: one pass of all pairwise joints under `calib`, reused for the MST
    # weights and the conditional parameters.
    joints <- matrix(list(), d, d)
    W <- matrix(0, d, d)
    for (i in seq_len(d - 1L)) for (j in (i + 1L):d) {
      tab <- dp_joint_counts(codes[[i]], codes[[j]], nbins[i], nbins[j])
      noisy <- pmax(matrix(calib$add_noise(as.vector(tab)),
                           nbins[i], nbins[j]), 0)
      joints[[i, j]] <- noisy
      W[i, j] <- W[j, i] <- dp_mutual_information(noisy)
    }
    edges <- dp_max_spanning_tree(W)
    cond <- lapply(edges, function(e) {
      p <- e[1L]; ch <- e[2L]
      joint <- if (p < ch) joints[[p, ch]] else t(joints[[ch, p]])  # rows=parent
      t(apply(joint, 1L, dp_normalise))               # P(child | parent)
    })
  } else {
    # Budget-efficient: a cheap all-pairs scan under `calib_struct` picks the
    # tree; only its edges are then re-measured under `calib`.
    W <- matrix(0, d, d)
    for (i in seq_len(d - 1L)) for (j in (i + 1L):d) {
      tab <- dp_joint_counts(codes[[i]], codes[[j]], nbins[i], nbins[j])
      rough <- pmax(matrix(calib_struct$add_noise(as.vector(tab)),
                           nbins[i], nbins[j]), 0)
      W[i, j] <- W[j, i] <- dp_mutual_information(rough)
    }
    edges <- dp_max_spanning_tree(W)
    cond <- lapply(edges, function(e) {
      p <- e[1L]; ch <- e[2L]                          # parent as rows
      tab <- dp_joint_counts(codes[[p]], codes[[ch]], nbins[p], nbins[ch])
      noisy <- pmax(matrix(calib$add_noise(as.vector(tab)),
                           nbins[p], nbins[ch]), 0)
      t(apply(noisy, 1L, dp_normalise))                # P(child | parent)
    })
  }

  # Root marginal = the one-way of node 1 (the Prim root). `pairwise_mi` is the
  # full (noisy) mutual-information matrix over the modelled variables; it is
  # post-processing of already-measured marginals, so it is free to reuse — the DP
  # longitudinal engine consults it to pick cross-variable transition parents.
  dimnames(W) <- list(vars, vars)
  list(kind = "tree", vars = vars, nbins = nbins, edges = edges,
       cond = cond, root = 1L, root_marginal = dp_normalise(c1[[1L]]),
       marginals = lapply(c1, dp_normalise), n_est = n_est, pairwise_mi = W)
}

# Fit a flat tree release with Private-PGM reconciliation (MST). Measures the
# same marginals as the local single-pass tree fitter -- the d one-way counts and
# all C(d, 2) pairwise joints under `calib`, the pairs giving the Chow-Liu MST
# weights -- so the (eps, delta) is identical. It then keeps the chosen tree's
# d - 1 edges (as a treewidth-1 junction tree, root = node 1) plus the one-ways
# and reconciles that set into one consistent distribution (dp_pgm_reconcile),
# rather than using each edge's raw noisy 2-way locally. Returns a junction-kind
# model so the existing junction sampler draws it unchanged.
dp_fit_model_pgm <- function(codes, nbins, dp, calib) {
  d <- length(codes)
  vars <- names(codes)
  c1 <- lapply(seq_len(d), function(i)
    pmax(calib$add_noise(tabulate(codes[[i]], nbins[i])), 0))
  names(c1) <- vars
  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))

  joints <- matrix(list(), d, d)
  W <- matrix(0, d, d)
  for (i in seq_len(d - 1L)) for (j in (i + 1L):d) {
    tab <- dp_joint_counts(codes[[i]], codes[[j]], nbins[i], nbins[j])
    noisy <- pmax(matrix(calib$add_noise(as.vector(tab)),
                         nbins[i], nbins[j]), 0)
    joints[[i, j]] <- noisy
    W[i, j] <- W[j, i] <- dp_mutual_information(noisy)
  }
  edges <- dp_max_spanning_tree(W)                 # directed (parent, child)

  cliques <- vector("list", d)
  cliques[[1L]] <- list(vars = 1L, sep = integer(0), new = 1L)
  meas <- lapply(seq_len(d), function(i) list(vars = i, y = c1[[i]]))
  for (k in seq_along(edges)) {
    e <- edges[[k]]; p <- e[1L]; ch <- e[2L]
    cliques[[k + 1L]] <- list(vars = sort(c(p, ch)), sep = p, new = ch)
    # The stored joint has the smaller variable index as rows (fastest), which is
    # exactly the sorted-variable cell order dp_pgm_reconcile() expects.
    J <- joints[[min(p, ch), max(p, ch)]]
    meas[[length(meas) + 1L]] <- list(vars = sort(c(p, ch)), y = as.vector(J))
  }
  rec <- dp_pgm_reconcile(cliques, meas, nbins, vars, n_est)
  list(kind = "junction", vars = vars, nbins = nbins, cliques = rec$cliques,
       marginals = lapply(c1, dp_normalise), n_est = n_est, estimator = "pgm")
}

# Draw one synthetic dataset of `n` rows (matrix of cell codes, n x d).
dp_sample_codes <- function(model, n) {
  if (identical(model$kind, "pgm")) return(dp_sample_codes_pgm(model, n))
  if (model$kind %in% c("junction", "bayes")) return(dp_sample_codes_junction(model, n))
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
  if (identical(dp$select, "adaptive") && length(st$nested) > 0L) {
    stop(paste0("select = \"adaptive\" (AIM-style) is currently flat-table only; ",
                "this structure is longitudinal. Use the default select = \"fixed\" ",
                "for a longitudinal DP release."), call. = FALSE)
  }
  if (dp$degree > 1L && length(st$nested) > 0L) {
    stop(paste0("degree > 1 (a Bayesian network) is currently flat-table only; ",
                "this structure is longitudinal. Use degree = 1 for a ",
                "longitudinal DP release."), call. = FALSE)
  }
  if (identical(dp$estimator, "pgm") && length(st$nested) > 0L) {
    stop(paste0("estimator = \"pgm\" (Private-PGM reconciliation) is currently ",
                "flat-table only; this structure is longitudinal. Use the ",
                "default estimator = \"local\" for a longitudinal DP release."),
         call. = FALSE)
  }
  if (identical(dp$select, "aim") && length(st$nested) > 0L) {
    stop(paste0("select = \"aim\" (Full AIM) is currently flat-table only; this ",
                "structure is longitudinal. Use the default select = \"fixed\" ",
                "for a longitudinal DP release."), call. = FALSE)
  }
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

  d <- length(vars)
  # AIM-style adaptive selection: measure d one-ways, then greedily select
  # `d - treewidth` cliques by the exponential mechanism (see dp-adaptive.R).
  # Needs at least two variables; the treewidth is capped to d - 1.
  adaptive <- identical(dp$select, "adaptive") && d >= 2L
  w_eff <- if (adaptive) min(dp$treewidth, d - 1L) else NA_integer_
  n_cliques <- if (adaptive) d - w_eff else 0L
  # Full AIM: measure d one-ways, then a fixed number of exponential-mechanism
  # rounds over *loopy* pairwise marginals (bounded by the treewidth after
  # triangulation), reconciled and sampled with Private-PGM (see dp-aim.R). The
  # round count is data-independent: min(all pairs, treewidth * (d - 1)).
  aim <- identical(dp$select, "aim") && d >= 2L
  w_aim <- if (aim) min(dp$treewidth, d - 1L) else NA_integer_
  n_rounds_aim <- if (aim)
    as.integer(min(d * (d - 1L) / 2L, w_aim * (d - 1L))) else 0L
  # GreedyBayes: a degree-k Bayesian network (each node conditions on up to
  # `degree` already-generated predecessors, chosen by the exponential mechanism).
  # Generalises the Chow-Liu tree (degree 1); needs at least two variables and the
  # degree is capped to d - 1 (see dp-adaptive.R).
  bayes <- !adaptive && !aim && dp$degree > 1L && dp$dependence == "tree" && d >= 2L
  deg_eff <- if (bayes) min(dp$degree, d - 1L) else NA_integer_
  n_pairs <- if (!adaptive && !aim && !bayes && dp$dependence == "tree" && d > 1L)
    d * (d - 1L) / 2L else 0L
  # Budget-efficient structure learning: a separate cheap pairwise scan selects
  # the tree, then only its d-1 edges (plus d one-ways) are re-measured. Trivial
  # for < 3 variables, so it is inert there (and for the independent model).
  use_learn <- !adaptive && !aim && !bayes && !is.null(dp$structure_frac) &&
    dp$dependence == "tree" && d >= 3L
  # Private-PGM reconciliation (MST) of the flat tree's measured marginals. Same
  # measured set (d one-ways + all pairs) and budget as the local tree, so it is
  # a drop-in for the plain tree fitter; only the model built from them changes.
  pgm <- identical(dp$estimator, "pgm")
  pgm_tree <- pgm && !adaptive && !aim && !bayes && dp$dependence == "tree" && d >= 2L
  # A degree-k network measures d one-ways + (d - 1) family joints = 2d - 1.
  # Full AIM measures d one-ways + n_rounds_aim selected pairwise marginals.
  n_marginals <- if (adaptive) d + n_cliques
    else if (aim) d + n_rounds_aim
    else if (bayes) 2L * d - 1L
    else if (use_learn) n_pairs + (2L * d - 1L) else d + n_pairs

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

  # Split the marginal budget between the (cheap) structure scan and the
  # (concentrated) parameter re-measurement, or spend it all in one pass. Both
  # slices compose exactly into the same (eps, delta) as the domain slice: pure
  # eps adds; zCDP rho adds.
  learn_info <- NULL
  adapt_info <- NULL
  bayes_info <- NULL
  aim_info <- NULL
  sel_eps_bayes <- NULL
  sel_eps_aim <- NULL
  pool_meas <- NULL
  pool_sel <- NULL
  if (aim) {
    # Full AIM: the marginal budget splits into a selection slice (the
    # exponential-mechanism rounds over loopy pairs) and a measurement slice
    # (the d one-ways + n_rounds_aim selected pairs); both compose exactly into
    # the same (eps, delta) as the domain slice (zCDP rho adds; pure eps adds).
    sf <- dp$select_frac
    calib_struct <- NULL
    calib <- dp_calibrate(dp, n_marginals, cap,
                          budget_frac = marg_frac * (1 - sf))
    sel_eps_aim <- dp_select_eps(dp, n_rounds_aim, marg_frac * sf)
    aim_info <- list(treewidth = w_aim, n_rounds = n_rounds_aim,
                     select_frac = sf, select_eps = sel_eps_aim,
                     meas_noise = if (calib$mechanism == "laplace")
                       calib$scale else calib$sigma)
  } else if (adaptive) {
    # Marginal budget splits into a selection slice (the exponential-mechanism
    # rounds) and a measurement slice (the d one-ways + d - w cliques); both
    # compose exactly into the same (eps, delta) as the domain slice.
    sel_frac <- dp$select_frac
    calib_struct <- NULL
    calib <- dp_calibrate(dp, n_marginals, cap,
                          budget_frac = marg_frac * (1 - sel_frac))
    if (isTRUE(dp$anneal)) {
      # Annealed schedule: hand the measurement and selection pools straight to
      # the engine, which spends them over a data-adaptive number of rounds.
      total <- if (dp$mechanism == "gaussian")
        zcdp_rho_for(dp$epsilon, dp$delta) else dp$epsilon
      pool_meas <- marg_frac * (1 - sel_frac) * total
      pool_sel  <- marg_frac * sel_frac * total
    } else {
      sel_eps <- dp_select_eps(dp, n_cliques, marg_frac * sel_frac)
      adapt_info <- list(treewidth = w_eff, n_cliques = n_cliques,
                         select_frac = sel_frac, select_eps = sel_eps,
                         meas_noise = if (calib$mechanism == "laplace")
                           calib$scale else calib$sigma)
    }
  } else if (bayes) {
    # Marginal budget splits into a selection slice (the d - 1 exponential-
    # mechanism parent-set picks) and a measurement slice (the 2d - 1 measured
    # marginals); both compose exactly into the same (eps, delta) as the domain
    # slice (pure eps adds; zCDP rho adds).
    sf <- dp$select_frac
    calib_struct <- NULL
    calib <- dp_calibrate(dp, n_marginals, cap,
                          budget_frac = marg_frac * (1 - sf))
    sel_eps_bayes <- dp_select_eps(dp, d - 1L, marg_frac * sf)
    bayes_info <- list(degree = deg_eff, n_nodes = d, n_families = d - 1L,
                       n_select = d - 1L, select_frac = sf,
                       select_eps = sel_eps_bayes,
                       meas_noise = if (calib$mechanism == "laplace")
                         calib$scale else calib$sigma)
    if (dp$mechanism == "gaussian") {
      rho_total <- zcdp_rho_for(dp$epsilon, dp$delta)
      bayes_info$rho_meas <- marg_frac * (1 - sf) * rho_total
      bayes_info$rho_sel  <- marg_frac * sf * rho_total
    } else {
      bayes_info$eps_meas <- marg_frac * (1 - sf) * dp$epsilon
      bayes_info$eps_sel  <- marg_frac * sf * dp$epsilon
    }
  } else if (use_learn) {
    n_struct <- n_pairs
    n_param  <- 2L * d - 1L
    calib_struct <- dp_calibrate(dp, n_struct, cap,
                                 budget_frac = marg_frac * dp$structure_frac)
    calib <- dp_calibrate(dp, n_param, cap,
                          budget_frac = marg_frac * (1 - dp$structure_frac))
    learn_info <- list(frac = dp$structure_frac, n_struct = n_struct,
                       n_param = n_param,
                       struct_noise = if (calib_struct$mechanism == "laplace")
                         calib_struct$scale else calib_struct$sigma)
  } else {
    calib_struct <- NULL
    calib <- dp_calibrate(dp, n_marginals, cap, budget_frac = marg_frac)
  }

  dom <- dp_build_domain(cdata, vars, dp, est_bounds)
  nbins <- vapply(vars, function(v) dom[[v]]$nbin, integer(1))
  codes <- stats::setNames(
    lapply(vars, function(v) dp_encode(dom[[v]], cdata[[v]])), vars)

  model <- if (aim) {
    dp_fit_model_aim(codes, nbins, dp, calib, w_aim, sel_eps_aim, cap)
  } else if (adaptive) {
    if (isTRUE(dp$anneal)) {
      fit <- dp_fit_model_adaptive_anneal(codes, nbins, dp, w_eff, cap,
                                          pool_meas, pool_sel)
      adapt_info <- fit$anneal
      fit$model
    } else {
      dp_fit_model_adaptive(codes, nbins, dp, calib, w_eff,
                            adapt_info$select_eps, cap, reconcile = pgm)
    }
  } else if (bayes) {
    dp_fit_model_bayes(codes, nbins, dp, calib, dp$degree, sel_eps_bayes, cap)
  } else if (pgm_tree) {
    dp_fit_model_pgm(codes, nbins, dp, calib)
  } else dp_fit_model(codes, nbins, dp, calib, calib_struct)
  # The annealed path measures a data-adaptive number of histograms (d one-ways +
  # the realised clique rounds); report that actual count.
  n_marg_report <- if (adaptive && isTRUE(dp$anneal))
    d + adapt_info$n_rounds else n_marginals
  n_syn <- if (!is.null(tuning$k)) max(1L, as.integer(round(tuning$k)))
           else model$n_est

  make_one <- function() {
    out <- dp_decode_frame(dp_sample_codes(model, n_syn), dom)
    out[[id]] <- seq_len(n_syn)
    out[names(data)]
  }
  syns <- lapply(seq_len(m), function(i) make_one())
  syn <- if (m == 1L) syns[[1L]] else syns

  acct <- new_dp_accounting(dp, calib, cap, n_marg_report, vars, bounded$dropped,
                            domain_info, learn = learn_info, adaptive = adapt_info,
                            bayes = bayes_info, aim = aim_info,
                            estimator = if (pgm && !aim) "pgm" else "local")
  new_synth_result(
    syn = syn, m = as.integer(m), n = nrow(data),
    method = stats::setNames(rep(paste0("dp-", dp$mechanism), length(vars)), vars),
    structure = structure, visit_sequence = vars, fixed = character(0),
    subject = character(0), privacy = acct, seed = seed, call = NULL
  )
}
