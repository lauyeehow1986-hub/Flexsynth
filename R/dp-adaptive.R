# Internal: AIM-style adaptive marginal selection for the flat Track B engine.
# Not exported.
#
# Instead of measuring a predetermined set of marginals (all one-ways, plus all
# pairwise marginals for a Chow-Liu tree), this selects the model one marginal at
# a time. After measuring the one-way marginals it grows a bounded-treewidth
# junction tree greedily: at each round the *exponential mechanism* privately
# picks the clique whose true marginal the model-so-far fits worst (a genuine,
# budget-spending, model-error-guided choice), then that clique is measured under
# the main mechanism. Both the selection and the measurement compose exactly into
# the same (eps, delta) as every other Track B slice (zCDP for Gaussian, pure eps
# for Laplace); the number of rounds is fixed in advance from the variable count
# and the treewidth, so the accounting is data-independent.
#
# The measured cliques form a valid junction tree by construction (each new
# clique attaches its w-variable separator to an existing clique, so the running
# intersection property holds), which means the generative model is a rooted
# junction-tree sampler -- the exact generalisation of the tree sampler's
# `root_marginal + P(child | parent)` -- and needs no iterative (IPF / PGM)
# inference. At treewidth 1 it reduces to a spanning tree; at treewidth 2 its
# cliques are triangles, so it can hold three-way interactions a tree cannot.

# Per-round exponential-mechanism epsilon for `n_rounds` selections sharing a
# `budget_frac` slice of the budget. Laplace: the pure-eps slice splits evenly and
# adds. Gaussian: work in zCDP -- allocate `budget_frac` of rho, split it evenly,
# and invert the conservative pure-eps -> (eps^2 / 2)-zCDP bound
# (eps = sqrt(2 * rho_per_round)); this never under-charges the true (eps^2 / 8)-zCDP
# cost of a bounded-range exponential mechanism, matching `dp_quantile_eps()`.
dp_select_eps <- function(dp, n_rounds, budget_frac) {
  n_rounds <- as.integer(n_rounds)
  if (n_rounds <= 0L) return(0)
  if (dp$mechanism == "laplace") {
    (budget_frac * dp$epsilon) / n_rounds
  } else {
    rho_sel <- budget_frac * zcdp_rho_for(dp$epsilon, dp$delta)
    sqrt(2 * (rho_sel / n_rounds))
  }
}

# Exponential-mechanism pick of one index from a vector of `scores` (higher =
# better) with score sensitivity `sens` and privacy parameter `eps`. The pick
# probability is proportional to exp(eps * score / (2 * sens)); computed in a
# shifted, normalised way for numerical stability. eps <= 0 falls back to a
# uniform draw (no budget for selection).
dp_exp_select <- function(scores, eps, sens) {
  k <- length(scores)
  if (k <= 1L) return(1L)
  if (eps <= 0 || sens <= 0) return(sample.int(k, 1L))
  w <- (eps / (2 * sens)) * scores
  w <- w - max(w)
  p <- exp(w)
  s <- sum(p)
  if (!is.finite(s) || s <= 0) return(sample.int(k, 1L))
  sample.int(k, 1L, prob = p / s)
}

# Integer count array (dim = nbins[vars]) of the true joint over `vars` (in the
# given order), NA-paired rows dropped. Used only inside the exponential-mechanism
# score, which is why it reads the un-noised data -- that read is exactly what the
# selection budget pays for (its sensitivity is `cap`).
dp_true_joint_array <- function(codes, vars, nbins) {
  comps <- lapply(vars, function(v) codes[[v]])
  ok <- Reduce(`&`, lapply(comps, function(x) !is.na(x)))
  comps <- lapply(comps, function(x) x[ok])
  dims <- nbins[vars]
  idx <- dp_combo_index(comps, dims)
  array(tabulate(idx, nbins = prod(dims)), dim = dims)
}

# Marginalise a stored clique count array `A` (dims over `all_vars`, in that
# order) down to the subset `S`, returning a probability vector in column-major
# order of `S` (first element of `S` varies fastest). All-zero -> uniform.
dp_array_margin <- function(A, all_vars, S) {
  pos <- match(S, all_vars)
  m <- if (length(pos) == length(dim(A))) A else apply(A, pos, sum)
  v <- as.vector(m)
  s <- sum(v)
  if (s <= 0) rep(1 / length(v), length(v)) else v / s
}

# Decode a column-major linear index (1..prod(dims)) back to per-dimension codes
# (first dimension varies fastest); inverse of dp_combo_index(). Returns an
# integer matrix with length(idx) rows and length(dims) columns.
dp_decode_combo <- function(idx, dims) {
  m <- length(dims)
  out <- matrix(0L, length(idx), m)
  r <- idx - 1L
  for (i in seq_len(m)) {
    out[, i] <- as.integer(r %% dims[i]) + 1L
    r <- r %/% dims[i]
  }
  out
}

# Fit the adaptive junction-tree model. `calib` calibrates the measurement noise;
# `sel_eps` is the per-round exponential-mechanism epsilon; `w` is the treewidth
# (clique size - 1, already capped to d - 1); `cap` is the per-person row cap (the
# selection score's sensitivity). Returns a model with kind = "junction".
dp_fit_model_adaptive <- function(codes, nbins, dp, calib, w, sel_eps, cap) {
  d <- length(codes)
  vars <- names(codes)
  idx_all <- seq_len(d)

  # One-way marginals, measured under the shared measurement noise.
  c1 <- lapply(idx_all, function(i) pmax(calib$add_noise(tabulate(codes[[i]], nbins[i])), 0))
  names(c1) <- vars
  p1 <- lapply(c1, dp_normalise)
  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))

  # Expected absolute measurement noise per cell -- a data-independent penalty
  # that discourages selecting a large, sparse clique whose measurement would be
  # swamped by noise (the AIM quality-score correction).
  noise_abs <- if (calib$mechanism == "laplace") calib$scale else calib$sigma * sqrt(2 / pi)

  # Score a candidate clique over variables `cvars` whose *new* variable set (not
  # yet jointly represented by the model) is `new_vars`; `ref_prob` is the model's
  # current probability over `cvars` in the SAME order (a product of the already
  # measured joint over the covered part and the one-ways of the new part). The
  # score is the L1 gap between the true counts and the model's expected counts,
  # minus the measurement-noise penalty. L1 is order-invariant, so `cvars` need
  # not be globally sorted -- only consistent between true and reference.
  score_candidate <- function(cvars, ref_prob) {
    true_cnt <- as.vector(dp_true_joint_array(codes, cvars, nbins))
    ref_cnt <- ref_prob * n_est
    ncell <- length(true_cnt)
    sum(abs(true_cnt - ref_cnt)) - ncell * noise_abs
  }

  covered <- logical(d)
  cliques <- list()       # each: list(vars = sorted ints, sep, new, root joint/cond)
  arrays  <- list()       # measured (noisy, clipped) count arrays, dims over vars

  # Measure one clique (variables `cvars`, any order) under the shared noise, in
  # sorted-variable order, and return its count array.
  measure_clique <- function(cvars) {
    sv <- sort(cvars)
    comps <- lapply(sv, function(v) codes[[v]])
    ok <- Reduce(`&`, lapply(comps, function(x) !is.na(x)))
    comps <- lapply(comps, function(x) x[ok])
    dims <- nbins[sv]
    cnt <- tabulate(dp_combo_index(comps, dims), nbins = prod(dims))
    noisy <- pmax(calib$add_noise(cnt), 0)
    array(noisy, dim = dims)
  }

  # --- Seed clique: the best (w + 1)-subset, model reference = product of one-ways.
  seed_sets <- utils::combn(d, w + 1L, simplify = FALSE)
  seed_scores <- vapply(seed_sets, function(cs) {
    ref <- Reduce(function(a, b) as.vector(outer(a, b)), lapply(cs, function(v) p1[[v]]))
    score_candidate(cs, ref)
  }, numeric(1))
  seed <- sort(seed_sets[[dp_exp_select(seed_scores, sel_eps, cap)]])
  A <- measure_clique(seed)
  arrays[[1L]] <- A
  cliques[[1L]] <- list(vars = seed, sep = integer(0), new = seed,
                        joint = as.vector(A) / max(sum(A), 1e-12))
  covered[seed] <- TRUE

  # --- Grow: each round attaches one uncovered variable to a w-subset that lies
  # inside an existing clique, chosen by the exponential mechanism.
  while (any(!covered)) {
    # Enumerate the w-variable separators available inside the current cliques.
    sep_list <- list()
    for (ci in seq_along(cliques)) {
      cv <- cliques[[ci]]$vars
      subs <- if (length(cv) == w) list(cv) else utils::combn(cv, w, simplify = FALSE)
      for (s in subs) sep_list[[length(sep_list) + 1L]] <- list(S = sort(s), ci = ci)
    }
    # Deduplicate separators (keep first hosting clique for the model marginal).
    keys <- vapply(sep_list, function(e) paste(e$S, collapse = "-"), character(1))
    sep_list <- sep_list[!duplicated(keys)]

    uncovered <- which(!covered)
    cand <- list()
    cand_scores <- numeric(0)
    for (e in sep_list) {
      S <- e$S
      pS <- dp_array_margin(arrays[[e$ci]], cliques[[e$ci]]$vars, S)
      for (v in uncovered) {
        cvars <- c(S, v)                              # consistent true/ref order
        ref <- as.vector(outer(pS, p1[[v]]))
        cand[[length(cand) + 1L]] <- list(S = S, v = v)
        cand_scores <- c(cand_scores, score_candidate(cvars, ref))
      }
    }
    pick <- cand[[dp_exp_select(cand_scores, sel_eps, cap)]]
    S <- pick$S; v <- pick$v
    cvars <- sort(c(S, v))
    A <- measure_clique(cvars)
    # P(v | S): reorganise the sorted-clique array into [S-cells (col-major), v].
    posS <- match(S, cvars); posV <- match(v, cvars)
    B <- apply(A, c(posS, posV), sum)                 # dims: S... then v
    cond <- matrix(as.vector(B), prod(nbins[S]), nbins[v])
    cond <- t(apply(cond, 1L, dp_normalise))
    arrays[[length(arrays) + 1L]] <- A
    cliques[[length(cliques) + 1L]] <- list(vars = cvars, sep = S, new = v, cond = cond)
    covered[v] <- TRUE
  }

  list(kind = "junction", vars = vars, nbins = nbins, cliques = cliques,
       marginals = p1, n_est = n_est, treewidth = w)
}

# Inverse-variance ("precision-weighted") combination of two independent noisy
# measurements of the same array. Both `A1` and `A2` estimate the same true
# counts, with per-cell noise parameters `p1` and `p2` (Gaussian sd, or Laplace
# scale -- either way the variance is proportional to the square, so the optimal
# weights are 1 / param^2). Measuring the same marginal twice and combining this
# way is exactly equivalent, in privacy cost, to one measurement with the summed
# budget (rho or eps adds), which is what lets the annealer spend surplus budget
# by re-measuring a clique. Returns the combined estimate and its effective noise
# parameter (`sd`), which satisfies 1 / sd^2 = 1 / p1^2 + 1 / p2^2.
dp_combine_gaussian <- function(A1, p1, A2, p2) {
  w1 <- 1 / p1^2
  w2 <- 1 / p2^2
  list(value = (w1 * A1 + w2 * A2) / (w1 + w2),
       sd = sqrt(1 / (w1 + w2)))
}

# Fit the adaptive junction-tree model with AIM-style budget annealing: a
# data-adaptive round schedule instead of the fixed `d - w` rounds. The d one-way
# marginals are measured at a fair fixed share; the clique measurements and their
# exponential-mechanism selections start at a small per-round quantum (large
# noise) and *double* whenever a round's measured signal fails to beat its noise
# floor (AIM's sigma-halving rule). After the mandatory `d - w` spanning cliques
# (which guarantee every variable is covered, so the sampler stays PGM-free) any
# surplus budget is spent on extra rounds that re-measure the worst-fit existing
# clique (inverse-variance combined). The final round absorbs the exact remainder
# of both the measurement and the selection pools, so the total spend is exactly
# the pre-committed (eps, delta) over a variable number of rounds -- every
# per-round decision is post-processing of already-privatised outputs, so this is
# sound by adaptive composition. `pool_meas` / `pool_sel` are the measurement and
# selection budgets (zCDP rho for Gaussian, pure eps for Laplace). Returns
# `list(model, anneal)`: the junction model (as dp_fit_model_adaptive) plus an
# accounting record of the realised schedule.
dp_fit_model_adaptive_anneal <- function(codes, nbins, dp, w, cap,
                                         pool_meas, pool_sel) {
  d <- length(codes)
  vars <- names(codes)
  gauss <- dp$mechanism == "gaussian"

  # Budget quantum -> per-cell noise parameter, and per-round selection eps.
  noise_from   <- function(step) if (gauss) cap / sqrt(2 * step) else cap / step
  sel_eps_from <- function(step) if (gauss) sqrt(2 * step) else step
  mean_abs     <- function(np) if (gauss) np * sqrt(2 / pi) else np
  add_noise    <- function(cnt, np) if (gauss)
    cnt + stats::rnorm(length(cnt), sd = np) else cnt + rlaplace(length(cnt), np)

  # A clipped probability vector over subset S of a stored (possibly unclipped)
  # count array A whose dims are over `allv`.
  marg_prob <- function(A, allv, S) {
    pos <- match(S, allv)
    m <- if (length(pos) == length(dim(A))) A else apply(A, pos, sum)
    dp_normalise(as.vector(m))
  }
  # Measure a clique (variables in any order) at noise parameter `np`; the raw
  # (unclipped) noisy count array in sorted-variable order, so repeated
  # measurements combine without clipping bias.
  measure_clique <- function(cvars, np) {
    sv <- sort(cvars)
    comps <- lapply(sv, function(v) codes[[v]])
    ok <- Reduce(`&`, lapply(comps, function(x) !is.na(x)))
    comps <- lapply(comps, function(x) x[ok])
    dims <- nbins[sv]
    cnt <- tabulate(dp_combo_index(comps, dims), nbins = prod(dims))
    array(add_noise(cnt, np), dim = dims)
  }

  # --- Fixed, fair share for the one-way marginals; the rest anneals. --------
  n_oneway <- d
  n_span   <- d - w
  pool_one    <- pool_meas * n_oneway / (n_oneway + n_span)
  pool_clique <- pool_meas - pool_one
  q_one <- pool_one / n_oneway

  meas_rounds <- rep(q_one, n_oneway)          # per-round measurement budgets
  sel_rounds  <- numeric(0)                    # per-round selection budgets
  noise_params <- numeric(0)                   # per measurement round
  n_anneal_steps <- 0L

  np_one <- noise_from(q_one)
  c1 <- lapply(seq_len(d), function(i) add_noise(tabulate(codes[[i]], nbins[i]), np_one))
  names(c1) <- vars
  noise_params <- c(noise_params, rep(np_one, n_oneway))
  p1 <- lapply(c1, function(x) dp_normalise(x))
  n_est <- max(1, round(mean(vapply(c1, function(x) sum(pmax(x, 0)), numeric(1)))))

  # Annealed clique quanta (start small -> noisy -> the signal test can ramp).
  q_clique0 <- pool_clique / (2 * n_span)
  q_sel0    <- pool_sel   / (2 * n_span)
  q_clique  <- q_clique0
  q_sel     <- q_sel0
  clique_left <- pool_clique
  sel_left    <- pool_sel

  score_candidate <- function(cvars, ref_prob, noise_abs) {
    true_cnt <- as.vector(dp_true_joint_array(codes, cvars, nbins))
    ncell <- length(true_cnt)
    sum(abs(true_cnt - ref_prob * n_est)) - ncell * noise_abs
  }
  # AIM signal test: did the measured marginal beat its expected noise floor? If
  # not, the quantum is too small -> double both quanta for later rounds.
  signal_test <- function(new_counts, prior_counts, np) {
    floor <- mean_abs(np) * length(new_counts)
    if (sum(abs(new_counts - prior_counts)) <= floor) {
      q_clique <<- q_clique * 2
      q_sel    <<- q_sel * 2
      n_anneal_steps <<- n_anneal_steps + 1L
    }
  }

  covered <- logical(d)
  cliques <- list()
  arrays  <- list()         # raw (unclipped) combined count arrays per clique
  clique_np <- numeric(0)   # current effective noise parameter per clique

  # --- Seed clique: best (w + 1)-subset, reference = product of one-ways. -----
  np_r <- noise_from(min(q_clique, clique_left))
  noise_abs <- mean_abs(np_r)
  seed_sets <- utils::combn(d, w + 1L, simplify = FALSE)
  seed_scores <- vapply(seed_sets, function(cs) {
    ref <- Reduce(function(a, b) as.vector(outer(a, b)), lapply(cs, function(v) p1[[v]]))
    score_candidate(cs, ref, noise_abs)
  }, numeric(1))
  # Mandatory reserve: keep q_clique0 for each remaining mandatory clique round.
  q_c <- min(q_clique, clique_left - (n_span - 1L) * q_clique0)
  q_s <- min(q_sel, sel_left - (n_span - 1L) * q_sel0)
  np_r <- noise_from(q_c)
  sel_e <- sel_eps_from(q_s)
  seed <- sort(seed_sets[[dp_exp_select(seed_scores, sel_e, cap)]])
  ref_seed <- Reduce(function(a, b) as.vector(outer(a, b)),
                     lapply(seed, function(v) p1[[v]]))
  A <- measure_clique(seed, np_r)
  signal_test(as.vector(A), ref_seed * n_est, np_r)
  arrays[[1L]] <- A; clique_np[1L] <- np_r
  cliques[[1L]] <- list(vars = seed, sep = integer(0), new = seed,
                        joint = dp_normalise(as.vector(A)))
  covered[seed] <- TRUE
  meas_rounds <- c(meas_rounds, q_c); sel_rounds <- c(sel_rounds, q_s)
  noise_params <- c(noise_params, np_r)
  clique_left <- clique_left - q_c; sel_left <- sel_left - q_s

  # --- Grow the mandatory spanning cliques (one uncovered var per round). -----
  span_done <- 1L
  while (any(!covered)) {
    remaining_mand <- n_span - span_done - 1L   # mandatory clique rounds AFTER this
    q_c <- min(q_clique, clique_left - max(remaining_mand, 0L) * q_clique0)
    q_s <- min(q_sel, sel_left - max(remaining_mand, 0L) * q_sel0)
    np_r <- noise_from(q_c); sel_e <- sel_eps_from(q_s); noise_abs <- mean_abs(np_r)

    sep_list <- list()
    for (ci in seq_along(cliques)) {
      cv <- cliques[[ci]]$vars
      subs <- if (length(cv) == w) list(cv) else utils::combn(cv, w, simplify = FALSE)
      for (s in subs) sep_list[[length(sep_list) + 1L]] <- list(S = sort(s), ci = ci)
    }
    keys <- vapply(sep_list, function(e) paste(e$S, collapse = "-"), character(1))
    sep_list <- sep_list[!duplicated(keys)]

    uncovered <- which(!covered)
    cand <- list(); cand_scores <- numeric(0); cand_ref <- list()
    for (e in sep_list) {
      pS <- marg_prob(arrays[[e$ci]], cliques[[e$ci]]$vars, e$S)
      for (v in uncovered) {
        ref <- as.vector(outer(pS, p1[[v]]))
        cand[[length(cand) + 1L]] <- list(S = e$S, v = v)
        cand_ref[[length(cand_ref) + 1L]] <- ref
        cand_scores <- c(cand_scores, score_candidate(c(e$S, v), ref, noise_abs))
      }
    }
    pick_i <- dp_exp_select(cand_scores, sel_e, cap)
    S <- cand[[pick_i]]$S; v <- cand[[pick_i]]$v
    cvars <- sort(c(S, v))
    A <- measure_clique(cvars, np_r)
    posS <- match(S, cvars); posV <- match(v, cvars)
    B <- apply(A, c(posS, posV), sum)
    cond <- t(apply(matrix(as.vector(B), prod(nbins[S]), nbins[v]), 1L, dp_normalise))
    signal_test(as.vector(A), cand_ref[[pick_i]] * n_est, np_r)
    arrays[[length(arrays) + 1L]] <- A; clique_np[length(clique_np) + 1L] <- np_r
    cliques[[length(cliques) + 1L]] <- list(vars = cvars, sep = S, new = v, cond = cond)
    covered[v] <- TRUE
    meas_rounds <- c(meas_rounds, q_c); sel_rounds <- c(sel_rounds, q_s)
    noise_params <- c(noise_params, np_r)
    clique_left <- clique_left - q_c; sel_left <- sel_left - q_s
    span_done <- span_done + 1L
  }

  # --- Refinement: spend surplus budget re-measuring the worst-fit clique. ----
  # Recompute a clique's stored conditional / joint from its (updated) array.
  refit_clique <- function(ci) {
    cl <- cliques[[ci]]; A <- arrays[[ci]]
    if (length(cl$sep) == 0L) {
      cliques[[ci]]$joint <<- dp_normalise(as.vector(A))
    } else {
      posS <- match(cl$sep, cl$vars); posV <- match(cl$new, cl$vars)
      B <- apply(A, c(posS, posV), sum)
      cliques[[ci]]$cond <<- t(apply(matrix(as.vector(B), prod(nbins[cl$sep]),
                                            nbins[cl$new]), 1L, dp_normalise))
    }
  }
  max_refine <- 8L * d
  n_refine <- 0L
  repeat {
    if (clique_left <= 1e-9 || n_refine >= max_refine) break
    q_c <- min(q_clique, clique_left)
    last <- (clique_left - q_c < q_clique0) || (n_refine + 1L >= max_refine)
    if (last) q_c <- clique_left
    q_s <- if (last) sel_left else min(q_sel, max(sel_left - q_sel0, 0))
    q_s <- min(q_s, sel_left)
    np_r <- noise_from(q_c); sel_e <- sel_eps_from(q_s); noise_abs <- mean_abs(np_r)

    # EM-select the existing clique the model currently fits worst.
    scores <- vapply(seq_along(cliques), function(ci) {
      ref <- dp_normalise(as.vector(arrays[[ci]]))
      score_candidate(cliques[[ci]]$vars, ref, noise_abs)
    }, numeric(1))
    ci <- dp_exp_select(scores, sel_e, cap)
    cvars <- cliques[[ci]]$vars
    prior <- dp_normalise(as.vector(arrays[[ci]])) * n_est
    A_new <- measure_clique(cvars, np_r)
    signal_test(as.vector(A_new), prior, np_r)
    comb <- dp_combine_gaussian(arrays[[ci]], clique_np[ci], A_new, np_r)
    arrays[[ci]] <- comb$value; clique_np[ci] <- comb$sd
    refit_clique(ci)

    meas_rounds <- c(meas_rounds, q_c); sel_rounds <- c(sel_rounds, q_s)
    noise_params <- c(noise_params, np_r)
    clique_left <- clique_left - q_c; sel_left <- sel_left - q_s
    n_refine <- n_refine + 1L
    if (last) break
  }

  model <- list(kind = "junction", vars = vars, nbins = nbins, cliques = cliques,
                marginals = p1, n_est = n_est, treewidth = w)
  info <- list(anneal = TRUE, treewidth = w, n_cliques = n_span,
               n_rounds = length(sel_rounds), n_refine = n_refine,
               n_anneal_steps = n_anneal_steps,
               noise_min = min(noise_params), noise_max = max(noise_params),
               select_frac = dp$select_frac)
  if (gauss) {
    info$rho_meas_rounds <- meas_rounds
    info$rho_sel_rounds  <- sel_rounds
  } else {
    info$eps_meas_rounds <- meas_rounds
    info$eps_sel_rounds  <- sel_rounds
  }
  list(model = model, anneal = info)
}

# Draw one synthetic code matrix (n x d) from an adaptive junction-tree model, by
# forward sampling along the junction tree: the root clique jointly, then each
# child clique's new variable conditioned on its already-sampled separator.
dp_sample_codes_junction <- function(model, n) {
  vars <- model$vars
  nbins <- model$nbins
  d <- length(vars)
  out <- matrix(NA_integer_, n, d, dimnames = list(NULL, vars))

  root <- model$cliques[[1L]]
  cell <- dp_sample_cat(n, root$joint)
  dec <- dp_decode_combo(cell, nbins[root$vars])
  for (j in seq_along(root$vars)) out[, root$vars[j]] <- dec[, j]

  for (k in seq_along(model$cliques)[-1L]) {
    cl <- model$cliques[[k]]
    sep_codes <- lapply(cl$sep, function(s) out[, s])
    sep_idx <- dp_combo_index(sep_codes, nbins[cl$sep])
    newc <- integer(n)
    for (u in unique(sep_idx)) {
      rows <- which(sep_idx == u)
      newc[rows] <- dp_sample_cat(length(rows), cl$cond[u, ])
    }
    out[, cl$new] <- newc
  }
  out
}
