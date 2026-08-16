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
