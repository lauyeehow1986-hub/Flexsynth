# Internal: "Full AIM" adaptive marginal selection for the flat Track B engine.
# Not exported. See McKenna et al., "AIM: An Adaptive and Iterative Mechanism
# for Differentially Private Synthetic Data".
#
# The bounded-treewidth `select = "adaptive"` selector (dp-adaptive.R) is
# deliberately constrained: every new clique attaches its separator to an
# existing clique and covers one *uncovered* variable, so the measured cliques
# always form a valid junction tree and the model is forward-sampled with no
# inference. The price is that it can never measure a marginal between two
# variables that are *both already in the model* -- so it cannot represent a loop
# of pairwise dependence (a cycle A-B-C-A), which no tree-shaped junction
# structure holds at any budget.
#
# Full AIM removes that constraint on *selection*: the exponential mechanism may
# pick any pairwise marginal, including a loopy one. A loopy measured set has no
# forward junction-tree sampler, so Private-PGM inference (dp-pgm.R) stops being
# optional and becomes the estimator: the whole measured set (the d one-way
# marginals plus the selected pairs) is reconciled into one graphical-model
# distribution over a *triangulated* junction tree, and that model is sampled
# directly (dp_sample_codes_pgm). Every variable already carries a one-way
# marginal, so coverage is automatic and no forced spanning phase is needed.
#
# Budget: selection (the exponential-mechanism rounds) and measurement compose
# into exactly the same (eps, delta) as any other flat slice; the number of
# rounds is fixed in advance from the variable count and the treewidth, so the
# accounting is data-independent. Reconciliation is pure post-processing of the
# already-noised marginals, so it spends no extra budget.

# Maximal-clique triangulation of the interaction graph implied by a set of
# measured marginals. `d` is the variable count; `edges` a list of undirected
# c(i, j) pairs (the measured pairwise marginals; one-way marginals add no edge).
# Triangulates by a min-fill elimination ordering, extracts the maximal cliques
# of the resulting chordal graph, and builds their junction tree by a
# maximum-weight spanning tree on separator size (a valid JT satisfying the
# running-intersection property for a chordal graph). Every variable ends up in
# at least one clique -- an isolated variable becomes its own singleton clique.
# Returns list(cliques = list of sorted integer vectors, edges = junction-tree
# edges as c(a, b) clique-index pairs, width = maximum clique size in variables).
dp_triangulate <- function(d, edges) {
  G <- matrix(FALSE, d, d)
  for (e in edges) { a <- e[1L]; b <- e[2L]
    if (a != b) { G[a, b] <- TRUE; G[b, a] <- TRUE } }

  remaining <- rep(TRUE, d)
  rec <- list()
  while (any(remaining)) {
    verts <- which(remaining)
    # Fill count for eliminating each remaining vertex: non-adjacent neighbour
    # pairs that eliminating it would have to connect.
    fill <- vapply(verts, function(v) {
      nb <- which(G[v, ] & remaining)
      k <- length(nb)
      if (k <= 1L) return(0L)
      miss <- 0L
      for (a in seq_len(k - 1L)) for (b in (a + 1L):k)
        if (!G[nb[a], nb[b]]) miss <- miss + 1L
      miss
    }, integer(1))
    deg <- vapply(verts, function(v) sum(G[v, ] & remaining), integer(1))
    # Min-fill, tie-break on min degree then min index (verts is ascending).
    ord <- order(fill, deg, verts)
    v <- verts[ord[1L]]
    nb <- which(G[v, ] & remaining)
    rec[[length(rec) + 1L]] <- sort(c(v, nb))
    if (length(nb) >= 2L)                          # add the fill edges
      for (a in seq_len(length(nb) - 1L)) for (b in (a + 1L):length(nb)) {
        G[nb[a], nb[b]] <- TRUE; G[nb[b], nb[a]] <- TRUE }
    remaining[v] <- FALSE
  }

  # Keep only the maximal recorded cliques (drop exact duplicates, then any that
  # are a strict subset of another).
  keys <- vapply(rec, function(c) paste(c, collapse = "-"), character(1))
  rec <- rec[!duplicated(keys)]
  K0 <- length(rec)
  keep <- rep(TRUE, K0)
  if (K0 > 1L) for (i in seq_len(K0)) for (j in seq_len(K0))
    if (i != j && length(rec[[i]]) < length(rec[[j]]) &&
        all(rec[[i]] %in% rec[[j]])) { keep[i] <- FALSE; break }
  cliques <- rec[keep]

  K <- length(cliques)
  if (K <= 1L) {
    jt <- list()
  } else {
    W <- matrix(0, K, K)
    for (i in seq_len(K - 1L)) for (j in (i + 1L):K)
      W[i, j] <- W[j, i] <- length(intersect(cliques[[i]], cliques[[j]]))
    jt <- dp_max_spanning_tree(W)                  # undirected c(a, b) pairs
  }
  list(cliques = cliques, edges = jt,
       width = max(vapply(cliques, length, integer(1))))
}

# Candidate-scoring reference for the AIM selector: the count array a candidate
# pair is compared against when the exponential mechanism scores it. `scoring`
# picks between the two references, both read only from already-privatised
# quantities (so scoring stays pure post-processing -- the selection sensitivity
# and the (eps, delta) are identical either way):
#   "independence" -- the product of the two one-way marginals (n_est * p_i (x)
#     p_j). Cheap and data-independent-given-the-noisy-one-ways, but it ignores any
#     correlation the model already implies between i and j through other measured
#     marginals, so a loopy pair looks maximally surprising even once the model
#     explains it.
#   "model" -- the current reconciled Private-PGM model's own marginal over the
#     pair (n_est * project(i, j)). This is AIM's actual quality function: it stops
#     re-selecting a pair the model already fits and steers the budget to the
#     genuinely worst-fit interaction. Requires reconciling the measured-so-far set
#     each round (post-processing) and projecting it onto the candidate pair, which
#     for a not-yet-measured pair crosses cliques (see dp_pgm_project).
# Returns a function(i, j) giving the reference count vector in sorted-pair order.
dp_aim_reference <- function(scoring, d, nbins, idx_all, c1, p1, n_est,
                             measured, edges) {
  if (!identical(scoring, "model"))
    return(function(i, j) as.vector(outer(p1[[i]], p1[[j]])) * n_est)
  tri <- dp_triangulate(d, edges)
  meas <- c(
    lapply(idx_all, function(i) list(vars = i, y = pmax(c1[[i]], 0))),
    lapply(measured, function(m) list(vars = m$vars, y = as.vector(m$array)))
  )
  opt <- dp_pgm_optimize(tri$cliques, tri$edges, meas, nbins, max_iter = 200L)
  function(i, j)
    dp_pgm_project(tri$cliques, tri$edges, opt$bel, sort(c(i, j)), nbins)$vals *
      n_est
}

# Fit the Full AIM model. `calib` calibrates the measurement noise; `sel_eps` is
# the per-round exponential-mechanism epsilon; `w` is the treewidth (max clique
# size - 1) that bounds the triangulated model; `cap` is the per-person row cap
# (the selection score's sensitivity). Measures the d one-way marginals, then
# runs `min(choose(d, 2), w * (d - 1))` selection rounds -- each round the
# exponential mechanism picks a pairwise marginal, allowing loops, but rejecting
# any new pair whose triangulation would exceed the treewidth cap; when no
# cap-respecting new pair remains it re-measures the worst-fit existing pair
# (inverse-variance combined), so the round count -- and thus the accounting --
# stays exact. Then triangulates the measured set, reconciles it with the shared
# Private-PGM core, and returns a `kind = "pgm"` model sampled by
# dp_sample_codes_pgm().
dp_fit_model_aim <- function(codes, nbins, dp, calib, w, sel_eps, cap) {
  d <- length(codes)
  vars <- names(codes)
  idx_all <- seq_len(d)
  scoring <- if (is.null(dp$scoring)) "independence" else dp$scoring

  # One-way marginals under the shared measurement noise (clipped for n_est).
  c1 <- lapply(idx_all, function(i) pmax(calib$add_noise(tabulate(codes[[i]], nbins[i])), 0))
  names(c1) <- vars
  p1 <- lapply(c1, dp_normalise)
  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))
  base_np <- if (calib$mechanism == "laplace") calib$scale else calib$sigma
  noise_abs <- if (calib$mechanism == "laplace") calib$scale else calib$sigma * sqrt(2 / pi)

  # L1 gap between a candidate pair's true joint and a reference count array,
  # minus a noise penalty. The reference is `ref_of(i, j)` -- either the one-way
  # independence product or the current reconciled model's projection onto the
  # pair (see dp_aim_reference); the final PGM reconcile fits the whole set jointly.
  score_pair <- function(i, j, ref_of) {
    true_cnt <- as.vector(dp_true_joint_array(codes, c(i, j), nbins))
    sum(abs(true_cnt - ref_of(i, j))) - length(true_cnt) * noise_abs
  }
  # Raw (unclipped) noisy joint count array of a pair, sorted-variable order, so
  # repeated measurements inverse-variance combine without clipping bias.
  measure_pair <- function(i, j) {
    sv <- sort(c(i, j))
    comps <- lapply(sv, function(v) codes[[v]])
    ok <- Reduce(`&`, lapply(comps, function(x) !is.na(x)))
    comps <- lapply(comps, function(x) x[ok])
    dims <- nbins[sv]
    cnt <- tabulate(dp_combo_index(comps, dims), nbins = prod(dims))
    array(calib$add_noise(cnt), dim = dims)
  }

  all_pairs <- if (d >= 2L) utils::combn(d, 2L, simplify = FALSE) else list()
  n_rounds <- if (d >= 2L)
    as.integer(min(d * (d - 1L) / 2L, w * (d - 1L))) else 0L

  pair_key <- function(p) paste(sort(p), collapse = "-")
  measured <- list()                    # key -> list(vars, array, np)
  edges <- list()                       # measured undirected pairs

  tri_width_with <- function(extra) {
    e2 <- c(edges, list(sort(extra)))
    dp_triangulate(d, e2)$width
  }

  for (r in seq_len(n_rounds)) {
    have <- names(measured)
    # Reference the candidates are scored against (independence product, or the
    # model reconciled from the measured-so-far set -- refreshed each round).
    ref_of <- dp_aim_reference(scoring, d, nbins, idx_all, c1, p1, n_est,
                               measured, edges)
    # Candidate new pairs: not yet measured and keeping the model within cap.
    new_cand <- Filter(function(p) {
      k <- pair_key(p)
      !(k %in% have) && tri_width_with(p) <= w + 1L
    }, all_pairs)

    if (length(new_cand) > 0L) {
      scores <- vapply(new_cand, function(p) score_pair(p[1L], p[2L], ref_of),
                       numeric(1))
      pick <- new_cand[[dp_exp_select(scores, sel_eps, cap)]]
      i <- pick[1L]; j <- pick[2L]
      A <- measure_pair(i, j)
      measured[[pair_key(pick)]] <- list(vars = sort(c(i, j)), array = A, np = base_np)
      edges[[length(edges) + 1L]] <- sort(c(i, j))
    } else {
      # Structurally saturated: re-measure the worst-fit existing pair.
      keys <- names(measured)
      scores <- vapply(keys, function(k) {
        v <- measured[[k]]$vars; score_pair(v[1L], v[2L], ref_of)
      }, numeric(1))
      k <- keys[dp_exp_select(scores, sel_eps, cap)]
      v <- measured[[k]]$vars
      A_new <- measure_pair(v[1L], v[2L])
      comb <- dp_combine_gaussian(measured[[k]]$array, measured[[k]]$np, A_new, base_np)
      measured[[k]]$array <- array(comb$value, dim = nbins[v])
      measured[[k]]$np <- comb$sd
    }
  }

  tri <- dp_triangulate(d, edges)
  cv <- tri$cliques
  meas <- c(
    lapply(idx_all, function(i) list(vars = i, y = c1[[i]])),
    lapply(measured, function(m) list(vars = m$vars, y = as.vector(m$array)))
  )
  opt <- dp_pgm_optimize(cv, tri$edges, meas, nbins)

  list(kind = "pgm", vars = vars, nbins = nbins, cliques = cv, edges = tri$edges,
       beliefs = opt$bel, marginals = p1, n_est = n_est, treewidth = w,
       estimator = "pgm", n_rounds = n_rounds)
}

# Fit the Full AIM model with AIM-style budget annealing: a data-adaptive round
# schedule instead of the fixed `min(choose(d, 2), w * (d - 1))` rounds. The d
# one-way marginals take a fair fixed share; the pairwise measurements and their
# exponential-mechanism selections start at a small per-round quantum (large
# noise) and *double* whenever a round's measured signal fails to beat its noise
# floor (AIM's sigma-halving rule). A baseline of up to `n_target` treewidth-capped
# *new* pairs is selected first; then any surplus budget re-measures the worst-fit
# already-measured pair (inverse-variance combined). The final round absorbs the
# exact remainder of both the measurement and selection pools, so the total spend
# is exactly the pre-committed (eps, delta) over a variable number of rounds --
# sound by adaptive composition (every per-round decision is post-processing of
# already-privatised outputs). The measurement and selection quanta stay strictly
# proportional through every doubling and reserve, so the two pools deplete in
# lockstep and hit zero together (the exactness argument the fixed split relies on).
# Unlike the constrained annealer this measures *loopy* pairs (bounded by the
# treewidth after triangulation), so the result has no forward sampler: the whole
# measured set is triangulated and reconciled with the shared Private-PGM core and
# returned as a `kind = "pgm"` model. `pool_meas` / `pool_sel` are the measurement
# and selection budgets (zCDP rho for Gaussian, pure eps for Laplace). Returns
# `list(model, anneal)`: the pgm model plus a record of the realised schedule.
dp_fit_model_aim_anneal <- function(codes, nbins, dp, w, cap,
                                    pool_meas, pool_sel) {
  d <- length(codes)
  vars <- names(codes)
  idx_all <- seq_len(d)
  gauss <- dp$mechanism == "gaussian"
  scoring <- if (is.null(dp$scoring)) "independence" else dp$scoring

  # Budget quantum -> per-cell noise parameter, and per-round selection eps.
  noise_from   <- function(step) if (gauss) cap / sqrt(2 * step) else cap / step
  sel_eps_from <- function(step) if (gauss) sqrt(2 * step) else step
  mean_abs     <- function(np) if (gauss) np * sqrt(2 / pi) else np
  add_noise    <- function(cnt, np) if (gauss)
    cnt + stats::rnorm(length(cnt), sd = np) else cnt + rlaplace(length(cnt), np)

  all_pairs <- if (d >= 2L) utils::combn(d, 2L, simplify = FALSE) else list()
  # Baseline new-pair count: the fixed schedule's round count, used only as the
  # one-way fair-share denominator and the phase-A reserve horizon.
  n_target <- if (d >= 2L)
    as.integer(min(d * (d - 1L) / 2L, w * (d - 1L))) else 0L

  # --- Fixed, fair share for the one-way marginals; the rest anneals. --------
  n_oneway <- d
  pool_one  <- pool_meas * n_oneway / (n_oneway + n_target)
  pool_pair <- pool_meas - pool_one
  q_one <- pool_one / n_oneway
  np_one <- noise_from(q_one)

  meas_rounds  <- rep(q_one, n_oneway)         # per-round measurement budgets
  sel_rounds   <- numeric(0)                   # per-round selection budgets
  noise_params <- rep(np_one, n_oneway)        # per measurement round
  n_anneal_steps <- 0L

  c1 <- lapply(idx_all, function(i) add_noise(tabulate(codes[[i]], nbins[i]), np_one))
  names(c1) <- vars
  p1 <- lapply(c1, dp_normalise)
  n_est <- max(1, round(mean(vapply(c1, function(x) sum(pmax(x, 0)), numeric(1)))))

  # Annealed pair quanta, kept strictly proportional to the selection quanta.
  q_pair0 <- pool_pair / (2 * max(n_target, 1L))
  q_sel0  <- pool_sel  / (2 * max(n_target, 1L))
  q_pair  <- q_pair0
  q_sel   <- q_sel0
  pair_left <- pool_pair
  sel_left  <- pool_sel

  # L1 gap between the true pair joint and a reference count array (from `ref_of`:
  # the independence product or the reconciled model's projection; see
  # dp_aim_reference), minus a noise penalty.
  score_pair <- function(i, j, noise_abs, ref_of) {
    true_cnt <- as.vector(dp_true_joint_array(codes, c(i, j), nbins))
    sum(abs(true_cnt - ref_of(i, j))) - length(true_cnt) * noise_abs
  }
  # Raw (unclipped) noisy joint count array of a pair at noise parameter `np`, in
  # sorted-variable order, so repeated measurements combine without clipping bias.
  measure_pair <- function(i, j, np) {
    sv <- sort(c(i, j))
    comps <- lapply(sv, function(v) codes[[v]])
    ok <- Reduce(`&`, lapply(comps, function(x) !is.na(x)))
    comps <- lapply(comps, function(x) x[ok])
    dims <- nbins[sv]
    cnt <- tabulate(dp_combo_index(comps, dims), nbins = prod(dims))
    array(add_noise(cnt, np), dim = dims)
  }
  # AIM signal test: did the measured marginal beat its expected noise floor? If
  # not, the quantum is too small -> double both quanta for later rounds.
  signal_test <- function(new_counts, prior_counts, np) {
    floor <- mean_abs(np) * length(new_counts)
    if (sum(abs(new_counts - prior_counts)) <= floor) {
      q_pair <<- q_pair * 2
      q_sel  <<- q_sel * 2
      n_anneal_steps <<- n_anneal_steps + 1L
    }
  }

  pair_key <- function(p) paste(sort(p), collapse = "-")
  measured <- list()                    # key -> list(vars, array, np)
  edges <- list()                       # measured undirected pairs
  tri_width_with <- function(extra)
    dp_triangulate(d, c(edges, list(sort(extra))))$width

  # --- Phase A: baseline new-pair selection (up to n_target, cap-respecting). --
  n_new <- 0L
  saturated <- FALSE
  while (n_new < n_target && !saturated) {
    have <- names(measured)
    ref_of <- dp_aim_reference(scoring, d, nbins, idx_all, c1, p1, n_est,
                               measured, edges)
    new_cand <- Filter(function(p) {
      k <- pair_key(p)
      !(k %in% have) && tri_width_with(p) <= w + 1L
    }, all_pairs)
    if (length(new_cand) == 0L) { saturated <- TRUE; break }
    remaining <- n_target - n_new - 1L         # baseline rounds AFTER this one
    q_p <- min(q_pair, pair_left - max(remaining, 0L) * q_pair0)
    q_s <- min(q_sel,  sel_left  - max(remaining, 0L) * q_sel0)
    np_r <- noise_from(q_p); sel_e <- sel_eps_from(q_s); noise_abs <- mean_abs(np_r)

    scores <- vapply(new_cand,
                     function(p) score_pair(p[1L], p[2L], noise_abs, ref_of),
                     numeric(1))
    pick <- new_cand[[dp_exp_select(scores, sel_e, cap)]]
    i <- pick[1L]; j <- pick[2L]
    A <- measure_pair(i, j, np_r)
    # Signal test against the pair's prior estimate (the reference the selector
    # scored it by): did the measurement beat its noise floor?
    signal_test(as.vector(A), ref_of(i, j), np_r)
    measured[[pair_key(pick)]] <- list(vars = sort(c(i, j)), array = A, np = np_r)
    edges[[length(edges) + 1L]] <- sort(c(i, j))
    meas_rounds <- c(meas_rounds, q_p); sel_rounds <- c(sel_rounds, q_s)
    noise_params <- c(noise_params, np_r)
    pair_left <- pair_left - q_p; sel_left <- sel_left - q_s
    n_new <- n_new + 1L
  }

  # --- Phase B: spend surplus budget re-measuring the worst-fit existing pair. -
  max_refine <- 8L * d
  n_refine <- 0L
  repeat {
    if (pair_left <= 1e-9 || n_refine >= max_refine || length(measured) == 0L) break
    q_p <- min(q_pair, pair_left)
    last <- (pair_left - q_p < q_pair0) || (n_refine + 1L >= max_refine)
    if (last) q_p <- pair_left
    q_s <- if (last) sel_left else min(q_sel, max(sel_left - q_sel0, 0))
    q_s <- min(q_s, sel_left)
    np_r <- noise_from(q_p); sel_e <- sel_eps_from(q_s); noise_abs <- mean_abs(np_r)

    ref_of <- dp_aim_reference(scoring, d, nbins, idx_all, c1, p1, n_est,
                               measured, edges)
    keys <- names(measured)
    scores <- vapply(keys, function(k) {
      v <- measured[[k]]$vars; score_pair(v[1L], v[2L], noise_abs, ref_of)
    }, numeric(1))
    k <- keys[dp_exp_select(scores, sel_e, cap)]
    v <- measured[[k]]$vars
    A_new <- measure_pair(v[1L], v[2L], np_r)
    signal_test(as.vector(A_new), as.vector(measured[[k]]$array), np_r)
    comb <- dp_combine_gaussian(measured[[k]]$array, measured[[k]]$np, A_new, np_r)
    measured[[k]]$array <- array(comb$value, dim = nbins[v])
    measured[[k]]$np <- comb$sd

    meas_rounds <- c(meas_rounds, q_p); sel_rounds <- c(sel_rounds, q_s)
    noise_params <- c(noise_params, np_r)
    pair_left <- pair_left - q_p; sel_left <- sel_left - q_s
    n_refine <- n_refine + 1L
    if (last) break
  }

  # Triangulate the measured (loopy) edge set and reconcile the whole set with the
  # shared Private-PGM core, then sample it via dp_sample_codes_pgm().
  tri <- dp_triangulate(d, edges)
  cv <- tri$cliques
  meas <- c(
    lapply(idx_all, function(i) list(vars = i, y = pmax(c1[[i]], 0))),
    lapply(measured, function(m) list(vars = m$vars, y = as.vector(m$array)))
  )
  opt <- dp_pgm_optimize(cv, tri$edges, meas, nbins)

  model <- list(kind = "pgm", vars = vars, nbins = nbins, cliques = cv,
                edges = tri$edges, beliefs = opt$bel, marginals = p1,
                n_est = n_est, treewidth = w, estimator = "pgm",
                n_rounds = length(sel_rounds))
  info <- list(anneal = TRUE, treewidth = w, n_new = n_new, n_target = n_target,
               n_rounds = length(sel_rounds), n_refine = n_refine,
               n_anneal_steps = n_anneal_steps,
               noise_min = min(noise_params), noise_max = max(noise_params),
               select_frac = dp$select_frac, scoring = scoring)
  if (gauss) {
    info$rho_meas_rounds <- meas_rounds
    info$rho_sel_rounds  <- sel_rounds
  } else {
    info$eps_meas_rounds <- meas_rounds
    info$eps_sel_rounds  <- sel_rounds
  }
  list(model = model, anneal = info)
}

# Draw one synthetic code matrix (n x d) from a Full AIM (`kind = "pgm"`) model
# by forward sampling along its triangulated junction tree: root the tree at
# clique 1 and sample its joint belief; then for each child clique sample its
# *new* variables (those not shared with the parent) conditioned on its already-
# sampled separator. Unlike dp_sample_codes_junction(), separators and new-
# variable blocks are set-valued -- a clique may introduce several variables at
# once -- so cells are decoded into multiple columns.
dp_sample_codes_pgm <- function(model, n) {
  vars <- model$vars
  nbins <- model$nbins
  d <- length(vars)
  cv <- model$cliques
  bel <- model$beliefs
  K <- length(cv)
  out <- matrix(NA_integer_, n, d, dimnames = list(NULL, vars))

  # Sample a clique's new-variable block from a conditional built from its belief.
  emit <- function(c, sep, new) {
    A <- array(bel[[c]]$vals, dim = bel[[c]]$dims)
    posS <- match(sep, cv[[c]]); posV <- match(new, cv[[c]])
    if (length(sep) == 0L) {
      B <- if (length(posV) == length(dim(A))) A else apply(A, posV, sum)
      prob <- dp_normalise(as.vector(B))
      cell <- dp_sample_cat(n, prob)
      dec <- dp_decode_combo(cell, nbins[new])
      for (jj in seq_along(new)) out[, new[jj]] <<- dec[, jj]
    } else {
      B <- apply(A, c(posS, posV), sum)              # dims: sep..., then new...
      cond <- matrix(as.vector(B), prod(nbins[sep]), prod(nbins[new]))
      cond <- t(apply(cond, 1L, dp_normalise))
      sep_idx <- dp_combo_index(lapply(sep, function(s) out[, s]), nbins[sep])
      for (u in unique(sep_idx)) {
        rows <- which(sep_idx == u)
        cell <- dp_sample_cat(length(rows), cond[u, ])
        dec <- dp_decode_combo(cell, nbins[new])
        for (jj in seq_along(new)) out[rows, new[jj]] <<- dec[, jj]
      }
    }
  }

  # Root clique: sample its whole joint.
  emit(1L, integer(0), cv[[1L]])
  if (K == 1L) return(out)

  # BFS from clique 1 for a parent structure and a distribute order.
  adj <- vector("list", K)
  for (e in model$edges) { a <- e[1L]; b <- e[2L]
    adj[[a]] <- c(adj[[a]], b); adj[[b]] <- c(adj[[b]], a) }
  parent <- integer(K); visited <- logical(K)
  order <- integer(0); queue <- 1L; visited[1L] <- TRUE
  while (length(queue)) {
    a <- queue[1L]; queue <- queue[-1L]; order <- c(order, a)
    for (b in adj[[a]]) if (!visited[b]) {
      visited[b] <- TRUE; parent[b] <- a; queue <- c(queue, b) }
  }
  for (c in order[-1L]) {
    p <- parent[c]
    sep <- sort(intersect(cv[[c]], cv[[p]]))
    new <- sort(setdiff(cv[[c]], sep))
    emit(c, sep, new)
  }
  out
}
