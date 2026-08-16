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

  # One-way marginals under the shared measurement noise (clipped for n_est).
  c1 <- lapply(idx_all, function(i) pmax(calib$add_noise(tabulate(codes[[i]], nbins[i])), 0))
  names(c1) <- vars
  p1 <- lapply(c1, dp_normalise)
  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))
  base_np <- if (calib$mechanism == "laplace") calib$scale else calib$sigma
  noise_abs <- if (calib$mechanism == "laplace") calib$scale else calib$sigma * sqrt(2 / pi)

  # True and reference (independence) counts for a candidate pair, in sorted
  # order (first variable fastest) -- the same one-way-product reference the
  # adaptive selector uses; the final PGM reconcile fits the whole set jointly.
  score_pair <- function(i, j) {
    true_cnt <- as.vector(dp_true_joint_array(codes, c(i, j), nbins))
    ref <- as.vector(outer(p1[[i]], p1[[j]])) * n_est
    sum(abs(true_cnt - ref)) - length(true_cnt) * noise_abs
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
    # Candidate new pairs: not yet measured and keeping the model within cap.
    new_cand <- Filter(function(p) {
      k <- pair_key(p)
      !(k %in% have) && tri_width_with(p) <= w + 1L
    }, all_pairs)

    if (length(new_cand) > 0L) {
      scores <- vapply(new_cand, function(p) score_pair(p[1L], p[2L]), numeric(1))
      pick <- new_cand[[dp_exp_select(scores, sel_eps, cap)]]
      i <- pick[1L]; j <- pick[2L]
      A <- measure_pair(i, j)
      measured[[pair_key(pick)]] <- list(vars = sort(c(i, j)), array = A, np = base_np)
      edges[[length(edges) + 1L]] <- sort(c(i, j))
    } else {
      # Structurally saturated: re-measure the worst-fit existing pair.
      keys <- names(measured)
      scores <- vapply(keys, function(k) {
        v <- measured[[k]]$vars; score_pair(v[1L], v[2L])
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
