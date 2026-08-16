# Internal: Private-PGM inference for the flat Track B engine. Not exported.
#
# Every other Track B sampler is PGM-free: it measures noisy marginals and then
# uses each one *locally* -- the tree release, for instance, takes the root's
# one-way and each Chow-Liu edge's raw noisy 2-way as P(child | parent), and
# throws away both the other one-way marginals and the fact that overlapping
# edges (each measured independently under noise) disagree. This file adds the
# reconciliation those samplers omit, following McKenna et al.'s Private-PGM /
# MST: given the whole set of noisy, mutually-inconsistent measured marginals,
# find the single graphical-model distribution that best fits *all* of them at
# once (least squares), by belief propagation on a junction tree of the measured
# cliques plus entropic mirror descent, then read the reconciled clique marginals
# back into the existing junction sampler's shape (root joint + P(new | sep)).
#
# Reconciliation is pure post-processing of already-privatised marginals, so it
# spends no extra budget -- the (eps, delta) of a `estimator = "pgm"` release is
# identical to the same release with `estimator = "local"`; only the model built
# from the noisy marginals changes.
#
# A "factor" here is a labelled non-negative array over a sorted subset of the
# variable indices: list(vars = sorted ints, dims = nbins[vars], vals = numeric
# in column-major order with the first variable varying fastest, matching
# dp_combo_index() / dp_decode_combo()).

# Sum a factor down to the subset `keep` (any order; intersected with f's vars).
# The empty set yields the grand total as a length-1 factor.
dp_factor_marginalize <- function(f, keep, nbins) {
  keep <- sort(intersect(f$vars, keep))
  if (length(keep) == 0L)
    return(list(vars = integer(0), dims = integer(0), vals = sum(f$vals)))
  if (length(keep) == length(f$vars))
    return(list(vars = f$vars, dims = f$dims, vals = f$vals))
  pos <- match(keep, f$vars)
  A <- array(f$vals, dim = f$dims)
  m <- apply(A, pos, sum)                          # dims in increasing pos order
  list(vars = keep, dims = f$dims[pos], vals = as.vector(m))
}

# Pointwise product of two factors over the union of their variables (broadcast
# over the variables each lacks). A variable-less factor acts as a scalar.
dp_factor_multiply <- function(f, g, nbins) {
  if (length(f$vars) == 0L)
    return(list(vars = g$vars, dims = g$dims, vals = g$vals * f$vals))
  if (length(g$vars) == 0L)
    return(list(vars = f$vars, dims = f$dims, vals = f$vals * g$vals))
  uv <- sort(union(f$vars, g$vars))
  udims <- as.integer(nbins[uv])
  codes <- dp_decode_combo(seq_len(prod(udims)), udims)   # cells x length(uv)
  fpos <- match(f$vars, uv); gpos <- match(g$vars, uv)
  fi <- dp_combo_index(lapply(seq_along(f$vars), function(k) codes[, fpos[k]]),
                       f$dims)
  gi <- dp_combo_index(lapply(seq_along(g$vars), function(k) codes[, gpos[k]]),
                       g$dims)
  list(vars = uv, dims = udims, vals = f$vals[fi] * g$vals[gi])
}

# Broadcast a factor up to the (superset) variable set `target`, constant across
# the added variables. Used to lift a marginal residual back onto its clique.
dp_factor_expand <- function(f, target, nbins) {
  target <- sort(as.integer(target))
  extra <- setdiff(target, f$vars)
  if (length(extra) == 0L) return(f)
  ones <- list(vars = extra, dims = as.integer(nbins[extra]),
               vals = rep(1, prod(nbins[extra])))
  dp_factor_multiply(f, ones, nbins)
}

# Exact marginal inference on a junction tree by Shafer-Shenoy message passing
# (no division). `cv` is the list of clique variable sets, `edges` the undirected
# spanning-tree edges (c(a, b) clique indices), `theta` the per-clique log
# potentials (vectors over the clique cells). Returns one normalised belief
# factor (the clique marginal of the global distribution p ~ exp(sum_c theta_c))
# per clique. Messages and beliefs are normalised for numerical stability, which
# is harmless because only marginals (defined up to scale) are used downstream.
dp_pgm_infer <- function(cv, edges, theta, nbins) {
  K <- length(cv)
  psi <- lapply(seq_len(K), function(c) {
    t <- theta[[c]]; t <- t - max(t)               # recentre before exp
    list(vars = cv[[c]], dims = as.integer(nbins[cv[[c]]]), vals = exp(t))
  })
  normed <- function(f) { s <- sum(f$vals)
    if (s > 0) f$vals <- f$vals / s else
      f$vals <- rep(1 / length(f$vals), length(f$vals)); f }
  if (K == 1L) return(list(normed(psi[[1L]])))

  adj <- vector("list", K)
  for (e in edges) { a <- e[1L]; b <- e[2L]
    adj[[a]] <- c(adj[[a]], b); adj[[b]] <- c(adj[[b]], a) }
  sepvars <- function(a, b) sort(intersect(cv[[a]], cv[[b]]))
  msg <- lapply(seq_len(K), function(i) vector("list", K))

  # BFS from clique 1 for a parent structure and a collect/distribute order.
  order <- integer(0); parent <- integer(K); visited <- logical(K)
  queue <- 1L; visited[1L] <- TRUE
  while (length(queue)) {
    a <- queue[1L]; queue <- queue[-1L]; order <- c(order, a)
    for (b in adj[[a]]) if (!visited[b]) {
      visited[b] <- TRUE; parent[b] <- a; queue <- c(queue, b) }
  }
  compute_msg <- function(a, b) {
    prod_f <- psi[[a]]
    for (k in adj[[a]]) if (k != b && !is.null(msg[[k]][[a]]))
      prod_f <- dp_factor_multiply(prod_f, msg[[k]][[a]], nbins)
    mm <- dp_factor_marginalize(prod_f, sepvars(a, b), nbins)
    s <- sum(mm$vals); if (s > 0) mm$vals <- mm$vals / s
    mm
  }
  for (a in rev(order)) { p <- parent[a]
    if (p != 0L) msg[[a]][[p]] <- compute_msg(a, p) }        # collect
  for (a in order) for (b in adj[[a]]) if (parent[b] == a)
    msg[[a]][[b]] <- compute_msg(a, b)                       # distribute

  lapply(seq_len(K), function(c) {
    bf <- psi[[c]]
    for (k in adj[[c]]) if (!is.null(msg[[k]][[c]]))
      bf <- dp_factor_multiply(bf, msg[[k]][[c]], nbins)
    normed(bf)
  })
}

# Reconcile a set of noisy measured marginals onto a fixed junction tree by
# entropic mirror descent: infer the clique marginals of p ~ exp(sum_c theta_c),
# take each measurement's marginal residual (mu - y) as the mirror gradient
# (lifted to its hosting clique), step the log potentials, and backtrack the step
# so the least-squares loss decreases -- to a fixed point. `cv` is the list of
# (sorted) clique variable sets, `edges` the junction-tree edges, `measurements`
# a list of list(vars, y) noisy marginal histograms (y in factor cell order over
# sort(vars)); each measurement must be a subset of some clique. Returns the
# reconciled clique beliefs and the initial / final loss / iteration count. This
# is the shared inference core behind both dp_pgm_reconcile() (which writes the
# beliefs back into an ancestral sampler shape) and the AIM fitter (which samples
# the beliefs directly over a general triangulated junction tree; see dp-aim.R).
dp_pgm_optimize <- function(cv, edges, measurements, nbins,
                            max_iter = 500L, tol = 1e-12) {
  K <- length(cv)
  # Normalise each measurement to a distribution and bind it to a hosting clique.
  meas <- lapply(measurements, function(m) {
    mv <- sort(as.integer(m$vars))
    host <- NA_integer_
    for (c in seq_len(K)) if (all(mv %in% cv[[c]])) { host <- c; break }
    if (is.na(host))
      stop("PGM reconcile: a measured marginal is not covered by any clique.",
           call. = FALSE)
    list(vars = mv, y = dp_normalise(as.numeric(m$y)), host = host)
  })

  theta <- lapply(cv, function(v) rep(0, prod(nbins[v])))
  loss_of <- function(bel) {
    L <- 0
    for (m in meas) {
      mu <- dp_factor_marginalize(bel[[m$host]], m$vars, nbins)$vals
      L <- L + 0.5 * sum((mu - m$y)^2)
    }
    L
  }
  bel <- dp_pgm_infer(cv, edges, theta, nbins)
  L <- loss_of(bel); L0 <- L
  alpha <- 1; iters <- 0L

  for (it in seq_len(max_iter)) {
    iters <- it
    grad <- lapply(cv, function(v) rep(0, prod(nbins[v])))
    for (m in meas) {
      mu <- dp_factor_marginalize(bel[[m$host]], m$vars, nbins)$vals
      r <- list(vars = m$vars, dims = as.integer(nbins[m$vars]), vals = mu - m$y)
      grad[[m$host]] <- grad[[m$host]] +
        dp_factor_expand(r, cv[[m$host]], nbins)$vals
    }
    step <- alpha; accepted <- FALSE; L_new <- L; b_new <- bel
    for (bt in seq_len(40L)) {
      theta_try <- lapply(seq_len(K), function(c) theta[[c]] - step * grad[[c]])
      b_try <- dp_pgm_infer(cv, edges, theta_try, nbins)
      L_try <- loss_of(b_try)
      if (L_try < L - 1e-16) {
        accepted <- TRUE; theta_new <- theta_try; b_new <- b_try; L_new <- L_try
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    rel <- (L - L_new) / max(L, 1e-12)
    theta <- lapply(theta_new, function(t) t - mean(t))   # recentre, no effect
    bel <- b_new; L <- L_new
    alpha <- min(step * 1.5, 1e3)
    if (rel < tol) break
  }
  list(bel = bel, loss0 = L0, loss = L, iters = iters)
}

# Reconcile a set of noisy measured marginals onto the clique structure of a
# fitted junction / tree model. `cliques` is the fitter's ancestral clique list
# (each list(vars, sep, new); the first has sep = integer(0)); `measurements` is
# a list of list(vars, y) noisy marginal histograms (y in the factor cell order
# over sort(vars)). Builds the junction tree over the cliques (maximum-weight
# spanning tree by shared-variable count -- a valid JT for the triangulated
# clique sets the fitters produce), runs the shared mirror-descent core, then
# writes each clique's reconciled joint (root) or P(new | sep) (child) back in
# place. Returns the updated cliques, the reconciled clique beliefs, and the
# initial / final loss.
dp_pgm_reconcile <- function(cliques, measurements, nbins, vars, n_est,
                             max_iter = 500L, tol = 1e-12) {
  K <- length(cliques)
  cv <- lapply(cliques, function(cl) sort(as.integer(cl$vars)))

  if (K == 1L) {
    edges <- list()
  } else {
    W <- matrix(0, K, K)
    for (i in seq_len(K - 1L)) for (j in (i + 1L):K)
      W[i, j] <- W[j, i] <- length(intersect(cv[[i]], cv[[j]]))
    edges <- dp_max_spanning_tree(W)
  }

  opt <- dp_pgm_optimize(cv, edges, measurements, nbins, max_iter, tol)
  bel <- opt$bel

  # Write the reconciled marginals back into the sampler's clique shape.
  cl_out <- cliques
  for (c in seq_len(K)) {
    cl <- cliques[[c]]
    if (length(cl$sep) == 0L) {
      cl_out[[c]]$joint <- dp_normalise(bel[[c]]$vals)
    } else {
      sep <- sort(as.integer(cl$sep)); new <- as.integer(cl$new)
      jf <- dp_factor_marginalize(bel[[c]], c(sep, new), nbins)
      posS <- match(sep, jf$vars); posV <- match(new, jf$vars)
      A <- array(jf$vals, dim = jf$dims)
      B <- apply(A, c(posS, posV), sum)              # dims: sep..., then new
      cond <- matrix(as.vector(B), prod(nbins[sep]), nbins[new])
      cl_out[[c]]$cond <- t(apply(cond, 1L, dp_normalise))
    }
  }
  list(cliques = cl_out, mu = bel, loss0 = opt$loss0, loss = opt$loss,
       iters = opt$iters)
}
