# Private-PGM inference (belief propagation + entropic mirror descent) on top of
# the flat Track B engine. Covers: the factor primitives; exact reconciliation of
# a set of consistent marginals; mutual-consistency enforcement on inconsistent
# noisy marginals (the payoff the local, per-marginal construction cannot give);
# monotone loss decrease; the estimator argument defaults to "local" and leaves
# the existing paths byte-identical; validation of the new argument; that the
# tree + PGM and adaptive + PGM releases run and stay budget-neutral (the same
# measured marginals, so the same (eps, delta) accounting); reproducibility; and
# the prints.

# ---- factor primitives -----------------------------------------------------

test_that("dp_factor_marginalize sums out the right axes", {
  nb <- c(2L, 3L)                                  # vars 1 (2 levels), 2 (3 levels)
  M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)       # M[a, b], col-major over (1, 2)
  f <- list(vars = c(1L, 2L), dims = nb, vals = as.vector(M))
  # marginalise to var 1 -> rowSums; to var 2 -> colSums.
  m1 <- flexsynth:::dp_factor_marginalize(f, 1L, nb)
  m2 <- flexsynth:::dp_factor_marginalize(f, 2L, nb)
  expect_identical(m1$vars, 1L)
  expect_equal(m1$vals, rowSums(M))
  expect_identical(m2$vars, 2L)
  expect_equal(m2$vals, colSums(M))
  # marginalise to the empty set -> the grand total.
  m0 <- flexsynth:::dp_factor_marginalize(f, integer(0), nb)
  expect_equal(m0$vals, sum(M))
})

test_that("dp_factor_multiply broadcasts over the union of variables", {
  nb <- c(2L, 2L, 2L)
  f <- list(vars = c(1L, 2L), dims = c(2L, 2L),
            vals = c(1, 2, 3, 4))                  # f[a, b]
  g <- list(vars = c(2L, 3L), dims = c(2L, 2L),
            vals = c(10, 20, 30, 40))              # g[b, c]
  h <- flexsynth:::dp_factor_multiply(f, g, nb)    # h[a, b, c] = f[a,b] * g[b,c]
  expect_identical(h$vars, c(1L, 2L, 3L))
  H <- array(h$vals, dim = c(2, 2, 2))
  fA <- array(f$vals, c(2, 2)); gB <- array(g$vals, c(2, 2))
  for (a in 1:2) for (b in 1:2) for (cc in 1:2)
    expect_equal(H[a, b, cc], fA[a, b] * gB[b, cc])
})

# ---- reconciliation core ---------------------------------------------------

# A first-order Markov chain A -> B -> C over binary variables, as an exact joint.
chain_joint <- function() {
  pa  <- c(0.6, 0.4)                               # P(A)
  pbA <- rbind(c(0.8, 0.2), c(0.2, 0.8))           # P(B | A) rows = A
  pcB <- rbind(c(0.7, 0.3), c(0.3, 0.7))           # P(C | B) rows = B
  T <- array(0, dim = c(2, 2, 2))
  for (a in 1:2) for (b in 1:2) for (cc in 1:2)
    T[a, b, cc] <- pa[a] * pbA[a, b] * pcB[b, cc]
  T
}
# The chain's clique decomposition in the sampler's shape: root {A}, then
# {A,B} (P(B|A)) and {B,C} (P(C|B)).
chain_cliques <- function() list(
  list(vars = 1L,        sep = integer(0), new = 1L),
  list(vars = c(1L, 2L), sep = 1L,         new = 2L),
  list(vars = c(2L, 3L), sep = 2L,         new = 3L)
)
# Implied joint P(A,B,C) reconstructed from a reconciled junction model.
implied_joint <- function(cl, nb) {
  T <- array(0, dim = nb)
  root <- cl[[1]]$joint                            # P(A)
  c2 <- cl[[2]]$cond                               # P(B|A) rows = A
  c3 <- cl[[3]]$cond                               # P(C|B) rows = B
  for (a in 1:nb[1]) for (b in 1:nb[2]) for (cc in 1:nb[3])
    T[a, b, cc] <- root[a] * c2[a, b] * c3[b, cc]
  T
}

test_that("PGM recovers a consistent set of marginals exactly", {
  T <- chain_joint(); nb <- c(2L, 2L, 2L)
  PA <- apply(T, 1, sum); PB <- apply(T, 2, sum); PC <- apply(T, 3, sum)
  PAB <- apply(T, c(1, 2), sum); PBC <- apply(T, c(2, 3), sum)
  meas <- list(
    list(vars = 1L, y = PA), list(vars = 2L, y = PB), list(vars = 3L, y = PC),
    list(vars = c(1L, 2L), y = as.vector(PAB)),
    list(vars = c(2L, 3L), y = as.vector(PBC))
  )
  res <- flexsynth:::dp_pgm_reconcile(chain_cliques(), meas, nb,
                                      c("A", "B", "C"), n_est = 1000)
  # Consistent inputs -> the loss is driven to (essentially) zero ...
  expect_lt(res$loss, 1e-8)
  expect_lt(res$loss, res$loss0)
  # ... and the model implies exactly the true chain joint.
  expect_equal(implied_joint(res$cliques, nb), T, tolerance = 1e-5)
})

test_that("PGM reconciles inconsistent marginals into one consistent model", {
  T <- chain_joint(); nb <- c(2L, 2L, 2L)
  PA <- apply(T, 1, sum); PC <- apply(T, 3, sum)
  PAB <- apply(T, c(1, 2), sum); PBC <- apply(T, c(2, 3), sum)
  # A one-way for B that disagrees with both edge marginals (as noise would make
  # it): the raw measurements are mutually inconsistent.
  PB_bad <- c(0.3, 0.7)
  PB_from_ab <- colSums(PAB); PB_from_bc <- rowSums(PBC)
  expect_gt(sum(abs(PB_bad - PB_from_ab)), 0.1)    # the test is meaningful
  meas <- list(
    list(vars = 1L, y = PA), list(vars = 2L, y = PB_bad), list(vars = 3L, y = PC),
    list(vars = c(1L, 2L), y = as.vector(PAB)),
    list(vars = c(2L, 3L), y = as.vector(PBC))
  )
  res <- flexsynth:::dp_pgm_reconcile(chain_cliques(), meas, nb,
                                      c("A", "B", "C"), n_est = 1000)
  # The reconciled model is internally consistent: P(B) implied by the {A,B}
  # clique equals P(B) implied by the {B,C} clique (a single joint distribution).
  muAB <- res$mu[[2]]; muBC <- res$mu[[3]]
  pb_ab <- flexsynth:::dp_factor_marginalize(muAB, 2L, nb)$vals
  pb_bc <- flexsynth:::dp_factor_marginalize(muBC, 2L, nb)$vals
  expect_equal(pb_ab, pb_bc, tolerance = 1e-6)
  # And the implied joint is a valid distribution.
  Timp <- implied_joint(res$cliques, nb)
  expect_equal(sum(Timp), 1, tolerance = 1e-8)
  expect_true(all(Timp >= -1e-10))
})

# ---- estimator argument + back-compat --------------------------------------

test_that("estimator defaults to local and leaves the tree path byte-identical", {
  expect_identical(dp_control(epsilon = 1, delta = 1e-6,
                              mechanism = "gaussian")$estimator, "local")
  df <- data.frame(id = 1:400,
                   x = factor(sample(0:1, 400, TRUE)),
                   y = factor(sample(0:2, 400, TRUE)),
                   z = factor(sample(0:1, 400, TRUE)))
  d0 <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian")
  d1 <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                   estimator = "local")
  r0 <- as.data.frame(synth(df, structure = ~ id, privacy = d0, seed = 5))
  r1 <- as.data.frame(synth(df, structure = ~ id, privacy = d1, seed = 5))
  expect_identical(r0, r1)
})

test_that("dp_control validates estimator", {
  # pgm needs a tree (or adaptive) to reconcile -- not the independent model.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          dependence = "independent", estimator = "pgm"),
               "estimator|pgm|tree")
  # pgm is defined for a degree-1 tree; the Bayesian-network degree is separate.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          degree = 2, estimator = "pgm"), "estimator|pgm|degree")
  # structure_frac restructures the tree budget; not (yet) combined with pgm.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          structure_frac = 0.3, estimator = "pgm"),
               "estimator|pgm|structure")
  # annealed adaptive already refines by re-measuring; not combined with pgm.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          select = "adaptive", anneal = TRUE, estimator = "pgm"),
               "estimator|pgm|anneal")
  # the two valid homes for pgm.
  expect_identical(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                              estimator = "pgm")$estimator, "pgm")
  expect_identical(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                              select = "adaptive", estimator = "pgm")$estimator,
                   "pgm")
})

# ---- budget neutrality (pure post-processing) ------------------------------

test_that("tree + PGM measures the same marginals as the local tree (budget-neutral)", {
  df <- data.frame(id = 1:600,
                   x = factor(sample(0:1, 600, TRUE)),
                   y = factor(sample(0:2, 600, TRUE)),
                   z = factor(sample(0:1, 600, TRUE)))
  base <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian")
  pgm  <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                     estimator = "pgm")
  a_base <- synth(df, structure = ~ id, privacy = base, seed = 1)$privacy
  a_pgm  <- synth(df, structure = ~ id, privacy = pgm,  seed = 1)$privacy
  expect_equal(a_pgm$n_marginals, a_base$n_marginals)
  expect_equal(a_pgm$epsilon, a_base$epsilon)
  expect_equal(a_pgm$rho, a_base$rho, tolerance = 1e-10)
  expect_identical(a_pgm$estimator, "pgm")
})

test_that("adaptive + PGM measures the same marginals as local adaptive", {
  df <- data.frame(id = 1:600,
                   x = factor(sample(0:1, 600, TRUE)),
                   y = factor(sample(0:2, 600, TRUE)),
                   z = factor(sample(0:1, 600, TRUE)))
  base <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                     select = "adaptive", treewidth = 1)
  pgm  <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                     select = "adaptive", treewidth = 1, estimator = "pgm")
  a_base <- synth(df, structure = ~ id, privacy = base, seed = 2)$privacy
  a_pgm  <- synth(df, structure = ~ id, privacy = pgm,  seed = 2)$privacy
  expect_equal(a_pgm$n_marginals, a_base$n_marginals)
  expect_identical(a_pgm$estimator, "pgm")
})

# ---- end-to-end release ----------------------------------------------------

test_that("tree + PGM produces a well-formed release", {
  df <- data.frame(id = 1:500,
                   x = factor(sample(c("lo", "hi"), 500, TRUE)),
                   y = factor(sample(0:2, 500, TRUE)),
                   z = factor(sample(0:1, 500, TRUE)))
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   estimator = "pgm")
  out <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 3))
  expect_identical(sort(names(out)), sort(names(df)))
  expect_true(all(levels(out$x) == levels(df$x)))
  expect_true(all(as.character(out$y) %in% as.character(0:2)))
})

test_that("PGM synthesis is reproducible", {
  df <- data.frame(id = 1:400,
                   x = factor(sample(0:1, 400, TRUE)),
                   y = factor(sample(0:1, 400, TRUE)),
                   z = factor(sample(0:2, 400, TRUE)))
  dp <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                   estimator = "pgm")
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 42))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 42))
  expect_identical(a, b)
})

# ---- gates -----------------------------------------------------------------

test_that("estimator = pgm is refused on longitudinal and linked DP releases", {
  df <- data.frame(id = rep(1:200, each = 2),
                   visit = rep(1:2, times = 200),
                   x = factor(sample(0:1, 400, TRUE)),
                   y = factor(sample(0:1, 400, TRUE)))
  dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                   estimator = "pgm", max_rows_per_person = 2)
  expect_error(synth(df, structure = ~ id / visit, privacy = dp, seed = 1),
               "estimator|pgm|PGM")

  patients <- data.frame(id = 1:60, sex = factor(sample(c("F", "M"), 60, TRUE)))
  adm <- do.call(rbind, lapply(patients$id, function(pid)
    data.frame(id = pid, aid = seq_len(1 + rpois(1, 0.6)),
               los = factor(sample(0:2, 1 + 0, TRUE)))))
  adm$los <- factor(sample(0:2, nrow(adm), TRUE))
  expect_error(
    synth_linked(tables = list(patients = patients, admissions = adm),
                 structures = list(patients = ~ id, admissions = ~ id / aid),
                 keys = list(patients = "id", admissions = c("id", "aid")),
                 privacy = dp, seed = 1),
    "estimator|pgm|PGM")
})

# ---- prints ----------------------------------------------------------------

test_that("prints surface the PGM estimator", {
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   estimator = "pgm")
  expect_output(print(dp), "[Pp]rivate-?PGM|reconcil")

  df <- data.frame(id = 1:400,
                   x = factor(sample(0:1, 400, TRUE)),
                   y = factor(sample(0:1, 400, TRUE)),
                   z = factor(sample(0:2, 400, TRUE)))
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 1)$privacy
  expect_output(print(ac), "[Pp]rivate-?PGM|reconcil")
})
