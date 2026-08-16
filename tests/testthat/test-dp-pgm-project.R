# dp_pgm_project(): project a reconciled graphical model (clique beliefs over a
# junction tree) onto an arbitrary variable subset -- the primitive behind
# model-projection candidate scoring for select = "aim". A candidate *new* pair is
# by definition not yet an edge, so its two variables usually live in different
# cliques; projecting the model onto such a pair needs a junction-tree traversal,
# not a single-clique marginalisation. These tests pin the traversal against
# analytic marginals.

# A labelled non-negative factor (the shape dp-pgm.R uses): list(vars, dims, vals)
# with vals in column-major order, first variable fastest.
mk_factor <- function(vars, dims, vals)
  list(vars = as.integer(vars), dims = as.integer(dims), vals = as.numeric(vals))

test_that("projection onto a pair inside one clique equals marginalising it", {
  # Single triangle clique {1,2,3}; project {1,3}.
  set.seed(1)
  vals <- runif(8); vals <- vals / sum(vals)
  bel <- list(mk_factor(1:3, c(2, 2, 2), vals))
  nbins <- c(2L, 2L, 2L)
  got <- flexsynth:::dp_pgm_project(list(1:3), list(), bel, c(1L, 3L), nbins)
  # Reference: marginalise the clique array over variable 2.
  A <- array(vals, dim = c(2, 2, 2))
  want <- as.vector(apply(A, c(1, 3), sum))
  expect_equal(got$vars, c(1L, 3L))
  expect_equal(got$vals, want, tolerance = 1e-10)
})

test_that("projection across a chain matches the analytic A-C marginal", {
  # Two cliques {A,B}={1,2} and {B,C}={2,3}, joined on separator B. Beliefs are a
  # consistent Markov chain; the model's marginal over the non-adjacent pair (A,C)
  # is sum_B P(A,B) P(C | B).
  pAB <- c(0.40, 0.10, 0.10, 0.40)          # [A,B], A fastest; B-margin = (.5,.5)
  pBC <- c(0.45, 0.05, 0.05, 0.45)          # [B,C], B fastest; B-margin = (.5,.5)
  bel <- list(mk_factor(c(1, 2), c(2, 2), pAB),
              mk_factor(c(2, 3), c(2, 2), pBC))
  nbins <- c(2L, 2L, 2L)
  got <- flexsynth:::dp_pgm_project(list(c(1L, 2L), c(2L, 3L)), list(c(1L, 2L)),
                                    bel, c(1L, 3L), nbins)
  # Analytic P(A,C), A fastest: computed by hand from the tables above.
  want <- c(0.37, 0.13, 0.13, 0.37)
  expect_equal(got$vars, c(1L, 3L))
  expect_equal(got$vals, want, tolerance = 1e-10)
})

test_that("projection across an empty separator is the independent product", {
  # Disconnected measured components: cliques {1,2} and {3,4} share no variable,
  # so the junction-tree edge carries an empty separator and the model implies
  # 1 _||_ 3. The projected pair marginal must be P(1) (x) P(3).
  p12 <- c(0.30, 0.20, 0.10, 0.40)          # [1,2], var1 fastest
  p34 <- c(0.15, 0.35, 0.25, 0.25)          # [3,4], var3 fastest
  bel <- list(mk_factor(c(1, 2), c(2, 2), p12),
              mk_factor(c(3, 4), c(2, 2), p34))
  nbins <- c(2L, 2L, 2L, 2L)
  got <- flexsynth:::dp_pgm_project(list(c(1L, 2L), c(3L, 4L)), list(c(1L, 2L)),
                                    bel, c(1L, 3L), nbins)
  p1 <- c(sum(p12[c(1, 3)]), sum(p12[c(2, 4)]))   # marginal of var 1
  p3 <- c(sum(p34[c(1, 3)]), sum(p34[c(2, 4)]))   # marginal of var 3
  want <- as.vector(outer(p1, p3))                # var1 fastest
  expect_equal(got$vals, want, tolerance = 1e-10)
})
