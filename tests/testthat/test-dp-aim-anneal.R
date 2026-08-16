# Annealed Full AIM (select = "aim" with anneal = TRUE): the AIM loopy-marginal
# selector + Private-PGM reconciliation driven by the data-adaptive sigma-halving
# budget schedule instead of the fixed round count. Covers control acceptance,
# exact budget composition over a *variable* round count (both mechanisms), a
# well-formed / reproducible / fully-covered release, the loop still being
# captured, that the signal test and refinement fire, and the prints.

# Three binary variables with a loop of pairwise dependence (every pair
# correlated, no tree edge-set holds all three). Local copy so this file does not
# depend on another test file's top-level definitions.
aim_loop_data <- function(n = 4000L, seed = 11L) {
  set.seed(seed)
  a <- sample(0:1, n, TRUE)
  b <- ifelse(runif(n) < 0.85, a, 1L - a)
  cc <- ifelse(runif(n) < 0.85, b, 1L - b)
  a2 <- ifelse(runif(n) < 0.7, cc, a)
  data.frame(id = seq_len(n),
             A = factor(a2), B = factor(b), C = factor(cc))
}

# ---- control acceptance ----------------------------------------------------

test_that("dp_control accepts aim + anneal and stores both", {
  dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  expect_identical(dp$select, "aim")
  expect_true(isTRUE(dp$anneal))
})

test_that("anneal is still rejected without an adaptive or aim selector", {
  # Default select = "fixed" has no per-round schedule to anneal.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          anneal = TRUE),
               "anneal")
})

# ---- exact budget composition over a variable round count ------------------

test_that("annealed AIM gaussian budget composes exactly over the round count", {
  df <- aim_loop_data(2500L, 4L)
  sf <- 0.3
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, select_frac = sf,
                   anneal = TRUE)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 2)$privacy
  rho_total <- flexsynth:::zcdp_rho_for(8, 1e-6)
  am <- ac$aim
  expect_true(isTRUE(am$anneal))
  expect_equal(sum(am$rho_meas_rounds), (1 - sf) * rho_total, tolerance = 1e-6)
  expect_equal(sum(am$rho_sel_rounds), sf * rho_total, tolerance = 1e-6)
  expect_equal(sum(am$rho_meas_rounds) + sum(am$rho_sel_rounds),
               rho_total, tolerance = 1e-6)
})

test_that("annealed AIM laplace budget composes exactly over the round count", {
  df <- aim_loop_data(2500L, 6L)
  sf <- 0.25; eps <- 9
  dp <- dp_control(epsilon = eps, mechanism = "laplace",
                   select = "aim", treewidth = 2, select_frac = sf,
                   anneal = TRUE)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 3)$privacy
  am <- ac$aim
  expect_true(isTRUE(am$anneal))
  expect_equal(sum(am$eps_meas_rounds), (1 - sf) * eps, tolerance = 1e-6)
  expect_equal(sum(am$eps_sel_rounds), sf * eps, tolerance = 1e-6)
  expect_equal(sum(am$eps_meas_rounds) + sum(am$eps_sel_rounds),
               eps, tolerance = 1e-6)
})

# ---- release shape ---------------------------------------------------------

test_that("annealed AIM release is well-formed, covered and reproducible", {
  df <- aim_loop_data(2000L, 8L)
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  expect_named(a, names(df))
  expect_false(anyNA(a))
  expect_equal(a$id, seq_len(nrow(a)))
  expect_identical(a, b)                          # reproducible
})

test_that("annealed AIM accounting reports a pgm model and the realised schedule", {
  df <- aim_loop_data(2000L, 9L)
  dp <- dp_control(epsilon = 10, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 5)$privacy
  am <- ac$aim
  expect_false(is.null(am))
  expect_true(isTRUE(am$anneal))
  # total selection rounds = distinct new pairs + refinement rounds
  expect_equal(am$n_rounds, am$n_new + am$n_refine)
  # d one-ways + the realised selection rounds were measured
  expect_equal(ac$n_marginals, length(ac$variables) + am$n_rounds)
  expect_equal(ac$epsilon, 10)
})

# ---- structure / utility ---------------------------------------------------

test_that("annealed AIM at treewidth 2 reconciles the loop into one triangle", {
  df <- aim_loop_data(4000L, 11L)
  codes <- lapply(c("A", "B", "C"), function(v) as.integer(df[[v]]))
  names(codes) <- c("A", "B", "C")
  nbins <- c(A = 2L, B = 2L, C = 2L)
  dp <- dp_control(epsilon = 40, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  total <- flexsynth:::zcdp_rho_for(40, 1e-6)
  set.seed(7)
  fit <- flexsynth:::dp_fit_model_aim_anneal(codes, nbins, dp, w = 2L, cap = 1,
                                             pool_meas = 0.75 * total,
                                             pool_sel = 0.25 * total)
  model <- fit$model
  expect_identical(model$kind, "pgm")
  # a complete graph on 3 vars triangulates to a single triangle clique
  expect_length(model$cliques, 1L)
  expect_setequal(model$cliques[[1L]], 1:3)
  # the reconciled model reproduces every pairwise marginal (a tree cannot)
  bel <- model$beliefs[[1L]]
  A <- array(bel$vals, dim = bel$dims)
  pac <- apply(A, c(1, 3), sum)
  tab <- function(i, j) prop.table(table(factor(codes[[i]], 1:2),
                                          factor(codes[[j]], 1:2)))
  expect_equal(as.numeric(pac), as.numeric(tab(1, 3)), tolerance = 0.05)
})

test_that("with surplus budget annealed AIM adds refinement rounds", {
  df <- aim_loop_data(3000L, 13L)
  codes <- lapply(c("A", "B", "C"), function(v) as.integer(df[[v]]))
  names(codes) <- c("A", "B", "C")
  nbins <- c(A = 2L, B = 2L, C = 2L)
  dp <- dp_control(epsilon = 30, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  total <- flexsynth:::zcdp_rho_for(30, 1e-6)
  set.seed(8)
  fit <- flexsynth:::dp_fit_model_aim_anneal(codes, nbins, dp, w = 2L, cap = 1,
                                             pool_meas = 0.75 * total,
                                             pool_sel = 0.25 * total)
  am <- fit$anneal
  expect_gte(am$n_refine, 1L)
  expect_equal(am$n_rounds, am$n_new + am$n_refine)
})

test_that("the sigma-halving signal test fires under a tiny budget", {
  df <- data.frame(id = 1:4000,
                   v1 = factor(sample(1:4, 4000, TRUE)),
                   v2 = factor(sample(1:4, 4000, TRUE)),
                   v3 = factor(sample(1:4, 4000, TRUE)),
                   v4 = factor(sample(1:4, 4000, TRUE)))
  dp <- dp_control(epsilon = 0.3, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  am <- synth(df, structure = ~ id, privacy = dp, seed = 7)$privacy$aim
  expect_true(isTRUE(am$anneal))
  expect_gte(am$n_anneal_steps, 1L)
})

# ---- prints ----------------------------------------------------------------

test_that("prints surface annealed AIM", {
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, anneal = TRUE)
  expect_output(print(dp), "annealed")

  df <- aim_loop_data(1500L, 15L)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 1)$privacy
  expect_output(print(ac), "annealed")
  expect_output(print(ac), "round\\(s\\)")
})
