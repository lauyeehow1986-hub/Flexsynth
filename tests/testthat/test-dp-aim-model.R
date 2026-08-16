# Model-projection candidate scoring for select = "aim". Full AIM scores each
# candidate marginal by how badly the *current reconciled model* predicts it,
# rather than against the one-way independence product. The reference is read from
# already-privatised marginals (the reconciled model), so it is pure
# post-processing: the exponential-mechanism sensitivity and the exact (eps,delta)
# are unchanged. Covers control validation, the core scoring difference, a
# well-formed release (fixed and annealed), the accounting invariance, and prints.

# Three binary variables in a chain A-B-C: A~B~C strongly, so the true A-C pair is
# correlated even though no A-C edge is measured. Local copy so this file is
# self-contained.
aim_chain_data <- function(n = 4000L, seed = 21L) {
  set.seed(seed)
  a <- sample(0:1, n, TRUE)
  b <- ifelse(runif(n) < 0.9, a, 1L - a)
  cc <- ifelse(runif(n) < 0.9, b, 1L - b)
  data.frame(id = seq_len(n), A = factor(a), B = factor(b), C = factor(cc))
}

# ---- control validation ----------------------------------------------------

test_that("dp_control accepts scoring = 'model' with select = 'aim'", {
  dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, scoring = "model")
  expect_identical(dp$select, "aim")
  expect_identical(dp$scoring, "model")
})

test_that("scoring defaults to independence and is stored on every control", {
  dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian")
  expect_identical(dp$scoring, "independence")
})

test_that("scoring = 'model' is rejected without select = 'aim'", {
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          scoring = "model"),
               "scoring")
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          select = "adaptive", scoring = "model"),
               "scoring")
})

# ---- the core scoring difference -------------------------------------------

test_that("the model reference explains an induced pair better than independence", {
  # Measure the chain edges A-B and B-C (exact, noise-free counts). The true A-C
  # pair is correlated through B; the independence product misses that, while the
  # model reconciled from {AB, BC} projects the induced A-C correlation.
  df <- aim_chain_data(6000L, 3L)
  codes <- list(A = as.integer(df$A), B = as.integer(df$B), C = as.integer(df$C))
  nbins <- c(A = 2L, B = 2L, C = 2L)
  idx_all <- 1:3
  n_est <- nrow(df)
  c1 <- lapply(idx_all, function(i) tabulate(codes[[i]], nbins[i]))
  p1 <- lapply(c1, flexsynth:::dp_normalise)
  joint <- function(i, j)
    as.vector(flexsynth:::dp_true_joint_array(codes, c(i, j), nbins))
  measured <- list(
    "1-2" = list(vars = c(1L, 2L), array = array(joint(1, 2), c(2, 2)), np = 1),
    "2-3" = list(vars = c(2L, 3L), array = array(joint(2, 3), c(2, 2)), np = 1))
  edges <- list(c(1L, 2L), c(2L, 3L))

  indep_ref <- flexsynth:::dp_aim_reference("independence", 3L, nbins, idx_all,
                                            c1, p1, n_est, measured, edges)
  model_ref <- flexsynth:::dp_aim_reference("model", 3L, nbins, idx_all,
                                            c1, p1, n_est, measured, edges)
  true_ac <- joint(1, 3)
  err_indep <- sum(abs(true_ac - indep_ref(1, 3)))
  err_model <- sum(abs(true_ac - model_ref(1, 3)))
  expect_lt(err_model, err_indep)
  # And on a directly-measured edge the model reference is near-exact (counts sum
  # to n_est ~ 6000; the reconciled model reproduces the measured edge to << 1).
  expect_lt(sum(abs(joint(1, 2) - model_ref(1, 2))), 1)
})

# ---- release shape ---------------------------------------------------------

test_that("model-scored Full AIM release is well-formed and reproducible", {
  df <- aim_chain_data(2500L, 5L)
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, scoring = "model")
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  expect_named(a, names(df))
  expect_false(anyNA(a))
  expect_identical(a, b)
})

test_that("model scoring also drives the annealed AIM schedule", {
  df <- aim_chain_data(2500L, 7L)
  sf <- 0.3
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, select_frac = sf,
                   anneal = TRUE, scoring = "model")
  res <- synth(df, structure = ~ id, privacy = dp, seed = 6)
  a <- as.data.frame(res)
  expect_false(anyNA(a))
  # Budget still composes exactly over the (data-adaptive) round count.
  am <- res$privacy$aim
  rho_total <- flexsynth:::zcdp_rho_for(8, 1e-6)
  expect_true(isTRUE(am$anneal))
  expect_identical(am$scoring, "model")
  expect_equal(sum(am$rho_meas_rounds) + sum(am$rho_sel_rounds),
               rho_total, tolerance = 1e-6)
})

# ---- privacy accounting invariance (post-processing) -----------------------

test_that("model scoring leaves the (eps, delta) accounting unchanged", {
  df <- aim_chain_data(2000L, 9L)
  mk <- function(sc) dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                                select = "aim", treewidth = 2, scoring = sc)
  ind <- synth(df, structure = ~ id, privacy = mk("independence"), seed = 2)$privacy
  mod <- synth(df, structure = ~ id, privacy = mk("model"), seed = 2)$privacy
  expect_equal(mod$epsilon, ind$epsilon)
  expect_equal(mod$delta, ind$delta)
  expect_equal(mod$rho, ind$rho, tolerance = 1e-12)
  # Fixed-round AIM measures a data-independent marginal count either way.
  expect_equal(mod$n_marginals, ind$n_marginals)
  expect_identical(mod$aim$scoring, "model")
  expect_identical(ind$aim$scoring, "independence")
})

test_that("default select = 'aim' output is identical to explicit independence", {
  df <- aim_chain_data(1500L, 11L)
  d1 <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2)
  d2 <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, scoring = "independence")
  a <- as.data.frame(synth(df, structure = ~ id, privacy = d1, seed = 3))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = d2, seed = 3))
  expect_identical(a, b)
})

# ---- prints ----------------------------------------------------------------

test_that("prints surface model-projection scoring", {
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, scoring = "model")
  expect_output(print(dp), "model")

  df <- aim_chain_data(1200L, 13L)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 1)$privacy
  expect_output(print(ac), "model")
})
