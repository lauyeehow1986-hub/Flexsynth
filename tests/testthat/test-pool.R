# Tests for the valid-inference pooling layer (pool_synth / synth_glm).
# The arithmetic tests independently reproduce the synthpop (Raab/Nowok/Dibben)
# and Reiter (2003) fully-synthetic combining rules, so a regression in the
# formulae is caught without depending on synthpop being installed.

# --- pure combining arithmetic -------------------------------------------

test_that("synthpop rule, simple synthesis, population inference matches formula", {
  ests <- list(c(a = 2), c(a = 4), c(a = 6))      # qbar = 4, m = 3
  vars <- list(c(a = 1), c(a = 1), c(a = 1))      # ubar = 1
  out <- combine_estimates(ests, vars, n = 100, k = 100, proper = FALSE,
                           rule = "synthpop", population_inference = TRUE,
                           conf.level = 0.95)
  # vars_scaled = ubar * k/n = 1; Tf = vars_scaled * (1 + n/k/m) = 1 + 1/3
  expect_equal(out$estimate, 4)
  expect_equal(out$std.error, sqrt(1 + 1/3))
})

test_that("synthpop rule, proper synthesis uses (n/k + 1)/m", {
  ests <- list(c(a = 2), c(a = 4), c(a = 6))
  vars <- list(c(a = 1), c(a = 1), c(a = 1))
  out <- combine_estimates(ests, vars, n = 100, k = 100, proper = TRUE,
                           rule = "synthpop", population_inference = TRUE,
                           conf.level = 0.95)
  # Tf = 1 * (1 + (n/k + 1)/m) = 1 + 2/3
  expect_equal(out$std.error, sqrt(1 + 2/3))
})

test_that("synthpop rule scales within-variance by k/n", {
  ests <- list(c(a = 2), c(a = 4), c(a = 6))
  vars <- list(c(a = 1), c(a = 1), c(a = 1))
  out <- combine_estimates(ests, vars, n = 200, k = 100, proper = FALSE,
                           rule = "synthpop", population_inference = TRUE,
                           conf.level = 0.95)
  # vars_scaled = 1 * 100/200 = 0.5; Tf = 0.5 * (1 + (200/100)/3) = 0.5 * (1 + 2/3)
  expect_equal(out$std.error, sqrt(0.5 * (1 + 2/3)))
})

test_that("population_inference = FALSE gives the rescaled within-variance", {
  ests <- list(c(a = 2), c(a = 4), c(a = 6))
  vars <- list(c(a = 2), c(a = 2), c(a = 2))      # ubar = 2
  out <- combine_estimates(ests, vars, n = 100, k = 100, proper = FALSE,
                           rule = "synthpop", population_inference = FALSE,
                           conf.level = 0.95)
  expect_equal(out$std.error, sqrt(2))            # sqrt(ubar * k/n)
})

test_that("reiter rule matches (1 + 1/m) b_m - ubar", {
  ests <- list(c(a = 2), c(a = 4), c(a = 6))      # bm = var(c(2,4,6)) = 4
  vars <- list(c(a = 1), c(a = 1), c(a = 1))      # ubar = 1
  out <- combine_estimates(ests, vars, n = 100, k = 100, proper = FALSE,
                           rule = "reiter", population_inference = TRUE,
                           conf.level = 0.95)
  # Tf = (1 + 1/3) * 4 - 1
  expect_equal(out$std.error, sqrt((1 + 1/3) * 4 - 1))
})

test_that("reiter negative variance warns and returns NA se", {
  ests <- list(c(a = 3), c(a = 3), c(a = 3))      # bm = 0
  vars <- list(c(a = 1), c(a = 1), c(a = 1))      # Tf = 0 - 1 < 0
  expect_warning(
    out <- combine_estimates(ests, vars, n = 100, k = 100, proper = FALSE,
                             rule = "reiter", population_inference = TRUE,
                             conf.level = 0.95),
    "negative")
  expect_true(is.na(out$std.error))
})

# --- public API integration ----------------------------------------------

test_that("synth_glm recovers a known linear coefficient", {
  set.seed(1)
  n <- 800
  x <- rnorm(n)
  d <- data.frame(id = seq_len(n), x = x, y = 2 + 3 * x + rnorm(n))
  res <- synth(d, ~ id, m = 5, seed = 1)
  pooled <- synth_glm(res, y ~ x)
  slope <- pooled$estimates[pooled$estimates$term == "x", ]
  # the true slope (3) should sit inside the pooled 95% CI
  expect_true(slope$conf.low < 3 && slope$conf.high > 3)
  expect_s3_class(pooled, "flexsynth_pool")
})

test_that("pool_synth accepts a custom analysis returning estimate/variance", {
  set.seed(2)
  d <- data.frame(id = 1:500, y = rnorm(500, 10, 2))
  res <- synth(d, ~ id, m = 4, seed = 2)
  mean_analysis <- function(dat)
    list(estimate = c(mu = mean(dat$y)),
         variance = c(mu = stats::var(dat$y) / nrow(dat)))
  pooled <- pool_synth(res, mean_analysis)
  expect_equal(nrow(pooled$estimates), 1L)
  expect_true(is.finite(pooled$estimates$std.error))
})

test_that("pool_synth refuses DP (Track B) results", {
  set.seed(3)
  d <- data.frame(id = 1:400, y = rnorm(400))
  res <- synth(d, ~ id, privacy = dp_control(epsilon = 1, bounds = list(y = c(-4, 4))),
               m = 3, seed = 3)
  expect_error(pool_synth(res, function(dat) stats::lm(y ~ 1, dat)),
               "Track B|differentially private|DP")
})

test_that("reiter requires m >= 2; synthpop allows m = 1", {
  set.seed(4)
  d <- data.frame(id = 1:300, x = rnorm(300), y = rnorm(300))
  res1 <- synth(d, ~ id, m = 1, seed = 4)
  expect_error(synth_glm(res1, y ~ x, rule = "reiter"), "m >= 2|at least 2")
  expect_s3_class(synth_glm(res1, y ~ x, rule = "synthpop"), "flexsynth_pool")
})

test_that("pool result prints", {
  set.seed(5)
  d <- data.frame(id = 1:300, x = rnorm(300), y = rnorm(300))
  res <- synth(d, ~ id, m = 3, seed = 5)
  expect_output(print(synth_glm(res, y ~ x)), "flexsynth_pool|estimate|Beta")
})
