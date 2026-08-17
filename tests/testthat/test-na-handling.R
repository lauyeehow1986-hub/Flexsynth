# Missing-data robustness for the parametric numeric methods. `norm`/`normrank`
# fit a linear model, which previously crashed when NA rows were present
# (model.matrix dropped them but the response kept full length). They must now
# fit on complete cases and synthesise without error.

make_na_data <- function() {
  set.seed(11)
  d <- data.frame(id = 1:300,
                  age = round(rnorm(300, 60, 10)),
                  sbp = round(rnorm(300, 130, 15)))
  d$sbp[sample(300, 90)] <- NA          # 30% missing in the target
  d$age[sample(300, 30)] <- NA          # some missing in a predictor
  d
}

test_that("method = 'norm' does not crash and preserves the missingness rate", {
  d <- make_na_data()
  expect_no_error(res <- synth(d, ~ id, method = "norm",
                               tuning = synth_control(smoothing = FALSE), seed = 1))
  s <- as.data.frame(res)
  expect_equal(nrow(s), nrow(d))
  # The missingness model preserves sbp's ~30% NA rate; observed values finite.
  expect_equal(mean(is.na(s$sbp)), mean(is.na(d$sbp)), tolerance = 0.08)
  expect_true(all(is.finite(s$sbp[!is.na(s$sbp)])))
})

test_that("method = 'normrank' does not crash and preserves the missingness rate", {
  d <- make_na_data()
  expect_no_error(res <- synth(d, ~ id, method = "normrank",
                               tuning = synth_control(smoothing = FALSE), seed = 1))
  s <- as.data.frame(res)
  expect_equal(nrow(s), nrow(d))
  expect_equal(mean(is.na(s$sbp)), mean(is.na(d$sbp)), tolerance = 0.08)
  # normrank maps to observed values, so every non-missing output is a real value.
  fin <- !is.na(s$sbp)
  expect_true(all(s$sbp[fin] %in% d$sbp[!is.na(d$sbp)]))
})

test_that("norm still recovers a linear signal with some missingness", {
  set.seed(12)
  x <- rnorm(600)
  d <- data.frame(id = 1:600, x = x, y = 1 + 2 * x + rnorm(600))
  d$y[sample(600, 120)] <- NA
  res <- synth(d, ~ id, method = "norm", m = 5, seed = 12)
  pooled <- synth_glm(res, y ~ x)
  slope <- pooled$estimates[pooled$estimates$term == "x", ]
  expect_true(slope$conf.low < 2 && slope$conf.high > 2)
})
