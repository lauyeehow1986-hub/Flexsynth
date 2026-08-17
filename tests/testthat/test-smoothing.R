# Default smoothing: continuous numeric variables (many distinct values) are
# density-smoothed by default so CART draws are not confined to observed values
# (better marginal realism, less exact-value replication). Low-cardinality
# numeric codes are left alone, and smoothing can be turned off explicitly.

test_that("a continuous variable is smoothed by default (new values appear)", {
  set.seed(201)
  d <- data.frame(id = 1:500, x = rnorm(500), sbp = rnorm(500, 130, 15))
  s <- as.data.frame(synth(d, ~ id, method = "cart", seed = 1))
  # smoothing adds continuous noise -> synthetic values are not a subset of real
  expect_false(all(s$sbp %in% d$sbp))
  # marginal is still preserved
  expect_equal(mean(s$sbp), mean(d$sbp), tolerance = 3)
  expect_equal(sd(s$sbp),  sd(d$sbp),   tolerance = 4)
})

test_that("a low-cardinality numeric code is NOT smoothed by default", {
  set.seed(202)
  d <- data.frame(id = 1:500,
                  grade = sample(1:4, 500, TRUE),          # 4 distinct -> categorical
                  x = rnorm(500))
  s <- as.data.frame(synth(d, ~ id, method = "cart", seed = 1))
  expect_true(all(s$grade %in% d$grade))                   # stays on the observed set
})

test_that("smoothing = FALSE disables smoothing (values stay observed)", {
  set.seed(203)
  d <- data.frame(id = 1:500, sbp = round(rnorm(500, 130, 15) * 10))  # many distinct
  s <- as.data.frame(synth(d, ~ id, method = "cart",
                           tuning = synth_control(smoothing = FALSE), seed = 1))
  expect_true(all(s$sbp %in% d$sbp))
})

test_that("smoothing = FALSE is accepted by synth_control", {
  expect_s3_class(synth_control(smoothing = FALSE), "synth_control")
  expect_s3_class(synth_control(smoothing = "none"), "synth_control")
})
