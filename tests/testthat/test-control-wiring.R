# Phase 4: synth_control() knobs wired through the engine.

make_ctrl_df <- function(seed = 1, n = 300) {
  set.seed(seed)
  df <- data.frame(id = seq_len(n), x = as.integer(round(rnorm(n, 50, 10))))
  df$y <- round(1.5 * df$x + rnorm(n, 0, 4))
  df
}

test_that("smoothing jitters numeric draws and preserves integer type", {
  skip_if_not_installed("rpart")
  df <- make_ctrl_df()
  free   <- as.data.frame(synth(df, ~ id, method = "cart",
                                tuning = synth_control(smoothing = FALSE), seed = 5))
  smooth <- as.data.frame(synth(df, ~ id, method = "cart",
                                tuning = synth_control(smoothing = "x"), seed = 5))
  expect_gt(length(unique(smooth$x)), length(unique(free$x)))
  expect_true(is.integer(smooth$x))                    # integer target stays integer
  expect_gte(min(smooth$x), min(df$x))                 # clamped to observed range
  expect_lte(max(smooth$x), max(df$x))
})

test_that("smoothing = TRUE smooths every numeric column", {
  skip_if_not_installed("rpart")
  df <- make_ctrl_df()
  free   <- as.data.frame(synth(df, ~ id, method = "cart",
                                tuning = synth_control(smoothing = FALSE), seed = 9))
  smooth <- as.data.frame(synth(df, ~ id, method = "cart",
                                tuning = synth_control(smoothing = TRUE), seed = 9))
  expect_gt(length(unique(smooth$y)), length(unique(free$y)))
})

test_that("predictor_matrix severs a conditioning path", {
  skip_if_not_installed("rpart")
  df <- make_ctrl_df()
  base <- as.data.frame(synth(df, ~ id, method = "cart", seed = 2))
  expect_gt(cor(base$x, base$y), 0.6)                  # x drives y by default

  pm <- matrix(1, 2, 2, dimnames = list(c("x", "y"), c("x", "y")))
  pm["y", "x"] <- 0                                    # forbid x -> y
  cut <- as.data.frame(synth(df, ~ id, method = "cart",
                             tuning = synth_control(predictor_matrix = pm), seed = 2))
  expect_lt(abs(cor(cut$x, cut$y)), 0.3)               # relationship gone
})

test_that("forest honours its ntree hyperparameter", {
  skip_if_not_installed("rpart")
  df <- make_ctrl_df()
  s <- as.data.frame(synth(df, ~ id, method = "forest",
                           tuning = synth_control(forest = list(ntree = 3, mtry = 1)),
                           seed = 4))
  expect_equal(nrow(s), nrow(df))
  expect_identical(vapply(s, class, ""), vapply(df, class, ""))
})

test_that("synth_control validates the new knobs", {
  expect_error(synth_control(smoothing = 1L), "smoothing")
  expect_error(synth_control(constraint_max_tries = 0), "constraint_max_tries")
  ctrl <- synth_control(smoothing = "x", constraint_max_tries = 12)
  expect_equal(ctrl$constraint_max_tries, 12L)
  expect_output(print(ctrl), "constraint tries")
})
