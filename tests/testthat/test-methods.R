# Phase 4: method registry + built-in methods.

make_reg <- function(seed = 1, n = 300) {
  set.seed(seed)
  df <- data.frame(id = seq_len(n),
                   x = round(rnorm(n, 50, 10)),
                   g = sample(c("a", "b", "c"), n, replace = TRUE),
                   stringsAsFactors = FALSE)
  df$y <- round(2 * df$x + ifelse(df$g == "a", 0, 20) + rnorm(n, 0, 5))
  df
}

test_that("the registry lists the built-in methods", {
  m <- list_methods()
  expect_true(all(c("sample", "cart", "forest", "norm", "normrank", "ctree") %in% m))
})

test_that("register_method adds a method usable by synth()", {
  register_method(
    "const7",
    fit  = function(y, x, control) 7,
    draw = function(model, x, n, control) rep(model, n),
    categorical = FALSE, needs_predictors = FALSE
  )
  expect_true("const7" %in% list_methods())

  df <- make_reg()
  ctrl <- synth_control(method = stats::setNames(
    c("const7", "cart", "cart"), c("x", "g", "y")))
  s <- as.data.frame(synth(df, ~ id, tuning = ctrl, seed = 2))
  expect_true(all(s$x == 7))
})

test_that("unsupported methods are rejected", {
  df <- make_reg()
  expect_error(synth(df, ~ id, method = "no_such_method"),
               "unsupported method")
})

test_that("forest preserves a conditional relationship", {
  skip_if_not_installed("rpart")
  df <- make_reg()
  s <- as.data.frame(synth(df, ~ id, method = "forest",
                           tuning = synth_control(forest = list(ntree = 15)),
                           seed = 3))
  expect_gt(cor(s$x, s$y), 0.6)                 # ~0.87 in the real data
  expect_identical(vapply(s, class, ""), vapply(df, class, ""))
})

test_that("norm and normrank synthesise numeric variables in range", {
  df <- make_reg()
  for (mth in c("norm", "normrank")) {
    ctrl <- synth_control(method = stats::setNames(
      c(mth, "cart", mth), c("x", "g", "y")))
    s <- as.data.frame(synth(df, ~ id, tuning = ctrl, seed = 4))
    expect_gt(cor(s$x, s$y), 0.6)
    expect_gte(min(s$x), min(df$x))             # clamped / rank-mapped to range
    expect_lte(max(s$x), max(df$x))
  }
  # normrank returns only observed values (marginal preserved exactly)
  ctrl <- synth_control(method = stats::setNames(
    c("normrank", "cart", "normrank"), c("x", "g", "y")))
  s <- as.data.frame(synth(df, ~ id, tuning = ctrl, seed = 5))
  expect_true(all(s$x %in% df$x))
})

test_that("numeric-only methods reject categorical targets", {
  df <- make_reg()
  ctrl <- synth_control(method = c(g = "norm"))
  expect_error(synth(df, ~ id, tuning = ctrl, seed = 1),
               "cannot synthesise a categorical")
})

test_that("ctree works when partykit is available", {
  skip_if_not_installed("partykit")
  df <- make_reg()
  s <- as.data.frame(synth(df, ~ id, method = "ctree", seed = 6))
  expect_gt(cor(s$x, s$y), 0.6)
  expect_true(all(s$g %in% df$g))
})
