# A missingness model: variables with NAs get a companion indicator synthesised
# in sequence; the value model is fit on observed rows; missingness is imposed at
# the end and must not cascade into fully-observed downstream variables.

test_that("missingness rate is preserved without a warning (cart)", {
  set.seed(101)
  d <- data.frame(id = 1:400, age = round(rnorm(400, 60, 10)),
                  sbp = round(rnorm(400, 130, 15)))
  d$sbp[sample(400, 120)] <- NA                 # 30% missing
  w <- 0
  s <- withCallingHandlers(
    as.data.frame(synth(d, ~ id, method = "cart", seed = 1)),
    warning = function(x) { w <<- w + 1; invokeRestart("muffleWarning") })
  expect_equal(w, 0)                             # no more "split variable" warning
  expect_equal(mean(is.na(s$sbp)), 0.30, tolerance = 0.06)
})

test_that("no companion indicator columns leak into the output", {
  set.seed(102)
  d <- data.frame(id = 1:200, x = rnorm(200), y = rnorm(200))
  d$y[sample(200, 40)] <- NA
  s <- as.data.frame(synth(d, ~ id, seed = 1))
  expect_equal(sort(names(s)), sort(names(d)))
  expect_false(any(grepl("^\\.na_", names(s))))
})

test_that("missingness does not cascade into fully-observed variables", {
  set.seed(103)
  d <- data.frame(id = 1:400, a = rnorm(400), b = rnorm(400), c = rnorm(400))
  d$a[sample(400, 160)] <- NA                   # only `a` is missing
  s <- as.data.frame(synth(d, ~ id, method = "cart", seed = 1))
  expect_false(anyNA(s$b))                       # b, c had no NAs -> none introduced
  expect_false(anyNA(s$c))
  expect_true(mean(is.na(s$a)) > 0.2)            # a's missingness preserved
})

test_that("observed value distribution is unbiased by the placeholder", {
  set.seed(104)
  d <- data.frame(id = 1:600, age = round(rnorm(600, 60, 10)),
                  sbp = round(rnorm(600, 140, 12)))
  d$sbp[sample(600, 200)] <- NA
  s <- as.data.frame(synth(d, ~ id, method = "cart", seed = 1))
  # mean of the observed synthetic values should track the real observed mean,
  # i.e. the imputed placeholder did not leak into the value model.
  expect_equal(mean(s$sbp, na.rm = TRUE), mean(d$sbp, na.rm = TRUE),
               tolerance = 4)
})

test_that("missingness conditional on a predictor is captured", {
  set.seed(105)
  n <- 1500
  age <- round(rnorm(n, 60, 12))
  sbp <- round(rnorm(n, 130, 15))
  # sbp is missing much more often for older subjects (MAR on age)
  p_miss <- plogis((age - 60) / 6)
  sbp[runif(n) < p_miss] <- NA
  d <- data.frame(id = seq_len(n), age = age, sbp = sbp)
  s <- as.data.frame(synth(d, ~ id, method = "cart", seed = 1))
  # the real association: mean age is higher among missing rows
  real_gap <- mean(d$age[is.na(d$sbp)]) - mean(d$age[!is.na(d$sbp)])
  syn_gap  <- mean(s$age[is.na(s$sbp)]) - mean(s$age[!is.na(s$sbp)])
  expect_gt(real_gap, 1)                         # sanity: signal exists
  expect_gt(syn_gap, 0.5)                         # and it survives synthesis
})

test_that("norm preserves missingness rate without cascade", {
  set.seed(106)
  d <- data.frame(id = 1:500, x = rnorm(500), y = rnorm(500))
  d$y[sample(500, 150)] <- NA                   # only y missing; x fully observed
  s <- as.data.frame(synth(d, ~ id, method = "norm", seed = 1))
  expect_false(anyNA(s$x))
  expect_equal(mean(is.na(s$y)), 0.30, tolerance = 0.07)
})
