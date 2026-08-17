# C8: specific-utility comparison. Fit the same analysis on real and synthetic
# data and compare the estimates via confidence-interval overlap (Karr et al.)
# and the standardised difference -- "does my regression come out the same?".

test_that("high overlap when synthetic reproduces the real relationship", {
  set.seed(301)
  n <- 1500
  x <- rnorm(n)
  d <- data.frame(id = seq_len(n), x = x, y = 1 + 2 * x + rnorm(n))
  res <- synth(d, ~ id, m = 10, seed = 1)
  cmp <- compare_estimates(d, res, function(dat) lm(y ~ x, dat))
  slope <- cmp$estimates[cmp$estimates$term == "x", ]
  expect_gt(slope$overlap, 0.7)                 # CIs largely overlap
  expect_lt(slope$std_diff, 0.5)                # estimates close on the real SE scale
  expect_s3_class(cmp, "flexsynth_utility")
})

test_that("compare_estimates works on a single synthetic data.frame too", {
  set.seed(302)
  d <- data.frame(id = 1:800, x = rnorm(800))
  d$y <- 0.5 + 1.5 * d$x + rnorm(800)
  syn <- as.data.frame(synth(d, ~ id, seed = 2))
  cmp <- compare_estimates(d, syn, function(dat) lm(y ~ x, dat))
  expect_true(all(c("term", "est_real", "est_syn", "overlap", "std_diff")
                  %in% names(cmp$estimates)))
  expect_output(print(cmp), "overlap|flexsynth_utility")
})

test_that("overlap is low when synthetic breaks the relationship", {
  set.seed(303)
  n <- 1200
  d <- data.frame(id = seq_len(n), x = rnorm(n), y = rnorm(n))
  d$y <- 3 * d$x + rnorm(n)                      # strong real slope ~ 3
  # synthesise x and y independently (destroys the x-y link)
  syn <- data.frame(x = rnorm(n), y = 3 * 0 + rnorm(n))
  cmp <- compare_estimates(d, syn, function(dat) lm(y ~ x, dat))
  slope <- cmp$estimates[cmp$estimates$term == "x", ]
  expect_lt(slope$overlap, 0.5)
})
