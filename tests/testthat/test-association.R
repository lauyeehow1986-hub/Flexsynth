# C9: categorical association (Cramer's V) view in diagnose().

test_that("diagnose reports a categorical association view", {
  set.seed(401)
  n <- 1500
  a <- sample(c("x", "y", "z"), n, TRUE)
  # b is strongly associated with a
  b <- ifelse(a == "x", sample(c("p", "q"), n, TRUE, c(0.9, 0.1)),
              sample(c("p", "q"), n, TRUE, c(0.2, 0.8)))
  d <- data.frame(id = seq_len(n), a = factor(a), b = factor(b))
  res <- synth(d, ~ id, seed = 1)
  dg <- diagnose(d, res)
  expect_false(is.null(dg$association))
  expect_true(dg$association$mean_abs_diff < 0.15)   # association roughly preserved
  expect_output(print(dg), "Cramer")
})

test_that("cramers_v is ~0 for independent factors and high for a copy", {
  set.seed(402)
  n <- 3000
  a <- factor(sample(letters[1:3], n, TRUE))
  indep <- factor(sample(letters[1:3], n, TRUE))
  expect_lt(cramers_v(a, indep), 0.1)
  expect_gt(cramers_v(a, a), 0.9)                    # perfect association
})

test_that("association view is absent with fewer than two categorical vars", {
  set.seed(403)
  d <- data.frame(id = 1:300, x = rnorm(300), g = factor(sample(c("a", "b"), 300, TRUE)))
  dg <- diagnose(d, as.data.frame(synth(d, ~ id, seed = 1)))
  expect_null(dg$association)
})
