# B5: covariate-conditional rows-per-unit count model. When unit size depends on
# a subject covariate (e.g. sicker patients have more visits), the default
# marginal count model loses that link; count_model = "conditional" preserves it.

make_sizedep <- function(n_units = 400) {
  set.seed(9)
  sev <- sample(c("low", "high"), n_units, TRUE)
  sizes <- ifelse(sev == "high", 3L + rpois(n_units, 3), 1L + rpois(n_units, 0))
  do.call(rbind, lapply(seq_len(n_units), function(i) {
    m <- sizes[i]
    data.frame(id = i, visit = seq_len(m), sev = sev[i],
               sbp = round(rnorm(m, 130, 10)), stringsAsFactors = FALSE)
  }))
}

# mean unit size (rows per id) among high- vs low-severity synthetic units
size_gap <- function(s) {
  sz <- as.integer(table(s$id))
  sv <- tapply(s$sev, s$id, function(z) z[1])
  mean(sz[sv == "high"]) - mean(sz[sv == "low"])
}

test_that("conditional count model preserves size-covariate dependence", {
  d <- make_sizedep()
  real_gap <- size_gap(d)
  expect_gt(real_gap, 2)                                   # sanity: real signal
  s <- as.data.frame(synth(d, ~ id / visit, method = "cart",
                           tuning = synth_control(count_model = "conditional"),
                           seed = 1))
  expect_gt(size_gap(s), 1.5)                              # dependence survives
})

test_that("marginal count model (default) loses the size-covariate dependence", {
  d <- make_sizedep()
  s <- as.data.frame(synth(d, ~ id / visit, method = "cart", seed = 1))
  expect_lt(abs(size_gap(s)), 1.2)                         # size ~ independent of sev
})

test_that("count_model is validated and defaults to marginal", {
  expect_identical(synth_control()$count_model, "marginal")
  expect_s3_class(synth_control(count_model = "conditional"), "synth_control")
  expect_error(synth_control(count_model = "bogus"), "count_model|arg")
})
