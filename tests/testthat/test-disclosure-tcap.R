# TCAP (Target Correct Attribution Probability): given quasi-identifier keys, how
# well can an attacker read a sensitive target off the synthetic data, versus a
# marginal-only baseline. High TCAP well above baseline = attribute-disclosure
# risk; TCAP near baseline = the synthetic data leaks no more than the margin.

make_attr_data <- function(n = 1200, leak = TRUE) {
  set.seed(7)
  region <- sample(c("N", "S", "E", "W"), n, TRUE)
  ageband <- sample(c("young", "mid", "old"), n, TRUE)
  # sensitive target strongly determined by the key combination when leak=TRUE
  base <- interaction(region, ageband, drop = TRUE)
  dx <- if (leak) ifelse(as.integer(base) %% 2 == 0, "yes", "no")
        else sample(c("yes", "no"), n, TRUE)
  data.frame(id = seq_len(n), region = region, ageband = ageband,
             dx = dx, stringsAsFactors = FALSE)
}

test_that("an exact copy yields high TCAP well above baseline", {
  d <- make_attr_data(leak = TRUE)
  # syn = real is the worst case for attribute disclosure
  r <- disclosure_risk(d, d, quasi = c("region", "ageband", "dx"),
                       target = "dx", seed = 1)
  a <- r$attribute
  expect_false(is.null(a))
  expect_gt(a$tcap, 0.9)
  expect_gt(a$tcap - a$baseline, 0.2)          # real lift over the marginal attacker
  expect_true(a$coverage > 0.9)
})

test_that("a target independent of the keys gives TCAP near baseline", {
  d <- make_attr_data(leak = FALSE)             # dx independent of region/ageband
  r <- disclosure_risk(d, d, quasi = c("region", "ageband", "dx"),
                       target = "dx", seed = 1)
  a <- r$attribute
  expect_lt(abs(a$tcap - a$baseline), 0.1)      # no more than the margin reveals
})

test_that("attribute disclosure is skipped unless a target is given", {
  d <- make_attr_data()
  r <- disclosure_risk(d, d, quasi = c("region", "ageband"), seed = 1)
  expect_null(r$attribute)
  expect_output(print(r), "not run|Attribute")
})

test_that("target must be present and leave at least one key", {
  d <- make_attr_data()
  expect_error(disclosure_risk(d, d, quasi = c("region", "ageband"),
                               target = "nope"), "target")
  expect_error(disclosure_risk(d, d, quasi = "dx", target = "dx"), "key")
})
