# Phase 5: disclosure-risk diagnostics.

make_risk_df <- function(seed = 1, n = 200) {
  set.seed(seed)
  data.frame(
    id  = seq_len(n),
    age = round(rnorm(n, 60, 10)),
    sex = sample(c("F", "M"), n, replace = TRUE),
    sbp = rnorm(n, 130, 15),                     # continuous: exact pairs rare
    stringsAsFactors = FALSE
  )
}

test_that("copying the real data is maximally risky", {
  df <- make_risk_df()
  r  <- disclosure_risk(df, df, quasi = c("age", "sex"), seed = 1)
  expect_s3_class(r, "flexsynth_disclosure")
  expect_equal(r$replicated_uniques$prop_exact_syn, 1)   # every syn row copies
  expect_equal(r$dcr$prop_syn_zero, 1)                    # DCR all zero
  expect_output(print(r), "Replicated uniques")
})

test_that("genuinely synthesised data has lower identity risk than a copy", {
  df   <- make_risk_df()
  qi   <- c("age", "sex", "sbp")
  copy <- disclosure_risk(df, df, quasi = qi, seed = 1)
  res  <- synth(df, ~ id, seed = 3)
  r    <- disclosure_risk(df, res, quasi = qi, seed = 1)

  expect_equal(copy$replicated_uniques$prop_exact_syn, 1)  # a copy reproduces all
  expect_lt(r$replicated_uniques$prop_exact_syn, 1)        # synthesis: far fewer
  expect_lt(r$dcr$prop_syn_zero, copy$dcr$prop_syn_zero)   # fewer exact matches
  expect_gt(r$dcr$median_syn, 0)                           # and non-zero DCR
})

test_that("membership inference separates members from a distant holdout", {
  df   <- make_risk_df()
  hold <- data.frame(
    id  = 1:100,
    age = round(rnorm(100, 200, 5)),                      # far from training
    sex = sample(c("F", "M"), 100, replace = TRUE),
    stringsAsFactors = FALSE
  )
  r <- disclosure_risk(df, df, quasi = c("age", "sex"), holdout = hold, seed = 1)
  expect_false(is.null(r$membership))
  expect_gt(r$membership$auc, 0.9)                        # members are copied
  expect_gt(r$membership$advantage, 0.8)
})

test_that("membership inference is skipped without a holdout", {
  df <- make_risk_df()
  r  <- disclosure_risk(df, df, quasi = c("age", "sex"), seed = 1)
  expect_null(r$membership)
  expect_output(print(r), "not run")
})

test_that("unknown quasi-identifier columns error", {
  df <- make_risk_df()
  expect_error(disclosure_risk(df, df, quasi = c("age", "nope")),
               "not present")
})

test_that("disclosure_risk dispatches over a list of tables", {
  df <- make_risk_df()
  rl <- disclosure_risk(list(a = df, b = df), list(a = df, b = df),
                        quasi = c("age", "sex"), seed = 1)
  expect_s3_class(rl, "flexsynth_disclosure_list")
  expect_named(rl, c("a", "b"))
})
