# Phase 5: utility diagnostics.

make_diag_df <- function(seed = 1, n = 250) {
  set.seed(seed)
  df <- data.frame(
    id  = seq_len(n),
    age = round(rnorm(n, 60, 10)),
    sex = sample(c("F", "M"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  df$sbp <- round(1.2 * df$age + rnorm(n, 40, 8))
  df
}

test_that("diagnose reports a distance for every shared variable", {
  df <- make_diag_df()
  d  <- diagnose(df, df)                          # real vs itself
  expect_s3_class(d, "flexsynth_diagnostics")
  expect_setequal(d$univariate$variable, names(df))
  expect_true(all(d$univariate$distance >= 0 & d$univariate$distance <= 1,
                  na.rm = TRUE))
  # identical data: marginals match exactly
  expect_true(all(d$univariate$distance < 1e-8))
})

test_that("identical data is near-indistinguishable by pMSE", {
  df <- make_diag_df()
  d  <- diagnose(df, df)
  expect_lt(d$pmse$pmse, 1e-6)                    # propensities collapse to c
  expect_lt(d$correlation$frobenius, 1e-8)
})

test_that("a shifted synthetic set is flagged as distinguishable", {
  df  <- make_diag_df()
  bad <- df; bad$sbp <- bad$sbp + 60             # crude, systematically off
  d   <- diagnose(df, bad)
  expect_gt(max(d$univariate$distance), 0.3)     # KS picks up the shift
  expect_gt(d$pmse$ratio, 5)                      # very distinguishable
})

test_that("diagnose accepts a synth_result and honours vars", {
  df  <- make_diag_df()
  res <- synth(df, ~ id, seed = 2)
  d   <- diagnose(df, res, vars = c("age", "sbp"))
  expect_setequal(d$vars, c("age", "sbp"))
  expect_error(diagnose(df, res, vars = "nope"), "not present")
})

test_that("diagnose dispatches over a list of tables", {
  df <- make_diag_df()
  dl <- diagnose(list(a = df, b = df), list(a = df, b = df))
  expect_s3_class(dl, "flexsynth_diagnostics_list")
  expect_named(dl, c("a", "b"))
})

test_that("plot runs on a null device", {
  df <- make_diag_df()
  d  <- diagnose(df, df)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(d))
})
