test_that("synth runs Track A, and Track B routes to the DP engine", {
  df <- data.frame(id = 1:40, visit = 1, y = rnorm(40))
  expect_s3_class(synth(df, structure = ~ id / visit, seed = 1), "synth_result")
  # Track B (DP) is now built: it returns a differentially private result.
  res <- suppressWarnings(
    synth(df, structure = ~ id / visit, privacy = dp_control(epsilon = 1), seed = 1))
  expect_s3_class(res, "synth_result")
  expect_s3_class(res$privacy, "dp_accounting")
})

test_that("synth rejects structure variables absent from data", {
  df <- data.frame(id = 1:4, y = rnorm(4))
  expect_error(synth(df, structure = ~ id / visit), "not found in `data`")
})

test_that("synth rejects bad structure and data", {
  df <- data.frame(id = 1:4, visit = c(1, 2, 1, 2))
  expect_error(synth("not a df", structure = ~ id), "data.frame")
  expect_error(synth(df, structure = "id / visit"), "one-sided formula")
  expect_error(synth(df, structure = ~ id, m = 0), "positive integer")
})

test_that("synth_linked requires aligned named lists", {
  adm <- data.frame(id = 1:2, admission_id = 1:2)
  expect_error(
    synth_linked(tables = list(adm), structures = list(), keys = list()),
    "named"
  )
  expect_error(
    synth_linked(
      tables = list(admissions = adm),
      structures = list(admissions = ~ id / admission_id),
      keys = list(wrong = "id")
    ),
    "keys"
  )
})
