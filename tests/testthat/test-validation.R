test_that("synth validates inputs and signals not-yet-implemented", {
  df <- data.frame(id = 1:4, visit = c(1, 2, 1, 2), y = rnorm(4))
  expect_error(synth(df, structure = ~ id / visit), "not implemented yet")
  # Track B path is validated too
  expect_error(
    synth(df, structure = ~ id / visit, privacy = dp_control(epsilon = 1)),
    "Track B"
  )
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
