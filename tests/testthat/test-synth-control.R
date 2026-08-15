test_that("synth_control returns a validated object with defaults", {
  ctrl <- synth_control()
  expect_s3_class(ctrl, "synth_control")
  expect_false(ctrl$proper)
  expect_false(ctrl$parallel)
  expect_null(ctrl$k)
})

test_that("synth_control validates its arguments", {
  expect_error(synth_control(proper = "yes"), "proper")
  expect_error(synth_control(k = 0), "k")
  expect_error(synth_control(parallel = NA), "parallel")
  expect_error(synth_control(cart = "not a list"), "cart")
})

test_that("synth_control has a print method", {
  expect_output(print(synth_control()), "synth_control")
})
