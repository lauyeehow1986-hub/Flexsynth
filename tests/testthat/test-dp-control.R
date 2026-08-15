test_that("dp_control builds a person-level DP spec by default", {
  dp <- dp_control(epsilon = 1)
  expect_s3_class(dp, "dp_control")
  expect_equal(dp$unit, "person")
  expect_equal(dp$mechanism, "laplace")
})

test_that("dp_control validates the privacy budget", {
  expect_error(dp_control(), "epsilon")
  expect_error(dp_control(epsilon = 0), "epsilon")
  expect_error(dp_control(epsilon = 1, delta = 1), "delta")
  expect_error(dp_control(epsilon = 1, delta = 0, mechanism = "gaussian"), "gaussian")
})

test_that("gaussian mechanism accepts delta > 0", {
  dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian")
  expect_equal(dp$mechanism, "gaussian")
})

test_that("dp_control has a print method", {
  expect_output(print(dp_control(epsilon = 1)), "Track B")
})
