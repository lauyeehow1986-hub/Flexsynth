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

test_that("dp_control defaults to a tree model with numeric binning", {
  dp <- dp_control(epsilon = 1)
  expect_equal(dp$dependence, "tree")
  expect_equal(dp$bins, 12L)
  expect_null(dp$max_rows_per_person)
})

test_that("dp_control validates the new marginal-model parameters", {
  expect_error(dp_control(epsilon = 1, bins = 1), "bins")
  expect_error(dp_control(epsilon = 1, bins = 2.5), "bins")
  expect_error(dp_control(epsilon = 1, max_rows_per_person = 0), "max_rows_per_person")
  expect_error(dp_control(epsilon = 1, max_rows_per_person = 2.5), "max_rows_per_person")
  expect_error(dp_control(epsilon = 1, bounds = list(c(0, 1))), "named list")
  expect_error(dp_control(epsilon = 1, bounds = list(a = c(1, 0))), "lower < upper")
  expect_error(dp_control(epsilon = 1, dependence = "bogus"))
})

test_that("dp_control accepts public bounds and an explicit row cap", {
  dp <- dp_control(epsilon = 1, max_rows_per_person = 3,
                   bounds = list(age = c(0, 120)))
  expect_identical(dp$max_rows_per_person, 3L)
  expect_identical(dp$bounds$age, c(0, 120))
})
