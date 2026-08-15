# Phase 4: constraint / temporal-logic system.

test_that("rule() captures the expression, scope and label", {
  r <- rule(dbp <= sbp)
  expect_s3_class(r, "flexsynth_rule")
  expect_equal(r$scope, "row")
  expect_match(r$label, "dbp")
  u <- rule(all(diff(los) >= 0), scope = "unit")
  expect_equal(u$scope, "unit")
  expect_output(print(r), "flexsynth_rule")
})

test_that("a row constraint eliminates violations that otherwise occur", {
  set.seed(1); n <- 400
  df <- data.frame(id = seq_len(n),
                   a = round(rnorm(n, 50, 20)),
                   b = round(rnorm(n, 50, 20)))
  expect_gt(sum(df$a > df$b), 0)                      # real data violates

  free <- as.data.frame(synth(df, ~ id, method = "sample", seed = 7))
  expect_gt(sum(free$a > free$b), 0)                  # so does unconstrained synth

  con <- as.data.frame(synth(df, ~ id, method = "sample",
                             constraints = rule(a <= b), seed = 7))
  expect_equal(sum(con$a > con$b), 0)                 # constraint removes them
  expect_equal(nrow(con), nrow(df))                   # without dropping rows
})

test_that("a within-unit temporal constraint is enforced", {
  skip_if_not_installed("rpart")
  set.seed(2)
  adm <- do.call(rbind, lapply(1:120, function(i) {
    k <- 2L + rpois(1, 1)
    data.frame(id = i, visit = seq_len(k), los = sort(1L + rpois(k, 3)))
  }))
  res <- synth(adm, ~ id / visit,
               constraints = rule(all(diff(los) >= 0), scope = "unit"), seed = 3)
  s <- as.data.frame(res)
  decreasing <- tapply(s$los, s$id, function(v) any(diff(v) < 0))
  expect_false(any(decreasing))
})

test_that("multiple rules must all hold", {
  set.seed(4); n <- 300
  df <- data.frame(id = seq_len(n),
                   a = round(rnorm(n, 50, 15)),
                   b = round(rnorm(n, 50, 15)))
  s <- as.data.frame(synth(df, ~ id, method = "sample",
                           constraints = list(rule(a <= b), rule(b <= 90)),
                           seed = 8))
  expect_true(all(s$a <= s$b))
  expect_true(all(s$b <= 90))
})

test_that("an unsatisfiable constraint warns and returns the valid units only", {
  set.seed(5); n <- 100
  df <- data.frame(id = seq_len(n), a = round(rnorm(n, 50, 10)))
  expect_warning(
    res <- synth(df, ~ id, method = "sample",
                 constraints = rule(a > 1e6),
                 tuning = synth_control(constraint_max_tries = 3), seed = 1),
    "kept 0 of"
  )
  expect_equal(nrow(as.data.frame(res)), 0L)
})

test_that("constraints referencing unknown columns are rejected", {
  df <- data.frame(id = 1:10, a = 1:10)
  expect_error(synth(df, ~ id, constraints = rule(nope <= a)),
               "unknown column")
})

test_that("non-rule constraints are rejected", {
  df <- data.frame(id = 1:10, a = 1:10)
  expect_error(synth(df, ~ id, constraints = "a <= 5"),
               "must be a rule")
})
