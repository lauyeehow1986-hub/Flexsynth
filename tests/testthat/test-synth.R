# Track A single-table engine (Phase 1).

make_long <- function(n_subj = 40, seed = 42) {
  set.seed(seed)
  age <- round(rnorm(n_subj, 62, 9))
  sex <- sample(c("F", "M"), n_subj, replace = TRUE)
  data.frame(
    id    = rep(seq_len(n_subj), each = 2),
    visit = rep(1:2, times = n_subj),
    age   = rep(age, each = 2),
    sex   = rep(sex, each = 2),
    sbp   = round(rnorm(n_subj * 2, 132, 16)),
    drug  = sample(c("A", "B", "C"), n_subj * 2, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("synth returns a synth_result with the same schema", {
  df <- make_long()
  res <- synth(df, structure = ~ id / visit, seed = 1)
  expect_s3_class(res, "synth_result")
  syn <- as.data.frame(res)
  expect_named(syn, names(df))
  expect_identical(vapply(syn, class, character(1)),
                   vapply(df, class, character(1)))
  expect_gt(nrow(syn), 0)
})

test_that("a fixed seed makes synthesis reproducible", {
  df <- make_long()
  a <- as.data.frame(synth(df, structure = ~ id / visit, seed = 123))
  b <- as.data.frame(synth(df, structure = ~ id / visit, seed = 123))
  expect_equal(a, b)
})

test_that("synthetic categories stay within the observed support", {
  df <- make_long()
  syn <- as.data.frame(synth(df, structure = ~ id / visit, seed = 7))
  expect_true(all(syn$sex %in% unique(df$sex)))
  expect_true(all(syn$drug %in% unique(df$drug)))
  # CART draws real donor values, so numeric output stays within the real range
  expect_gte(min(syn$sbp), min(df$sbp))
  expect_lte(max(syn$sbp), max(df$sbp))
})

test_that("subject-invariant columns stay constant within a synthetic unit", {
  df <- make_long()
  res <- synth(df, structure = ~ id / visit, seed = 11)
  expect_setequal(res$subject, c("age", "sex"))  # detected as unit-level
  syn <- as.data.frame(res)
  per_unit <- tapply(seq_len(nrow(syn)), syn$id, function(i) {
    length(unique(syn$age[i])) == 1L && length(unique(syn$sex[i])) == 1L
  })
  expect_true(all(per_unit))
  # a genuinely time-varying column is allowed to vary within a unit
  varies <- tapply(seq_len(nrow(syn)), syn$id,
                   function(i) length(unique(syn$sbp[i])))
  expect_true(any(varies > 1))
})

test_that("m > 1 yields a list of distinct datasets", {
  df <- make_long()
  res <- synth(df, structure = ~ id / visit, m = 3, seed = 5)
  expect_equal(res$m, 3L)
  expect_length(res$syn, 3L)
  expect_error(as.data.frame(res), "m > 1")
  expect_false(isTRUE(all.equal(res$syn[[1]], res$syn[[2]])))
})

test_that("method = 'sample' works without any tree fitting", {
  df <- make_long()
  syn <- as.data.frame(synth(df, structure = ~ id / visit,
                             method = "sample", seed = 2))
  expect_named(syn, names(df))
  expect_true(all(syn$drug %in% unique(df$drug)))
})

test_that("the flat (non-nested) case behaves like row synthesis", {
  set.seed(9)
  df <- data.frame(rid = 1:100, x = rnorm(100), g = sample(c("a", "b"), 100, TRUE),
                   stringsAsFactors = FALSE)
  syn <- as.data.frame(synth(df, structure = ~ rid, seed = 3))
  expect_equal(nrow(syn), 100)
  expect_named(syn, names(df))
  expect_true(all(syn$g %in% c("a", "b")))
})

test_that("visit_sequence reorders which variable is synthesised first", {
  df <- make_long()
  ctrl <- synth_control(visit_sequence = c("sbp", "age"))
  res <- synth(df, structure = ~ id / visit, tuning = ctrl, seed = 1)
  expect_equal(res$visit_sequence[1:2], c("sbp", "age"))
})

test_that("unsupported methods are rejected", {
  df <- make_long()
  expect_error(synth(df, structure = ~ id / visit, method = "neural"),
               "unsupported method")
})

test_that("an empty constraints list is a no-op", {
  df <- make_long()
  expect_silent(synth(df, structure = ~ id / visit, constraints = list(),
                      seed = 1))
})
