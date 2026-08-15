# Phase 6: data.table fast-paths and parallel replicate generation. The
# fast-paths must return results identical to their base-R fallbacks, and the
# serial (default) code path must be unchanged.

# ---- data.table fast-path: rbind_rows -------------------------------------

test_that("rbind_rows equals do.call(rbind) and preserves types", {
  parts <- lapply(1:20, function(i) {
    data.frame(id = i, x = rnorm(3), g = letters[1:3],
               flag = c(TRUE, FALSE, TRUE), stringsAsFactors = FALSE)
  })
  fast <- flexsynth:::rbind_rows(parts)
  base <- do.call(rbind, parts)
  rownames(fast) <- NULL
  rownames(base) <- NULL
  expect_equal(fast, base)
  expect_identical(vapply(fast, class, character(1)),
                   vapply(base, class, character(1)))
})

test_that("rbind_rows preserves factor columns and their levels", {
  lv <- c("A", "B", "C")
  parts <- list(
    data.frame(f = factor("A", levels = lv), y = 1),
    data.frame(f = factor(c("C", "B"), levels = lv), y = 2:3)
  )
  fast <- flexsynth:::rbind_rows(parts)
  expect_s3_class(fast$f, "factor")
  expect_identical(levels(fast$f), lv)
  expect_identical(as.character(fast$f), c("A", "C", "B"))
})

test_that("rbind_rows handles NULL / empty / single inputs", {
  expect_null(flexsynth:::rbind_rows(list()))
  expect_null(flexsynth:::rbind_rows(list(NULL, NULL)))
  one <- data.frame(a = 1:2)
  expect_identical(flexsynth:::rbind_rows(list(one)), one)
  expect_identical(flexsynth:::rbind_rows(list(NULL, one)), one)
})

# ---- data.table fast-path: constant_within --------------------------------

test_that("constant_within flags subject-invariant columns", {
  data <- data.frame(
    baseline = rep(c(10, 20, 30), each = 2),   # constant within id
    varying  = 1:6,                            # changes within id
    label    = rep(c("x", "y", "z"), each = 2) # constant within id
  )
  gid <- rep(1:3, each = 2)
  res <- flexsynth:::constant_within(data, gid, c("baseline", "varying", "label"))
  expect_identical(res, c(baseline = TRUE, varying = FALSE, label = TRUE))
})

test_that("constant_within matches an explicit tapply computation", {
  set.seed(11)
  gid <- rep(1:15, each = 3)
  data <- data.frame(
    a = rep(rnorm(15), each = 3),                 # constant
    b = rnorm(45),                                # varying
    c = rep(sample(c("p", "q"), 15, TRUE), each = 3)  # constant
  )
  cols <- c("a", "b", "c")
  fast <- flexsynth:::constant_within(data, gid, cols)
  base <- vapply(cols, function(col) {
    all(tapply(data[[col]], gid, function(v) length(unique(v)) == 1L))
  }, logical(1))
  expect_identical(fast, base)
})

# ---- parallel plumbing -----------------------------------------------------

test_that("parallel_workers maps the control to a worker count", {
  expect_identical(flexsynth:::parallel_workers(FALSE), 1L)
  expect_identical(flexsynth:::parallel_workers(NULL), 1L)
  expect_identical(flexsynth:::parallel_workers(3L), 3L)
  expect_identical(flexsynth:::parallel_workers(2), 2L)
  expect_gte(flexsynth:::parallel_workers(TRUE), 1L)
})

test_that("run_replicates serial equals a plain lapply and returns m items", {
  ctrl <- synth_control(parallel = FALSE)
  counter <- 0L
  gen <- function() { counter <<- counter + 1L; counter }
  out <- flexsynth:::run_replicates(3L, gen, ctrl, seed = NULL)
  expect_length(out, 3L)
  expect_identical(unlist(out), 1:3)
})

test_that("synth_control accepts integer parallel and rejects bad values", {
  expect_identical(synth_control(parallel = 4L)$parallel, 4L)
  expect_identical(synth_control(parallel = 2)$parallel, 2L)
  expect_true(isTRUE(synth_control(parallel = TRUE)$parallel))
  expect_false(synth_control(parallel = FALSE)$parallel)
  expect_error(synth_control(parallel = 0), "positive integer")
  expect_error(synth_control(parallel = -1), "positive integer")
  expect_error(synth_control(parallel = 1.5), "positive integer")
  expect_error(synth_control(parallel = "x"), "positive integer")
})

test_that("m = 1 never opens a cluster even when parallel is requested", {
  # workers is capped at m, so a single replicate always runs in-process.
  df <- data.frame(id = 1:30, x = rnorm(30))
  res <- synth(df, structure = ~ id,
               tuning = synth_control(parallel = 4L), m = 1, seed = 1)
  expect_s3_class(res, "synth_result")
  expect_equal(nrow(as.data.frame(res)), 30)
})

test_that("parallel synthesis returns m well-formed replicates", {
  # Guard: only run when a cluster can be started and flexsynth is loadable on
  # the workers (true under R CMD check, where the package is installed).
  cl <- tryCatch(parallel::makePSOCKcluster(2L), error = function(e) NULL)
  skip_if(is.null(cl), "no parallel cluster available")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  loadable <- tryCatch(
    all(unlist(parallel::clusterEvalQ(
      cl, requireNamespace("flexsynth", quietly = TRUE)))),
    error = function(e) FALSE)
  parallel::stopCluster(cl)
  on.exit(NULL)
  skip_if_not(isTRUE(loadable), "flexsynth not loadable on workers")

  df <- data.frame(id = rep(1:40, each = 2),
                   visit = rep(1:2, 40),
                   age = rep(round(rnorm(40, 60, 8)), each = 2),
                   sbp = round(rnorm(80, 130, 15)))
  res <- synth(df, structure = ~ id / visit,
               tuning = synth_control(parallel = 2L), m = 3, seed = 7)
  expect_identical(res$m, 3L)
  expect_length(res$syn, 3L)
  for (s in res$syn) {
    expect_named(s, names(df))
    expect_gt(nrow(s), 0)
  }
})

test_that("a parallel run is reproducible for a fixed seed and worker count", {
  can_par <- {
    cl <- tryCatch(parallel::makePSOCKcluster(2L), error = function(e) NULL)
    if (is.null(cl)) FALSE else {
      ok <- tryCatch(all(unlist(parallel::clusterEvalQ(
        cl, requireNamespace("flexsynth", quietly = TRUE)))),
        error = function(e) FALSE)
      parallel::stopCluster(cl)
      isTRUE(ok)
    }
  }
  skip_if_not(can_par, "no reproducible parallel backend")

  df <- data.frame(id = 1:60, x = rnorm(60), g = sample(c("a", "b"), 60, TRUE),
                   stringsAsFactors = FALSE)
  ctrl <- synth_control(parallel = 2L)
  a <- synth(df, structure = ~ id, tuning = ctrl, m = 4, seed = 99)$syn
  b <- synth(df, structure = ~ id, tuning = ctrl, m = 4, seed = 99)$syn
  expect_equal(a, b)
})
