# Track B linked multi-table DP (synth_linked(privacy = dp_control(...))).

# A compact two-level hierarchy: patients (root) -> admissions (child).
lk_two <- function(np = 60, seed = 1) {
  set.seed(seed)
  patients <- data.frame(
    id  = seq_len(np),
    age = round(stats::rnorm(np, 60, 10)),
    sex = factor(sample(c("F", "M"), np, TRUE)),
    stringsAsFactors = FALSE)
  adm <- do.call(rbind, lapply(patients$id, function(pid) {
    n <- stats::rpois(1, 1.5)
    if (n == 0) return(NULL)
    data.frame(id = pid, admission_id = seq_len(n),
               los = 1L + stats::rpois(n, 3))
  }))
  list(tables = list(patients = patients, admissions = adm),
       structures = list(patients = ~ id, admissions = ~ id / admission_id),
       keys = list(patients = "id", admissions = c("id", "admission_id")))
}

pub2 <- list(age = c(20, 100), los = c(0, 40))

test_that("synth_linked dispatches to the DP engine and keeps referential integrity", {
  d <- lk_two()
  dp <- dp_control(epsilon = 8, mechanism = "laplace", dependence = "independent",
                   max_rows_per_person = 5, domain = "public", bounds = pub2)
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  expect_s3_class(res, "synth_linked_result")
  expect_s3_class(res$privacy, "dp_accounting")
  expect_false(is.null(res$privacy$linked))

  syn <- as.list(res)
  expect_setequal(names(syn), c("patients", "admissions"))
  expect_true(all(syn$admissions$id %in% syn$patients$id))   # no orphans
  # admission_id regenerated as a within-patient position (1..count, no gaps).
  by_pt <- split(syn$admissions$admission_id, syn$admissions$id)
  expect_true(all(vapply(by_pt, function(v) identical(sort(v), seq_along(v)),
                         logical(1))))
})

test_that("check_linkage passes on a DP linked result", {
  d <- lk_two()
  dp <- dp_control(epsilon = 6, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = 4, domain = "public", bounds = pub2)
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 3)
  cl <- check_linkage(res)
  expect_true(all(cl$duplicate_keys == FALSE))
  expect_true(all(cl$orphan_rows[!is.na(cl$orphan_rows)] == 0))
})

test_that("Laplace composition is exact over the whole release", {
  d <- lk_two()
  cap <- 5L; eps <- 8
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "independent",
                   max_rows_per_person = cap, domain = "public", bounds = pub2)
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  # independent + public bounds (marg_frac = 1):
  #   patients : 2 var * cap1                  = 2
  #   admissions: 1 var * cap + count sens 1   = cap + 1
  total_l1 <- 2 + cap + 1
  expect_equal(res$privacy$noise, total_l1 / eps)
  # histograms: 2 (patients) + 1 var + 1 count (admissions) = 4
  expect_equal(res$privacy$n_marginals, 4L)
})

test_that("Gaussian composition is exact (summed squared L2 over the release)", {
  d <- lk_two()
  cap <- 5L; eps <- 8; delta <- 1e-5
  dp <- dp_control(epsilon = eps, delta = delta, mechanism = "gaussian",
                   dependence = "independent",
                   max_rows_per_person = cap, domain = "public", bounds = pub2)
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  sum_sq <- 2 * 1 + (1 * cap^2 + 1)               # 3 + cap^2
  rho <- (sqrt(log(1 / delta) + eps) - sqrt(log(1 / delta)))^2
  expect_equal(res$privacy$noise, sqrt(sum_sq / (2 * rho)))
})

test_that("path caps and count sensitivities compound over three levels", {
  set.seed(7)
  np <- 40
  patients <- data.frame(id = seq_len(np), age = round(stats::rnorm(np, 50, 8)))
  adm <- do.call(rbind, lapply(patients$id, function(pid) {
    n <- 1L + stats::rpois(1, 1)
    data.frame(id = pid, admission_id = seq_len(n), los = 1L + stats::rpois(n, 2))
  }))
  labs <- do.call(rbind, lapply(seq_len(nrow(adm)), function(i) {
    n <- 1L + stats::rpois(1, 1)
    data.frame(id = adm$id[i], admission_id = adm$admission_id[i],
               lab_id = seq_len(n), value = round(stats::rnorm(n, 5, 1), 1))
  }))
  tables <- list(patients = patients, admissions = adm, labs = labs)
  structures <- list(patients = ~ id, admissions = ~ id / admission_id,
                     labs = ~ id / admission_id / lab_id)
  keys <- list(patients = "id", admissions = c("id", "admission_id"),
               labs = c("id", "admission_id", "lab_id"))
  eps <- 12
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "independent",
                   max_rows_per_person = list(admissions = 4, labs = 3),
                   domain = "public",
                   bounds = list(age = c(20, 90), los = c(0, 30), value = c(0, 10)))
  res <- synth_linked(tables, structures, keys, privacy = dp, seed = 1)

  info <- res$privacy$linked$tables
  by_name <- stats::setNames(info, vapply(info, function(x) x$name, character(1)))
  expect_equal(by_name$patients$path_cap, 1L)
  expect_equal(by_name$admissions$path_cap, 4L)
  expect_equal(by_name$labs$path_cap, 12L)          # 4 * 3
  expect_equal(by_name$admissions$count_sensitivity, 1L)  # path_cap[patients]
  expect_equal(by_name$labs$count_sensitivity, 4L)        # path_cap[admissions]

  # total L1 = patients(1var*1) + admissions(1var*4 + count1) + labs(1var*12 + count4)
  total_l1 <- (1 * 1) + (1 * 4 + 1) + (1 * 12 + 4)
  expect_equal(res$privacy$noise, total_l1 / eps)
  expect_true(all(res$syn$labs$id %in% res$syn$admissions$id))
})

test_that("domain = 'dp' estimates edges per table and composes with a budget slice", {
  d <- lk_two()
  cap <- 5L; eps <- 8
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "independent",
                   max_rows_per_person = cap, domain = "dp")  # no bounds -> estimate
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  total_l1 <- 2 + cap + 1
  expect_equal(res$privacy$noise, total_l1 / (0.9 * eps))     # marg_frac = 1 - 0.1
  expect_true(length(res$privacy$domain$vars) > 0)
  expect_true(any(grepl("age", res$privacy$domain$vars)))
})

test_that("missing branch cap, unit = row, and constraints are refused", {
  d <- lk_two()
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = dp_control(epsilon = 1, domain = "public", bounds = pub2),
                 seed = 1),
    "no branching cap")
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = dp_control(epsilon = 1, unit = "row",
                                      max_rows_per_person = 3,
                                      domain = "public", bounds = pub2),
                 seed = 1),
    "unit = \"row\"")
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 constraints = rule(los >= 0),
                 privacy = dp_control(epsilon = 1, max_rows_per_person = 3,
                                      domain = "public", bounds = pub2),
                 seed = 1),
    "constraints")
})

test_that("character columns and value-less roots are refused; deeper nesting too", {
  d <- lk_two()
  d$tables$patients$sex <- as.character(d$tables$patients$sex)   # bare character
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = dp_control(epsilon = 4, max_rows_per_person = 4,
                                      domain = "dp"),
                 seed = 1),
    "character")

  # A root with only its key: nothing to model.
  e <- lk_two()
  e$tables$patients <- e$tables$patients["id"]
  expect_error(
    synth_linked(e$tables, e$structures, e$keys,
                 privacy = dp_control(epsilon = 4, max_rows_per_person = 4,
                                      domain = "public", bounds = list(los = c(0, 40))),
                 seed = 1),
    "no modellable variable")

  # A within-table longitudinal index outside the key.
  f <- lk_two()
  f$tables$admissions$visit <- 1L
  f$structures$admissions <- ~ id / admission_id / visit
  expect_error(
    synth_linked(f$tables, f$structures, f$keys,
                 privacy = dp_control(epsilon = 4, max_rows_per_person = 4,
                                      domain = "public", bounds = pub2),
                 seed = 1),
    "outside its key")
})

test_that("a scalar cap applies to every child table", {
  d <- lk_two()
  dp <- dp_control(epsilon = 6, mechanism = "laplace", max_rows_per_person = 7,
                   domain = "public", bounds = pub2)
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  info <- res$privacy$linked$tables
  adm <- Filter(function(x) x$name == "admissions", info)[[1]]
  expect_equal(adm$local_cap, 7L)
})

test_that("m > 1 yields independent collections; a fixed seed reproduces", {
  d <- lk_two()
  dp <- dp_control(epsilon = 8, mechanism = "laplace", max_rows_per_person = 4,
                   domain = "public", bounds = pub2)
  res_m <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, m = 3, seed = 5)
  expect_length(res_m$syn, 3L)
  expect_s3_class(res_m$syn[[1]]$patients, "data.frame")

  a <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 11)
  b <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 11)
  expect_equal(as.list(a)$admissions, as.list(b)$admissions)
})

test_that("a named per-table cap is rejected by the flat synth() engine", {
  df <- data.frame(id = 1:30, y = stats::rnorm(30))
  expect_error(
    synth(df, structure = ~ id,
          privacy = dp_control(epsilon = 1, max_rows_per_person = c(admissions = 3)),
          seed = 1),
    "only for linked DP")
})

test_that("dp_control accepts and validates per-table caps", {
  expect_s3_class(dp_control(epsilon = 1, max_rows_per_person = list(a = 2, b = 5)),
                  "dp_control")
  expect_error(dp_control(epsilon = 1, max_rows_per_person = list(a = 0)),
               "positive integer")
  expect_error(dp_control(epsilon = 1, max_rows_per_person = list(a = 2.5)),
               "positive integer")
  # A single unnamed integer is still the scalar form.
  expect_identical(dp_control(epsilon = 1, max_rows_per_person = 4)$max_rows_per_person,
                   4L)
})
