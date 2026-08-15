# Phase 3: linked multi-table joint synthesis.

# A three-level synthetic hierarchy: patients -> admissions -> procedures.
make_linked <- function(seed = 1) {
  set.seed(seed)
  np <- 40
  patients <- data.frame(
    id  = seq_len(np),
    sex = sample(c("F", "M"), np, replace = TRUE),
    age = round(rnorm(np, 62, 10)),
    stringsAsFactors = FALSE
  )
  adm <- do.call(rbind, lapply(patients$id, function(pid) {
    na <- 1L + rpois(1, 0.7)
    data.frame(id = pid, admission_id = seq_len(na),
               los = 1L + rpois(na, 3), stringsAsFactors = FALSE)
  }))
  proc <- do.call(rbind, lapply(seq_len(nrow(adm)), function(i) {
    n <- rpois(1, 1.0)
    if (n == 0L) return(NULL)
    data.frame(id = adm$id[i], admission_id = adm$admission_id[i],
               procedure_number = seq_len(n),
               kind = sample(c("PCI", "angio", "echo"), n, replace = TRUE),
               stringsAsFactors = FALSE)
  }))
  list(patients = patients, admissions = adm, procedures = proc)
}

linked_structures <- list(
  patients   = ~ id,
  admissions = ~ id / admission_id,
  procedures = ~ id / admission_id / procedure_number
)
linked_keys <- list(
  patients   = "id",
  admissions = c("id", "admission_id"),
  procedures = c("id", "admission_id", "procedure_number")
)

test_that("linked synthesis preserves referential integrity and schema", {
  skip_if_not_installed("rpart")
  tabs <- make_linked()
  res <- synth_linked(tabs, linked_structures, linked_keys,
                      method = "cart", seed = 2)
  s <- res$syn

  # every child foreign key resolves to an existing parent key
  expect_true(all(s$admissions$id %in% s$patients$id))
  adm_key  <- paste(s$admissions$id, s$admissions$admission_id)
  proc_pk  <- paste(s$procedures$id, s$procedures$admission_id)
  expect_true(all(proc_pk %in% adm_key))

  # keys identify rows uniquely within each table
  expect_equal(anyDuplicated(paste(s$patients$id)), 0L)
  expect_equal(anyDuplicated(adm_key), 0L)
  expect_equal(anyDuplicated(paste(proc_pk, s$procedures$procedure_number)), 0L)

  # schema + classes preserved per table
  for (t in names(tabs)) {
    expect_named(s[[t]], names(tabs[[t]]))
    expect_identical(vapply(s[[t]], class, character(1)),
                     vapply(tabs[[t]], class, character(1)))
  }
  expect_true(all(s$procedures$kind %in% c("PCI", "angio", "echo")))

  # check_linkage agrees
  rep <- check_linkage(res, verbose = FALSE)
  expect_true(attr(rep, "ok"))
  expect_true(all(rep$orphan_rows == 0L | is.na(rep$orphan_rows)))
  expect_false(any(rep$duplicate_keys))
})

test_that("the count model reproduces zero-inflated children-per-parent", {
  skip_if_not_installed("rpart")
  set.seed(9)
  np <- 200
  patients <- data.frame(id = seq_len(np), stringsAsFactors = FALSE)
  adm <- do.call(rbind, lapply(patients$id, function(pid) {
    na <- rpois(1, 0.8)                       # can be zero admissions
    if (na == 0L) return(NULL)
    data.frame(id = pid, admission_id = seq_len(na),
               x = round(rnorm(na, 10, 2), 1), stringsAsFactors = FALSE)
  }))
  real_frac0 <- mean(!(patients$id %in% adm$id))
  expect_gt(real_frac0, 0)                    # the signal is genuinely there

  res <- synth_linked(
    list(patients = patients, admissions = adm),
    list(patients = ~ id, admissions = ~ id / admission_id),
    list(patients = "id", admissions = c("id", "admission_id")),
    seed = 3
  )
  s <- res$syn
  syn_frac0 <- mean(!(s$patients$id %in% s$admissions$id))
  expect_gt(syn_frac0, 0)                     # some synthetic patients have none
  expect_lt(abs(syn_frac0 - real_frac0), 0.12)
})

test_that("cross-table predictors carry into the child variables", {
  skip_if_not_installed("rpart")
  set.seed(5)
  np <- 250
  patients <- data.frame(id = seq_len(np),
                         grp = sample(c("lo", "hi"), np, replace = TRUE),
                         stringsAsFactors = FALSE)
  # one admission per patient; its score is strongly driven by the parent group
  adm <- data.frame(
    id = patients$id, admission_id = 1L,
    score = round(ifelse(patients$grp == "hi",
                         rnorm(np, 100, 5), rnorm(np, 50, 5))),
    stringsAsFactors = FALSE
  )
  res <- synth_linked(
    list(patients = patients, admissions = adm),
    list(patients = ~ id, admissions = ~ id / admission_id),
    list(patients = "id", admissions = c("id", "admission_id")),
    method = "cart", seed = 7
  )
  s <- res$syn
  g <- s$patients$grp[match(s$admissions$id, s$patients$id)]
  gap <- mean(s$admissions$score[g == "hi"]) - mean(s$admissions$score[g == "lo"])
  expect_gt(gap, 30)                          # ~50 in the real data; survives
})

test_that("check_linkage flags orphans and duplicate keys in raw tables", {
  tabs <- list(
    patients   = data.frame(id = 1:3),
    admissions = data.frame(id = c(1, 2, 9), admission_id = c(1, 1, 1))
  )
  rep <- check_linkage(
    tabs, keys = list(patients = "id", admissions = c("id", "admission_id")),
    verbose = FALSE
  )
  expect_false(attr(rep, "ok"))
  expect_equal(rep$orphan_rows[rep$table == "admissions"], 1L)

  dupe <- list(a = data.frame(k = c(1, 1, 2)))
  rep2 <- check_linkage(dupe, keys = list(a = "k"), verbose = FALSE)
  expect_true(rep2$duplicate_keys[rep2$table == "a"])
})

test_that("hierarchy validation rejects duplicate keys and bad shapes", {
  bad <- list(a = data.frame(k = c(1, 1), v = 1:2))
  expect_error(
    synth_linked(bad, list(a = ~ k), list(a = "k")),
    "duplicate rows for its key"
  )
})
