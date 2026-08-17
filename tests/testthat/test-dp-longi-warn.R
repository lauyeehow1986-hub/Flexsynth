# D10: a DP linked-longitudinal child whose transition histogram will be swamped
# by noise (deep nesting / large branching cap / fine bins / small epsilon) warns
# the operator that within-unit autocorrelation may not survive; a well-budgeted
# release does not.

make_linked_longi <- function(np = 150) {
  set.seed(7)
  patients <- data.frame(id = seq_len(np),
                         sex = factor(sample(c("F", "M"), np, TRUE)))
  adm <- do.call(rbind, lapply(patients$id, function(pid) {
    k <- 1 + rpois(1, 0.7)
    data.frame(id = pid, admission_id = seq_len(k))
  }))
  labs <- do.call(rbind, lapply(seq_len(nrow(adm)), function(i) {
    m <- 3 + rpois(1, 2)
    v <- numeric(m); v[1] <- rnorm(1, 140, 12)
    for (t in 2:m) v[t] <- 0.8 * v[t - 1] + 0.2 * 140 + rnorm(1, 0, 6)
    data.frame(id = adm$id[i], admission_id = adm$admission_id[i],
               time = seq_len(m), sbp = round(v))
  }))
  list(
    tables = list(patients = patients, admissions = adm, labs = labs),
    structures = list(patients   = ~ id,
                      admissions = ~ id / admission_id,
                      labs       = ~ id / admission_id / time),
    keys = list(patients   = "id",
                admissions = c("id", "admission_id"),
                labs        = c("id", "admission_id", "time")))
}

test_that("a noise-swamped longitudinal child warns the operator", {
  d <- make_linked_longi()
  dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                   max_rows_per_person = list(admissions = 6, labs = 12),
                   longitudinal = TRUE,
                   bounds = list(sbp = c(80, 220)), bins = 12L)
  expect_warning(
    synth_linked(tables = d$tables, structures = d$structures, keys = d$keys,
                 privacy = dp, seed = 1),
    "autocorrelation may not survive")
})

test_that("a well-budgeted longitudinal child does not raise the noise warning", {
  d <- make_linked_longi(np = 400)
  dp <- dp_control(epsilon = 40, delta = 1e-6, mechanism = "gaussian",
                   max_rows_per_person = list(admissions = 2, labs = 4),
                   longitudinal = "labs",
                   bounds = list(sbp = c(80, 220)), bins = 5L)
  saw <- FALSE
  withCallingHandlers(
    synth_linked(tables = d$tables, structures = d$structures, keys = d$keys,
                 privacy = dp, seed = 1),
    warning = function(w) {
      if (grepl("autocorrelation may not survive", conditionMessage(w)))
        saw <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_false(saw)
})
