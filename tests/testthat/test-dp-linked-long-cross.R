# Track B: combined cross-table + longitudinal DP on the same linked child table.
# dp_control(longitudinal = ..., cross_table = TRUE) makes a child's rows a
# within-unit DP Markov trajectory whose INITIAL STATE is cross-conditioned on the
# synthetic parent; the transition chain then carries that parent dependence
# forward. Only the initial-state parent-by-child joints cost extra budget.

# patients (root, `sex`) -> visits (longitudinal child, `status`+`grade`).
# The FIRST visit's status depends strongly on the patient's sex, then status is a
# sticky 2-state Markov chain, so sex -> status dependence is set at baseline and
# propagates through the trajectory. Every patient has >= 2 visits.
lk_lc <- function(np = 300, seed = 1, stick = 0.85, sex_effect = 0.85) {
  set.seed(seed)
  patients <- data.frame(id = seq_len(np),
                         sex = factor(sample(c("F", "M"), np, TRUE)),
                         stringsAsFactors = FALSE)
  visits <- do.call(rbind, lapply(seq_len(np), function(i) {
    pid <- patients$id[i]; male <- patients$sex[i] == "M"
    k  <- 2L + stats::rpois(1, 1.3)
    p0 <- if (male) sex_effect else 1 - sex_effect      # P(start "worse")
    st <- character(k)
    st[1L] <- if (stats::runif(1) < p0) "worse" else "stable"
    for (j in 2:k)
      st[j] <- if (stats::runif(1) < stick) st[j - 1L]
               else setdiff(c("stable", "worse"), st[j - 1L])
    data.frame(id = pid, visit_num = seq_len(k),
               status = factor(st, levels = c("stable", "worse")),
               grade  = factor(sample(c("A", "B", "C"), k, TRUE),
                               levels = c("A", "B", "C")),
               stringsAsFactors = FALSE)
  }))
  list(tables = list(patients = patients, visits = visits),
       structures = list(patients = ~ id, visits = ~ id / visit_num),
       keys = list(patients = "id", visits = c("id", "visit_num")))
}

by_tbl <- function(res) {
  info <- res$privacy$linked$tables
  stats::setNames(info, vapply(info, function(x) x$name, character(1)))
}

# lag-1 within-patient agreement rate of `status` (autocorrelation proxy).
agree_status <- function(v) {
  s <- split(as.character(v$status), v$id)
  num <- 0L; den <- 0L
  for (x in s) if (length(x) >= 2L) {
    num <- num + sum(x[-1L] == x[-length(x)]); den <- den + length(x) - 1L
  }
  if (den == 0L) NA_real_ else num / den
}

# P(status == "worse" | sex == "M") - P(... | sex == "F") over all rows, joining
# sex onto visits by id. A parent -> child dependence proxy (0 when independent).
wdiff <- function(patients, visits) {
  sex <- patients$sex[match(visits$id, patients$id)]
  worse <- as.character(visits$status) == "worse"
  mean(worse[sex == "M"]) - mean(worse[sex == "F"])
}

test_that("Laplace composition is exact for a combined long+cross child", {
  d <- lk_lc(np = 120)
  cap <- 4L; eps <- 10
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = cap),
                   longitudinal = "visits", cross_table = TRUE, domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  # patients (root): sex, tree needs >1 var -> 1 one-way at path cap 1     = 1
  # visits (longi + cross, pcp = 1, nC = 2, nP = 1):
  #   init (cross) = dp_child_nvarmarg(2,1,tree) = (2 + 1) + 2*1 = 5, at pcp 1
  #   2 transitions at pcp*(cap-1) = 3 each                        = 6
  #   count histogram at cs = 1                                    = 1
  total_l1 <- 1 + 5 + 6 + 1
  expect_equal(res$privacy$noise, total_l1 / eps)
  # histograms: 1 (patients) + (5 init + 2 tran) + 1 count = 9
  expect_equal(res$privacy$n_marginals, 9L)
  bt <- by_tbl(res)
  expect_true(bt$visits$longitudinal)
  expect_true(bt$visits$cross_init)
})

test_that("Gaussian composition is exact (summed squared L2) for a combined child", {
  d <- lk_lc(np = 120)
  cap <- 4L; eps <- 8; delta <- 1e-5
  dp <- dp_control(epsilon = eps, delta = delta, mechanism = "gaussian",
                   dependence = "independent",
                   max_rows_per_person = c(visits = cap),
                   longitudinal = "visits", cross_table = TRUE, domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  # patients: sex 1 one-way at cap 1 -> sq 1
  # visits (indep): init (cross) = dp_child_nvarmarg(2,1,indep) = 2 + 2 = 4 at pcp 1
  #   2 transitions at 3 -> 2*9 ; count at cs 1 -> 1
  sum_sq <- 1 + (4 * 1 + 2 * 9) + 1
  rho <- (sqrt(log(1 / delta) + eps) - sqrt(log(1 / delta)))^2
  expect_equal(res$privacy$noise, sqrt(sum_sq / (2 * rho)))
})

test_that("combining adds exactly the parent-by-child initial joints over longi-only", {
  d <- lk_lc(np = 120)
  eps <- 10; cap <- 4L
  longi <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "tree",
                      max_rows_per_person = c(visits = cap),
                      longitudinal = "visits", domain = "public")
  both  <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "tree",
                      max_rows_per_person = c(visits = cap),
                      longitudinal = "visits", cross_table = TRUE, domain = "public")
  rL <- synth_linked(d$tables, d$structures, d$keys, privacy = longi, seed = 1)
  rB <- synth_linked(d$tables, d$structures, d$keys, privacy = both,  seed = 1)
  # nC * nP = 2 * 1 = 2 extra initial-state joints, nothing else changes.
  expect_equal(rB$privacy$n_marginals - rL$privacy$n_marginals, 2L)
  expect_false(by_tbl(rL)$visits$cross_init)
  expect_true(by_tbl(rB)$visits$cross_init)
})

test_that("the combined model keeps BOTH parent dependence and within-unit autocorrelation", {
  d <- lk_lc(np = 340, seed = 4)
  real_diff <- wdiff(d$tables$patients, d$tables$visits)
  real_agree <- agree_status(d$tables$visits)
  expect_gt(real_diff, 0.4)                             # sex drives status
  expect_gt(real_agree, 0.75)                           # status is sticky

  both  <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                      max_rows_per_person = c(visits = 6L),
                      longitudinal = "visits", cross_table = TRUE,
                      domain = "public")
  longi <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                      max_rows_per_person = c(visits = 6L),
                      longitudinal = "visits", domain = "public")
  rB <- synth_linked(d$tables, d$structures, d$keys, privacy = both,  seed = 21)
  rL <- synth_linked(d$tables, d$structures, d$keys, privacy = longi, seed = 21)

  synB <- as.list(rB); synL <- as.list(rL)
  # parent -> child dependence: preserved with cross, destroyed without it.
  dB <- wdiff(synB$patients, synB$visits)
  dL <- wdiff(synL$patients, synL$visits)
  expect_gt(dB, 0.25)                                   # sex->status survives
  expect_lt(abs(dL), 0.15)                              # longi-only ~ independent
  expect_gt(dB - dL, 0.2)
  # autocorrelation preserved in both (the transition chain is unchanged).
  expect_gt(agree_status(synB$visits), 0.65)
})

test_that("combined children keep referential integrity, contiguous positions, and reproduce", {
  d <- lk_lc(np = 90, seed = 2)
  dp <- dp_control(epsilon = 6, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 5L),
                   longitudinal = "visits", cross_table = TRUE, domain = "public")
  a <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 9)
  b <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 9)
  expect_equal(as.list(a)$visits, as.list(b)$visits)    # deterministic under a seed

  syn <- as.list(a)
  expect_true(all(syn$visits$id %in% syn$patients$id))  # no orphans
  pos <- split(syn$visits$visit_num, syn$visits$id)
  expect_true(all(vapply(pos, function(v) identical(sort(v), seq_along(v)),
                         logical(1))))                   # 1..count, no gaps
  cl <- check_linkage(a)
  expect_true(all(cl$orphan_rows[!is.na(cl$orphan_rows)] == 0))
})

test_that("the accounting print flags a combined child as cross-conditioned initial state", {
  d <- lk_lc(np = 60, seed = 5)
  dp <- dp_control(epsilon = 8, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 4L),
                   longitudinal = "visits", cross_table = TRUE, domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  out <- paste(utils::capture.output(print(res$privacy)), collapse = "\n")
  expect_match(out, "DP Markov over rows \\(initial state cond\\. on parent\\)")
})
