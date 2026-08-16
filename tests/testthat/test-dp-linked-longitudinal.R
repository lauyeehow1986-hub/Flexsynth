# Track B linked + longitudinal DP: a within-unit DP Markov model for a linked
# child table (dp_control(longitudinal = ...)).

# patients (root) -> visits (longitudinal child). `status` is a sticky 2-state
# Markov chain within a patient; `grade` a 3-level factor. Every patient has >= 2
# visits so transitions are always identified.
lk_longi <- function(np = 120, seed = 1, stick = 0.85) {
  set.seed(seed)
  patients <- data.frame(id = seq_len(np),
                         sex = factor(sample(c("F", "M"), np, TRUE)),
                         stringsAsFactors = FALSE)
  visits <- do.call(rbind, lapply(patients$id, function(pid) {
    k  <- 2L + stats::rpois(1, 1.3)
    st <- character(k); st[1L] <- sample(c("stable", "worse"), 1)
    for (i in 2:k)
      st[i] <- if (stats::runif(1) < stick) st[i - 1L]
               else setdiff(c("stable", "worse"), st[i - 1L])
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

# lag-1 within-patient agreement rate of `status` (autocorrelation proxy).
agree_status <- function(v) {
  s <- split(as.character(v$status), v$id)
  num <- 0L; den <- 0L
  for (x in s) if (length(x) >= 2L) {
    num <- num + sum(x[-1L] == x[-length(x)]); den <- den + length(x) - 1L
  }
  if (den == 0L) NA_real_ else num / den
}

by_tbl <- function(res) {
  info <- res$privacy$linked$tables
  stats::setNames(info, vapply(info, function(x) x$name, character(1)))
}

test_that("Laplace composition is exact for a longitudinal child (init + transitions)", {
  d <- lk_longi()
  cap <- 4L; eps <- 10
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = cap),
                   longitudinal = "visits", domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  # patients (root): sex, tree needs >1 var -> 1 one-way at path cap 1  = 1
  # visits (longi, pcp = 1): 2 vars, tree -> n_init = 2 + 1 = 3 at pcp 1;
  #   2 transitions at pcp*(cap-1) = 3 each; + count histogram at cs = 1.
  #   var l1 = 3*1 + 2*3 = 9 ; count l1 = 1
  total_l1 <- 1 + 9 + 1
  expect_equal(res$privacy$noise, total_l1 / eps)
  # histograms: 1 (patients) + (3 init + 2 tran) + 1 count = 7
  expect_equal(res$privacy$n_marginals, 7L)
  expect_true(by_tbl(res)$visits$longitudinal)
  expect_false(by_tbl(res)$patients$longitudinal)
})

test_that("Gaussian composition is exact (summed squared L2) for a longitudinal child", {
  d <- lk_longi()
  cap <- 4L; eps <- 8; delta <- 1e-5
  dp <- dp_control(epsilon = eps, delta = delta, mechanism = "gaussian",
                   dependence = "independent",
                   max_rows_per_person = c(visits = cap),
                   longitudinal = "visits", domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  # patients: sex 1 one-way at cap 1 -> sq 1
  # visits (indep): n_init = 2 at pcp 1 -> 2*1 ; 2 transitions at 3 -> 2*9
  #   count at cs 1 -> 1
  sum_sq <- 1 + (2 * 1 + 2 * 9) + 1
  rho <- (sqrt(log(1 / delta) + eps) - sqrt(log(1 / delta)))^2
  expect_equal(res$privacy$noise, sqrt(sum_sq / (2 * rho)))
})

test_that("longitudinal modelling adds transition histograms over the exchangeable release", {
  d <- lk_longi()
  base <- dp_control(epsilon = 10, mechanism = "laplace", dependence = "tree",
                     max_rows_per_person = c(visits = 4L), domain = "public")
  long <- dp_control(epsilon = 10, mechanism = "laplace", dependence = "tree",
                     max_rows_per_person = c(visits = 4L),
                     longitudinal = "visits", domain = "public")
  r0 <- synth_linked(d$tables, d$structures, d$keys, privacy = base, seed = 1)
  rL <- synth_linked(d$tables, d$structures, d$keys, privacy = long, seed = 1)
  # exchangeable: patients 1 + visits(tree 2 vars = 3) + count 1 = 5
  expect_equal(r0$privacy$n_marginals, 5L)
  expect_equal(rL$privacy$n_marginals, 7L)
  expect_gt(rL$privacy$n_marginals, r0$privacy$n_marginals)
})

test_that("the DP Markov model recovers within-unit autocorrelation the exchangeable model destroys", {
  d <- lk_longi(np = 260, seed = 4)
  real_a <- agree_status(d$tables$visits)
  expect_gt(real_a, 0.75)                              # the data is sticky

  long <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                     max_rows_per_person = c(visits = 6L),
                     longitudinal = "visits", domain = "public")
  exch <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                     max_rows_per_person = c(visits = 6L), domain = "public")
  rL <- synth_linked(d$tables, d$structures, d$keys, privacy = long, seed = 21)
  r0 <- synth_linked(d$tables, d$structures, d$keys, privacy = exch, seed = 21)

  aL <- agree_status(as.list(rL)$visits)
  a0 <- agree_status(as.list(r0)$visits)
  expect_gt(aL, 0.68)                                  # autocorrelation preserved
  expect_lt(a0, 0.60)                                  # exchangeable ~ chance
  expect_gt(aL - a0, 0.15)
})

test_that("path caps scale the transition sensitivity when the longitudinal table is a grandchild", {
  set.seed(7)
  nh <- 6
  hospitals <- data.frame(hosp_id = seq_len(nh),
                          region = factor(sample(c("N", "S"), nh, TRUE)))
  patients <- do.call(rbind, lapply(hospitals$hosp_id, function(h) {
    n <- 2L + stats::rpois(1, 1)
    data.frame(hosp_id = h, patient_id = seq_len(n),
               sex = factor(sample(c("F", "M"), n, TRUE)))
  }))
  visits <- do.call(rbind, lapply(seq_len(nrow(patients)), function(i) {
    k  <- 2L + stats::rpois(1, 1)
    st <- character(k); st[1L] <- sample(c("stable", "worse"), 1)
    for (j in 2:k)
      st[j] <- if (stats::runif(1) < 0.85) st[j - 1L]
               else setdiff(c("stable", "worse"), st[j - 1L])
    data.frame(hosp_id = patients$hosp_id[i], patient_id = patients$patient_id[i],
               visit_num = seq_len(k),
               status = factor(st, levels = c("stable", "worse")))
  }))
  tables <- list(hospitals = hospitals, patients = patients, visits = visits)
  structures <- list(hospitals = ~ hosp_id, patients = ~ hosp_id / patient_id,
                     visits = ~ hosp_id / patient_id / visit_num)
  keys <- list(hospitals = "hosp_id", patients = c("hosp_id", "patient_id"),
               visits = c("hosp_id", "patient_id", "visit_num"))
  eps <- 12
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "independent",
                   max_rows_per_person = list(patients = 3, visits = 4),
                   longitudinal = "visits", domain = "public")
  res <- synth_linked(tables, structures, keys, privacy = dp, seed = 1)

  bt <- by_tbl(res)
  expect_equal(bt$patients$path_cap, 3L)
  expect_equal(bt$visits$path_cap, 12L)                # 3 * 4
  expect_equal(bt$visits$count_sensitivity, 3L)        # path_cap[patients]
  expect_true(bt$visits$longitudinal)

  # hospitals: region 1 var at cap 1                       = 1
  # patients : sex 1 var at path cap 3 + count (cs 1)      = 3 + 1
  # visits   : status init 1 at pcp 3 + 1 transition at
  #            pcp*(cap-1) = 3*3 = 9 + count (cs 3)         = 3 + 9 + 3
  total_l1 <- 1 + (3 + 1) + (3 + 9 + 3)
  expect_equal(res$privacy$noise, total_l1 / eps)
  expect_true(all(res$syn$visits$patient_id %in% res$syn$patients$patient_id))
})

test_that("longitudinal children keep referential integrity, contiguous positions, and reproduce", {
  d <- lk_longi(np = 80, seed = 2)
  dp <- dp_control(epsilon = 6, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 5L),
                   longitudinal = "visits", domain = "public")
  a <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 9)
  b <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 9)
  expect_equal(as.list(a)$visits, as.list(b)$visits)   # deterministic under a seed

  syn <- as.list(a)
  expect_true(all(syn$visits$id %in% syn$patients$id)) # no orphans
  pos <- split(syn$visits$visit_num, syn$visits$id)
  expect_true(all(vapply(pos, function(v) identical(sort(v), seq_along(v)),
                         logical(1))))                  # 1..count, no gaps
  cl <- check_linkage(a)
  expect_true(all(cl$orphan_rows[!is.na(cl$orphan_rows)] == 0))
})

test_that("over-cap longitudinal units are prefix-truncated in temporal order", {
  set.seed(3)
  patients <- data.frame(id = 1:5, sex = factor(rep("F", 5)))
  # one long trajectory (10 visits) with an increasing marker; cap will drop the tail
  visits <- data.frame(id = rep(1:5, each = 10),
                       visit_num = rep(1:10, times = 5),
                       status = factor(rep(c("stable", "worse"), length.out = 50),
                                       levels = c("stable", "worse")))
  d <- list(tables = list(patients = patients, visits = visits),
            structures = list(patients = ~ id, visits = ~ id / visit_num),
            keys = list(patients = "id", visits = c("id", "visit_num")))
  dp <- dp_control(epsilon = 8, mechanism = "laplace", dependence = "independent",
                   max_rows_per_person = c(visits = 4L),
                   longitudinal = "visits", domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  # 5 units * (10 - 4) = 30 rows dropped by prefix truncation
  expect_equal(by_tbl(res)$visits$rows_dropped, 30L)
})

test_that("longitudinal = TRUE + cross_table combine: initial state is cross-conditioned", {
  d <- lk_longi(np = 60, seed = 5)
  dp <- dp_control(epsilon = 8, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 4L),
                   longitudinal = TRUE, cross_table = TRUE, domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  bt <- by_tbl(res)
  expect_true(bt$visits$longitudinal)
  expect_true(bt$visits$cross_init)                    # init state conditions on parent
  expect_false(bt$visits$cross)                        # not plain (all-row) cross
  expect_false(bt$patients$longitudinal)               # root is never longitudinal
})

test_that("longitudinal names are validated by the engine", {
  d <- lk_longi(np = 40)
  # unknown table
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = dp_control(epsilon = 4, mechanism = "laplace",
                                      max_rows_per_person = c(visits = 4L),
                                      longitudinal = "nope", domain = "public"),
                 seed = 1),
    "unknown table")
  # a root cannot be longitudinal
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = dp_control(epsilon = 4, mechanism = "laplace",
                                      max_rows_per_person = c(visits = 4L),
                                      longitudinal = "patients", domain = "public"),
                 seed = 1),
    "is a root")
  # branching cap must be >= 2 to measure a transition
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = dp_control(epsilon = 4, mechanism = "laplace",
                                      max_rows_per_person = c(visits = 1L),
                                      longitudinal = "visits", domain = "public"),
                 seed = 1),
    "cap >= 2")
})

test_that("dp_control validates longitudinal and it is inert for a flat synth()", {
  expect_error(dp_control(epsilon = 1, longitudinal = 1),
               "single TRUE/FALSE")
  expect_error(dp_control(epsilon = 1, longitudinal = c("a", NA)),
               "single TRUE/FALSE")
  expect_s3_class(dp_control(epsilon = 1, longitudinal = TRUE), "dp_control")

  # A flat DP release ignores `longitudinal` (no parent table to model over).
  df <- data.frame(id = 1:60,
                   g = factor(sample(c("x", "y"), 60, TRUE)))
  res <- synth(df, structure = ~ id, method = "cart",
               privacy = dp_control(epsilon = 4, mechanism = "laplace",
                                    longitudinal = TRUE, domain = "public"),
               seed = 1)
  expect_s3_class(res, "synth_result")
  expect_s3_class(res$privacy, "dp_accounting")
})
