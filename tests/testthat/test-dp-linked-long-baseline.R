# Track B: baseline (held) columns and higher-order / cross-variable transitions
# on a longitudinally-modelled linked child - the flat DP Markov engine's two
# within-unit transition controls (dp_control(baseline / transition_order /
# transition_cross)) carried onto a linked child, applied per table.
#
# patients (root, `sex`) -> visits (longitudinal child) with three child columns:
#   * status    - 2-state sticky Markov chain whose FIRST value depends on sex,
#   * egfr_band - subject-invariant, one band per patient (a baseline covariate),
#   * grade     - a second time-varying factor.
# So visits has 3 variables, of which `egfr_band` is baseline (2 time-varying).
lk_bl <- function(np = 300, seed = 1, stick = 0.85, sex_effect = 0.85) {
  set.seed(seed)
  patients <- data.frame(id = seq_len(np),
                         sex = factor(sample(c("F", "M"), np, TRUE)),
                         stringsAsFactors = FALSE)
  bands <- c("G1", "G2", "G3")
  visits <- do.call(rbind, lapply(seq_len(np), function(i) {
    pid <- patients$id[i]; male <- patients$sex[i] == "M"
    k  <- 2L + stats::rpois(1, 1.3)
    p0 <- if (male) sex_effect else 1 - sex_effect      # P(start "worse")
    st <- character(k)
    st[1L] <- if (stats::runif(1) < p0) "worse" else "stable"
    for (j in 2:k)
      st[j] <- if (stats::runif(1) < stick) st[j - 1L]
               else setdiff(c("stable", "worse"), st[j - 1L])
    band <- sample(bands, 1L)                           # one band for the patient
    data.frame(id = pid, visit_num = seq_len(k),
               status    = factor(st, levels = c("stable", "worse")),
               egfr_band = factor(rep(band, k), levels = bands),
               grade     = factor(sample(c("A", "B", "C"), k, TRUE),
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

# rho for the (eps, delta) Gaussian releases below.
rho_for <- function(eps, delta)
  (sqrt(log(1 / delta) + eps) - sqrt(log(1 / delta)))^2

test_that("Laplace composition is exact for a longitudinal child with a baseline column", {
  d <- lk_bl(np = 150)
  cap <- 4L; eps <- 10
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = cap),
                   longitudinal = "visits", baseline = "egfr_band",
                   domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  # patients (root): sex, tree needs >1 var -> 1 one-way at path cap 1        = 1
  # visits (longi, pcp = 1, nC = 3, baseline egfr_band -> nT = 2):
  #   init = dp_longi_n_init(3, tree) = 3 + 3 = 6, at pcp 1                    = 6
  #   2 transitions at pcp*(cap - 1) = 3 each                                 = 6
  #   count histogram at cs = 1                                               = 1
  total_l1 <- 1 + 6 + 6 + 1
  expect_equal(res$privacy$noise, total_l1 / eps)
  # histograms: 1 (patients) + (6 init + 2 tran) + 1 count = 10
  expect_equal(res$privacy$n_marginals, 10L)
  bt <- by_tbl(res)
  expect_true(bt$visits$longitudinal)
  expect_equal(bt$visits$baseline, "egfr_band")
})

test_that("declaring a baseline column drops exactly one transition histogram", {
  d <- lk_bl(np = 150)
  cap <- 4L; eps <- 10
  mk <- function(base) dp_control(
    epsilon = eps, mechanism = "laplace", dependence = "tree",
    max_rows_per_person = c(visits = cap), longitudinal = "visits",
    baseline = base, domain = "public")
  rFull <- synth_linked(d$tables, d$structures, d$keys, privacy = mk(NULL), seed = 1)
  rBase <- synth_linked(d$tables, d$structures, d$keys,
                        privacy = mk("egfr_band"), seed = 1)
  # one fewer transition histogram (nT 3 -> 2)
  expect_equal(rFull$privacy$n_marginals - rBase$privacy$n_marginals, 1L)
  # total L1 drops by exactly one transition sensitivity ts = pcp*(cap-1) = 3
  expect_equal((rFull$privacy$noise - rBase$privacy$noise) * eps, 3)
  expect_equal(by_tbl(rFull)$visits$baseline, character(0))
  expect_equal(by_tbl(rBase)$visits$baseline, "egfr_band")
})

test_that("a higher transition_order lowers the transition sensitivity (Gaussian sum_sq)", {
  d <- lk_bl(np = 150)
  cap <- 4L; eps <- 8; delta <- 1e-5
  mk <- function(ord) dp_control(
    epsilon = eps, delta = delta, mechanism = "gaussian", dependence = "tree",
    max_rows_per_person = c(visits = cap), longitudinal = "visits",
    transition_order = ord, domain = "public")
  r1 <- synth_linked(d$tables, d$structures, d$keys, privacy = mk(1L), seed = 1)
  r2 <- synth_linked(d$tables, d$structures, d$keys, privacy = mk(2L), seed = 1)
  rho <- rho_for(eps, delta)
  # order 1: init 6, 3 transitions at ts = 3 -> sq = 1 + (6 + 3*9) + 1 = 35
  expect_equal(r1$privacy$noise, sqrt((1 + (6 + 3 * 9) + 1) / (2 * rho)))
  # order 2: transitions at ts = cap - 2 = 2 -> sq = 1 + (6 + 3*4) + 1 = 20
  expect_equal(r2$privacy$noise, sqrt((1 + (6 + 3 * 4) + 1) / (2 * rho)))
  # same number of histograms (order does not add/remove histograms)
  expect_equal(r1$privacy$n_marginals, r2$privacy$n_marginals)
  expect_equal(by_tbl(r2)$visits$tran_order, 2L)
})

test_that("transition_cross is budget-neutral on a linked child", {
  d <- lk_bl(np = 150)
  cap <- 4L; eps <- 10
  mk <- function(crs) dp_control(
    epsilon = eps, mechanism = "laplace", dependence = "tree",
    max_rows_per_person = c(visits = cap), longitudinal = "visits",
    transition_cross = crs, domain = "public")
  r0 <- synth_linked(d$tables, d$structures, d$keys, privacy = mk(0L), seed = 1)
  r1 <- synth_linked(d$tables, d$structures, d$keys, privacy = mk(1L), seed = 1)
  expect_equal(r1$privacy$noise, r0$privacy$noise)          # same budget
  expect_equal(r1$privacy$n_marginals, r0$privacy$n_marginals)
  expect_equal(by_tbl(r1)$visits$tran_cross, 1L)
})

test_that("transition_order beyond the branching cap is refused (per table)", {
  d <- lk_bl(np = 60)
  dp <- dp_control(epsilon = 8, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 3L), longitudinal = "visits",
                   transition_order = 3L, domain = "public")   # needs cap >= 4
  expect_error(
    synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1),
    "transition_order")
})

test_that("a baseline column is byte-for-byte constant within each synthetic unit", {
  d <- lk_bl(np = 200, seed = 3)
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   dependence = "tree", max_rows_per_person = c(visits = 6L),
                   longitudinal = "visits", baseline = "egfr_band",
                   domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 7)
  syn <- as.list(res)
  nuniq <- tapply(as.character(syn$visits$egfr_band), syn$visits$id,
                  function(x) length(unique(x)))
  expect_true(all(nuniq == 1L))                    # held constant within a unit
  # a time-varying column is NOT forced constant (sanity: it does vary somewhere)
  vuniq <- tapply(as.character(syn$visits$status), syn$visits$id,
                  function(x) length(unique(x)))
  expect_gt(max(vuniq), 1L)
})

test_that("an order-2 linked child preserves within-unit autocorrelation and reproduces", {
  d <- lk_bl(np = 320, seed = 4)
  expect_gt(agree_status(d$tables$visits), 0.75)   # real status is sticky
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   dependence = "tree", max_rows_per_person = c(visits = 6L),
                   longitudinal = "visits", transition_order = 2L,
                   domain = "public")
  a <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 11)
  b <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 11)
  expect_equal(as.list(a)$visits, as.list(b)$visits)          # deterministic
  expect_gt(agree_status(as.list(a)$visits), 0.6)             # stickiness survives
  syn <- as.list(a)
  pos <- split(syn$visits$visit_num, syn$visits$id)
  expect_true(all(vapply(pos, function(v) identical(sort(v), seq_along(v)),
                         logical(1))))
})

test_that("baseline composes with a cross-conditioned (combined) longitudinal child", {
  d <- lk_bl(np = 150)
  cap <- 4L; eps <- 10
  dp <- dp_control(epsilon = eps, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = cap), longitudinal = "visits",
                   cross_table = TRUE, baseline = "egfr_band", domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)

  # patients: sex 1 one-way at cap 1                                          = 1
  # visits (combined cross + longi, nC = 3, nP = 1, baseline -> nT = 2):
  #   init (cross) = dp_child_nvarmarg(3, 1, tree) = (3 + 3) + 3*1 = 9, pcp 1  = 9
  #   2 transitions at pcp*(cap-1) = 3 each                                    = 6
  #   count at cs = 1                                                          = 1
  total_l1 <- 1 + 9 + 6 + 1
  expect_equal(res$privacy$noise, total_l1 / eps)
  bt <- by_tbl(res)
  expect_true(bt$visits$cross_init)               # initial state cond. on parent
  expect_equal(bt$visits$baseline, "egfr_band")   # and still holds the baseline
  # referential integrity holds
  cl <- check_linkage(res)
  expect_true(all(cl$orphan_rows[!is.na(cl$orphan_rows)] == 0))
})

test_that("the accounting print flags a linked child's baseline and deepened transitions", {
  d <- lk_bl(np = 80, seed = 5)
  dp <- dp_control(epsilon = 8, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 5L), longitudinal = "visits",
                   baseline = "egfr_band", transition_order = 2L,
                   transition_cross = 1L, domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  out <- paste(utils::capture.output(print(res$privacy)), collapse = "\n")
  expect_match(out, "baseline held: egfr_band")
  expect_match(out, "transitions: order 2 \\+ 1 cross-parent")
})
