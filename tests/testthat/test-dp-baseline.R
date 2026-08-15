# Track B: baseline (subject-invariant) columns held exactly constant within a
# synthetic unit under the flat DP Markov longitudinal engine.

# A longitudinal table with two genuinely subject-invariant baseline columns
# (`sex`, `age`, fixed per person) and one time-varying column (`sbp`, a random
# walk across visits). 2..4 visits per person.
mk_longi <- function(n = 60, seed = 1) {
  set.seed(seed)
  parts <- lapply(seq_len(n), function(i) {
    v    <- sample(2:4, 1L)
    sex  <- sample(c("F", "M"), 1L)
    age0 <- round(stats::rnorm(1L, 60, 8))
    sbp  <- round(130 + cumsum(stats::rnorm(v, 0, 8)))
    data.frame(id = i, visit = seq_len(v),
               sex = factor(sex, levels = c("F", "M")),
               age = age0, sbp = sbp)
  })
  do.call(rbind, parts)
}

longi_bounds <- list(age = c(18, 100), sbp = c(60, 240))

dp_longi <- function(baseline = NULL, mechanism = "laplace", dependence = "tree",
                     epsilon = 3, delta = if (mechanism == "gaussian") 1e-6 else 0) {
  dp_control(epsilon = epsilon, delta = delta, mechanism = mechanism,
             dependence = dependence, baseline = baseline,
             max_rows_per_person = 4L, bounds = longi_bounds, domain = "public")
}

test_that("declared baseline columns are exactly constant within every unit", {
  df <- mk_longi()
  res <- synth(df, structure = ~ id / visit,
               privacy = dp_longi(baseline = c("sex", "age")), seed = 11)
  out <- as.data.frame(res)
  by_unit <- split(out, out$id)
  const_sex <- vapply(by_unit, function(u) length(unique(u$sex)) == 1L, logical(1))
  const_age <- vapply(by_unit, function(u) length(unique(u$age)) == 1L, logical(1))
  expect_true(all(const_sex))
  expect_true(all(const_age))
})

test_that("without baseline the same columns are allowed to drift within a unit", {
  df <- mk_longi()
  res <- synth(df, structure = ~ id / visit,
               privacy = dp_longi(baseline = NULL), seed = 11)
  out <- as.data.frame(res)
  # At least one multi-row unit must show within-unit variation in age, proving
  # the transition model (not a held constant) is driving the column.
  by_unit <- split(out, out$id)
  multi <- by_unit[vapply(by_unit, nrow, integer(1)) > 1L]
  drift <- vapply(multi, function(u) length(unique(u$age)) > 1L, logical(1))
  expect_true(any(drift))
})

test_that("baseline drops those columns' transition histograms and reduces n_marginals", {
  df <- mk_longi()
  full <- synth(df, structure = ~ id / visit,
                privacy = dp_longi(baseline = NULL), seed = 1)$privacy
  held <- synth(df, structure = ~ id / visit,
                privacy = dp_longi(baseline = c("sex", "age")), seed = 1)$privacy
  # 3 variables (sex, age, sbp): full keeps 3 transitions, held keeps 1 (sbp).
  expect_equal(full$longitudinal$n_transitions, 3L)
  expect_equal(held$longitudinal$n_transitions, 1L)
  expect_identical(held$longitudinal$baseline, c("sex", "age"))
  expect_lt(held$n_marginals, full$n_marginals)
})

test_that("Laplace composition is exact and baseline lowers the noise scale", {
  df <- mk_longi()
  # dependence = independent so n_init_marg = nV exactly (no pairwise term).
  full <- synth(df, structure = ~ id / visit,
                privacy = dp_longi(baseline = NULL, dependence = "independent",
                                   epsilon = 3), seed = 1)$privacy
  held <- synth(df, structure = ~ id / visit,
                privacy = dp_longi(baseline = c("sex", "age"),
                                   dependence = "independent", epsilon = 3),
                seed = 1)$privacy
  cap <- 4; nV <- 3
  # total_l1 = 1 length + nV initial + nT*(cap-1) transitions; scale = total_l1/eps
  full_l1 <- 1 + nV + 3 * (cap - 1)          # 1 + 3 + 9 = 13
  held_l1 <- 1 + nV + 1 * (cap - 1)          # 1 + 3 + 3 = 7
  expect_equal(full$noise, full_l1 / 3, tolerance = 1e-9)
  expect_equal(held$noise, held_l1 / 3, tolerance = 1e-9)
  expect_lt(held$noise, full$noise)
})

test_that("Gaussian composition is exact under baseline (sum of squared L2)", {
  df <- mk_longi()
  held <- synth(df, structure = ~ id / visit,
                privacy = dp_longi(baseline = c("sex", "age"),
                                   dependence = "independent",
                                   mechanism = "gaussian", epsilon = 3,
                                   delta = 1e-6), seed = 1)$privacy
  cap <- 4; nV <- 3; nT <- 1
  sum_sq <- 1 + nV + nT * (cap - 1)^2        # 1 + 3 + 9 = 13
  rho <- zcdp_rho_for(3, 1e-6)
  expect_equal(held$noise, sqrt(sum_sq / (2 * rho)), tolerance = 1e-9)
})

test_that("all-baseline is a valid degenerate release (no transitions)", {
  df <- mk_longi()
  res <- synth(df, structure = ~ id / visit,
               privacy = dp_longi(baseline = c("sex", "age", "sbp")), seed = 5)
  acct <- res$privacy
  expect_equal(acct$longitudinal$n_transitions, 0L)
  out <- as.data.frame(res)
  by_unit <- split(out, out$id)
  all_const <- vapply(by_unit, function(u)
    all(vapply(c("sex", "age", "sbp"),
               function(v) length(unique(u[[v]])) == 1L, logical(1))),
    logical(1))
  expect_true(all(all_const))
})

test_that("a baseline name that matches no modelled column is a harmless no-op", {
  df <- mk_longi()
  a <- synth(df, structure = ~ id / visit,
             privacy = dp_longi(baseline = NULL), seed = 7)$privacy
  b <- synth(df, structure = ~ id / visit,
             privacy = dp_longi(baseline = c("nonexistent", "id", "visit")),
             seed = 7)$privacy
  # id / visit are not modelled vars; "nonexistent" matches nothing -> identical
  # release to no baseline at all.
  expect_equal(b$longitudinal$n_transitions, a$longitudinal$n_transitions)
  expect_equal(b$n_marginals, a$n_marginals)
})

test_that("baseline synthesis is reproducible and preserves shape", {
  df <- mk_longi()
  r1 <- as.data.frame(synth(df, structure = ~ id / visit,
                            privacy = dp_longi(baseline = c("sex", "age")),
                            seed = 3))
  r2 <- as.data.frame(synth(df, structure = ~ id / visit,
                            privacy = dp_longi(baseline = c("sex", "age")),
                            seed = 3))
  expect_identical(r1, r2)
  expect_identical(names(r1), names(df))
  expect_true(all(r1$age >= 18 & r1$age <= 100))
})

test_that("dp_control validates the baseline argument", {
  expect_silent(dp_control(epsilon = 1, baseline = NULL))
  expect_silent(dp_control(epsilon = 1, baseline = c("a", "b")))
  expect_error(dp_control(epsilon = 1, baseline = 1:3), "baseline")
  expect_error(dp_control(epsilon = 1, baseline = c("a", NA)), "baseline")
  expect_error(dp_control(epsilon = 1, baseline = c("a", "")), "baseline")
})

test_that("baseline is announced in the control and accounting prints", {
  ctrl <- capture.output(print(dp_control(epsilon = 1, baseline = c("sex", "age"))))
  expect_true(any(grepl("baseline", ctrl, ignore.case = TRUE)))

  df <- mk_longi()
  acct <- synth(df, structure = ~ id / visit,
                privacy = dp_longi(baseline = c("sex", "age")), seed = 1)$privacy
  txt <- capture.output(print(acct))
  expect_true(any(grepl("held constant within unit", txt)))
  expect_true(any(grepl("sex", txt)))
})
