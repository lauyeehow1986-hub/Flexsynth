# Phase 7: Track B differentially private synthesis. These tests check the
# accounting math, the discretisation / contribution-bounding layer, and
# end-to-end DP synthesis via synth(privacy = dp_control(...)).

# ---- accounting: zCDP <-> (eps, delta) ------------------------------------

test_that("zcdp_rho_for inverts the zCDP-to-(eps,delta) bound", {
  for (eps in c(0.5, 1, 3)) for (delta in c(1e-4, 1e-6)) {
    rho <- flexsynth:::zcdp_rho_for(eps, delta)
    # The rho we spend must map back to exactly eps under the standard bound.
    eps_back <- rho + 2 * sqrt(rho * log(1 / delta))
    expect_equal(eps_back, eps, tolerance = 1e-8)
    expect_gt(rho, 0)
  }
})

test_that("dp_calibrate sets Laplace scale and Gaussian sd from the budget", {
  dpl <- dp_control(epsilon = 2, mechanism = "laplace")
  cl <- flexsynth:::dp_calibrate(dpl, n_marginals = 5L, cap = 3)
  expect_equal(cl$mechanism, "laplace")
  expect_equal(cl$scale, (5 * 3) / 2)                     # n_marginals*cap/eps

  dpg <- dp_control(epsilon = 2, delta = 1e-6, mechanism = "gaussian")
  cg <- flexsynth:::dp_calibrate(dpg, n_marginals = 5L, cap = 3)
  rho <- flexsynth:::zcdp_rho_for(2, 1e-6)
  expect_equal(cg$sigma, 3 * sqrt(5 / (2 * rho)))         # cap*sqrt(m/(2 rho))
  expect_length(cg$add_noise(rep(0, 4)), 4L)
})

test_that("rlaplace has the right spread", {
  set.seed(1)
  x <- flexsynth:::rlaplace(2e5, scale = 2)
  expect_lt(abs(mean(x)), 0.1)
  expect_equal(var(x), 2 * 2^2, tolerance = 0.15)         # Var = 2 b^2
})

# ---- contribution bounding -------------------------------------------------

test_that("dp_contribution_bound caps rows per person", {
  df <- data.frame(id = c(1, 1, 1, 2, 2, 3), x = 1:6)
  res <- flexsynth:::dp_contribution_bound(df, "id", cap = 2L)
  tab <- table(res$data$id)
  expect_true(all(tab <= 2L))
  expect_equal(res$dropped, 1L)                           # person 1 had 3 rows

  none <- flexsynth:::dp_contribution_bound(df, "id", cap = 5L)
  expect_equal(none$dropped, 0L)
  expect_equal(nrow(none$data), 6L)
})

# ---- discretisation --------------------------------------------------------

test_that("numeric encode/decode stays within public bounds", {
  dom <- flexsynth:::dp_domain_column(c(10, 55, 90), nbin = 20L, bound = c(0, 100))
  expect_equal(dom$kind, "numeric")
  code <- flexsynth:::dp_encode(dom, c(-5, 50, 200))       # clamped to [0,100]
  expect_true(all(code >= 1L & code <= 20L))
  val <- flexsynth:::dp_decode(dom, code)
  expect_true(all(val >= 0 & val <= 100))
})

test_that("integer columns decode back to integers", {
  dom <- flexsynth:::dp_domain_column(1:10L, nbin = 5L, bound = c(1, 10))
  expect_true(dom$is_integer)
  val <- flexsynth:::dp_decode(dom, c(1L, 3L, 5L))
  expect_true(is.integer(val))
  expect_true(all(val >= 1L & val <= 10L))
})

test_that("factor columns round-trip through their levels", {
  f <- factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  dom <- flexsynth:::dp_domain_column(f, nbin = 20L)
  expect_equal(dom$kind, "factor")
  code <- flexsynth:::dp_encode(dom, f)
  expect_equal(code, 1:3)
  out <- flexsynth:::dp_decode(dom, code)
  expect_s3_class(out, "factor")
  expect_identical(levels(out), levels(f))
})

# ---- end to end ------------------------------------------------------------

flat_data <- function(n = 400) {
  data.frame(
    id  = seq_len(n),
    age = round(rnorm(n, 60, 8)),
    sbp = rnorm(n, 130, 15),
    sex = factor(sample(c("F", "M"), n, TRUE), levels = c("F", "M")),
    stringsAsFactors = FALSE
  )
}

test_that("synth() with dp_control returns a Track B result", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, bounds = list(age = c(20, 100), sbp = c(60, 220)))
  res <- synth(df, structure = ~ id, privacy = dp, seed = 1)
  expect_s3_class(res, "synth_result")
  expect_s3_class(res$privacy, "dp_accounting")
  syn <- as.data.frame(res)
  expect_named(syn, names(df))
  expect_gt(nrow(syn), 0)
  expect_s3_class(syn$sex, "factor")
  expect_true(all(syn$age >= 20 & syn$age <= 100))         # respects public bounds
  expect_true(all(syn$sbp >= 60 & syn$sbp <= 220))
})

test_that("DP synthesis is reproducible for a fixed seed", {
  df <- flat_data()
  dp <- dp_control(epsilon = 1, bounds = list(age = c(20, 100), sbp = c(60, 220)))
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 42))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 42))
  expect_equal(a, b)
})

test_that("independent and tree dependence both synthesise", {
  df <- flat_data()
  bnds <- list(age = c(20, 100), sbp = c(60, 220))
  for (dep in c("independent", "tree")) {
    dp <- dp_control(epsilon = 2, dependence = dep, bounds = bnds)
    res <- synth(df, structure = ~ id, privacy = dp, seed = 3)
    expect_s3_class(res, "synth_result")
    expect_equal(res$privacy$dependence, dep)
  }
})

test_that("the gaussian mechanism produces a valid (eps, delta) release", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, delta = 1e-5, mechanism = "gaussian",
                   bounds = list(age = c(20, 100), sbp = c(60, 220)))
  res <- synth(df, structure = ~ id, privacy = dp, seed = 5)
  expect_equal(res$privacy$mechanism, "gaussian")
  expect_gt(res$privacy$rho, 0)
})

test_that("m > 1 draws several datasets from one budget", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, bounds = list(age = c(20, 100), sbp = c(60, 220)))
  res <- synth(df, structure = ~ id, privacy = dp, m = 3, seed = 7)
  expect_equal(res$m, 3L)
  expect_length(res$syn, 3L)
  # Same fitted model, independent draws: the accounting is a single spend.
  expect_equal(res$privacy$epsilon, 2)
})

test_that("large epsilon and n preserves marginals (utility sanity)", {
  set.seed(99)
  df <- flat_data(n = 4000)
  dp <- dp_control(epsilon = 60, dependence = "independent",
                   bounds = list(age = c(20, 100), sbp = c(60, 220)))
  syn <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 8))
  expect_equal(mean(syn$age), mean(df$age), tolerance = 3)
  expect_equal(mean(syn$sex == "F"), mean(df$sex == "F"), tolerance = 0.05)
})

# ---- rigorous discretisation domain ---------------------------------------

test_that("the default domain='dp' estimates edges without warning, epsilon exact", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, delta = 1e-6, mechanism = "gaussian")
  expect_equal(dp$domain, "dp")
  expect_silent(res <- synth(df, structure = ~ id, privacy = dp, seed = 1))
  # No public bounds -> age & sbp are DP-estimated; the total spend is still eps.
  expect_equal(res$privacy$epsilon, 2)
  expect_setequal(res$privacy$domain$vars, c("age", "sbp"))
  expect_equal(res$privacy$domain$frac, 0.1)
  expect_gt(res$privacy$domain$eps_per_query, 0)
})

test_that("domain='public' requires bounds for every numeric variable", {
  df <- flat_data()
  expect_error(
    synth(df, structure = ~ id, privacy = dp_control(epsilon = 2, domain = "public"),
          seed = 1),
    "requires .bounds. for every numeric")
  # With full bounds it spends nothing on the domain.
  dp <- dp_control(epsilon = 2, domain = "public",
                   bounds = list(age = c(20, 100), sbp = c(60, 220)))
  res <- synth(df, structure = ~ id, privacy = dp, seed = 1)
  expect_length(res$privacy$domain$vars, 0L)
  expect_equal(res$privacy$domain$frac, 0)
})

test_that("domain='data' keeps the legacy warned, unaccounted behaviour", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, domain = "data")
  expect_warning(synth(df, structure = ~ id, privacy = dp, seed = 1),
                 "derived from the data")
})

test_that("rigorous modes refuse a bare character column", {
  df <- flat_data()
  df$grp <- sample(c("a", "b"), nrow(df), TRUE)             # character, not factor
  expect_error(
    synth(df, structure = ~ id, privacy = dp_control(epsilon = 2), seed = 1),
    "public category set")
  # Converting to a factor (public levels) makes it acceptable.
  df$grp <- factor(df$grp)
  expect_silent(synth(df, structure = ~ id,
                      privacy = dp_control(epsilon = 2, delta = 1e-6,
                                           mechanism = "gaussian"), seed = 1))
})

test_that("DP domain estimation is reproducible for a fixed seed", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, delta = 1e-6, mechanism = "gaussian")  # estimates edges
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 11))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 11))
  expect_equal(a, b)
})

# ---- domain-estimation math ------------------------------------------------

test_that("dp_exp_quantile is DP-smoothed but tracks the target quantile", {
  set.seed(2)
  x <- rnorm(5000, 50, 10)
  # Large eps_q -> the exponential mechanism concentrates on the true quantile.
  lo <- flexsynth:::dp_exp_quantile(x, 0.05, eps_q = 50, cap = 1)
  hi <- flexsynth:::dp_exp_quantile(x, 0.95, eps_q = 50, cap = 1)
  q <- stats::quantile(x, c(0.05, 0.95))
  expect_lt(abs(lo - q[[1]]), 1)
  expect_lt(abs(hi - q[[2]]), 1)
  expect_true(lo >= min(x) && hi <= max(x))                # outputs a data value
})

test_that("dp_estimate_bounds returns a valid increasing range", {
  b <- flexsynth:::dp_estimate_bounds(c(3, 3, 3, 3), eps_q = 5, cap = 1)  # constant
  expect_true(b[1] < b[2])
  b2 <- flexsynth:::dp_estimate_bounds(rnorm(200), eps_q = 5, cap = 1)
  expect_true(b2[1] < b2[2])
})

test_that("the budget split composes exactly to the total budget", {
  # Laplace: marginal eps + domain eps = total eps.
  dpl <- dp_control(epsilon = 4, domain_frac = 0.25)
  n_dom <- 4L                                              # 2 vars * 2 quantiles
  eps_q <- flexsynth:::dp_quantile_eps(dpl, n_dom, dpl$domain_frac)
  cl <- flexsynth:::dp_calibrate(dpl, n_marginals = 6L, cap = 2, budget_frac = 0.75)
  eps_marg <- (6 * 2) / cl$scale                          # invert scale = m*cap/eps
  expect_equal(eps_marg + n_dom * eps_q, 4)               # exact
  # Gaussian: marginal rho + domain rho = total rho (conservative eps^2/2 bound).
  dpg <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                    domain_frac = 0.25)
  eps_qg <- flexsynth:::dp_quantile_eps(dpg, n_dom, dpg$domain_frac)
  cg <- flexsynth:::dp_calibrate(dpg, n_marginals = 6L, cap = 2, budget_frac = 0.75)
  rho_total <- flexsynth:::zcdp_rho_for(4, 1e-6)
  rho_dom <- n_dom * (eps_qg^2 / 2)
  expect_equal(cg$rho_marginals + rho_dom, rho_total)     # exact
})

test_that("constraints are refused under DP", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, bounds = list(age = c(20, 100), sbp = c(60, 220)))
  expect_error(
    synth(df, structure = ~ id, privacy = dp,
          constraints = rule(sbp > 0), seed = 1),
    "not supported under differentially private")
})

test_that("linked DP is refused with a clear message", {
  patients <- data.frame(id = 1:10, sex = sample(c("F", "M"), 10, TRUE),
                         stringsAsFactors = FALSE)
  expect_error(
    synth_linked(
      tables = list(patients = patients),
      structures = list(patients = ~ id),
      keys = list(patients = "id"),
      privacy = dp_control(epsilon = 1)),
    "not available for linked")
})

test_that("dp_accounting prints its guarantee", {
  df <- flat_data()
  dp <- dp_control(epsilon = 2, bounds = list(age = c(20, 100), sbp = c(60, 220)))
  res <- synth(df, structure = ~ id, privacy = dp, seed = 1)
  expect_output(print(res$privacy), "Track B")
  expect_output(print(res), "differentially private")
})
