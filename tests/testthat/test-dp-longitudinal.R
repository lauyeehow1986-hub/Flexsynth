# Track B DP Markov longitudinal synthesis: within-unit temporal structure under
# differential privacy. Tests cover the contribution model, exact budget
# composition of the length + initial + transition release, structure
# regeneration, autocorrelation preservation, and reproducibility.

# Autocorrelated longitudinal cardiac data: `np` patients, several visits each,
# sbp follows an AR(1) within patient so there is real lag-1 structure to keep.
long_data <- function(np = 250, seed = 1) {
  set.seed(seed)
  parts <- lapply(seq_len(np), function(i) {
    nv <- sample(2:4, 1)
    sbp <- numeric(nv); sbp[1] <- rnorm(1, 130, 12)
    for (t in seq_len(nv)[-1L]) sbp[t] <- 0.85 * sbp[t - 1] + 0.15 * 130 + rnorm(1, 0, 5)
    data.frame(id = i, visit = seq_len(nv),
               age = round(rnorm(1, 62, 9)) + seq_len(nv) - 1L,
               sbp = round(sbp),
               sex = sample(c("F", "M"), 1),
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, parts)
  d$sex <- factor(d$sex, levels = c("F", "M"))
  d
}

# lag-1 within-person correlation of a numeric column
lag1_cor <- function(d, col = "sbp") {
  d <- d[order(d$id, d$visit), ]
  prev <- stats::ave(d[[col]], d$id, FUN = function(x) c(NA, utils::head(x, -1)))
  ok <- !is.na(prev)
  if (sum(ok) < 3) return(NA_real_)
  stats::cor(prev[ok], d[[col]][ok])
}

pub <- list(age = c(20, 110), sbp = c(60, 220))

# ---- contribution model ----------------------------------------------------

test_that("dp_truncate_prefix keeps each person's first `cap` rows in order", {
  df <- data.frame(id = c(1, 1, 1, 2, 2, 3), visit = c(1, 2, 3, 1, 2, 1), x = 1:6)
  r <- flexsynth:::dp_truncate_prefix(df, "id", cap = 2L)
  expect_equal(nrow(r$data), 5L)
  expect_equal(r$dropped, 1L)                 # person 1 lost their 3rd row
  expect_equal(r$data$x, c(1L, 2L, 4L, 5L, 6L))  # prefix, contiguous
})

# ---- dispatch / structure regeneration -------------------------------------

test_that("a nesting index triggers the DP Markov engine", {
  d <- long_data()
  dp <- dp_control(epsilon = 6, max_rows_per_person = 4, bounds = pub)
  res <- synth(d, structure = ~ id / visit, privacy = dp, seed = 1)
  expect_s3_class(res, "synth_result")
  expect_s3_class(res$privacy, "dp_accounting")
  expect_false(is.null(res$privacy$longitudinal))
  syn <- as.data.frame(res)
  expect_named(syn, names(d))
  # visit regenerated as a within-person position starting at 1
  expect_equal(min(syn$visit), 1L)
  expect_true(all(syn$visit <= 4L))
  # every person's visits are 1..len with no gaps
  by_id <- split(syn$visit, syn$id)
  ok <- vapply(by_id, function(v) identical(sort(v), seq_along(v)), logical(1))
  expect_true(all(ok))
})

test_that("longitudinal DP needs max_rows_per_person >= 2", {
  d <- long_data()
  expect_error(
    synth(d, structure = ~ id / visit, privacy = dp_control(epsilon = 4), seed = 1),
    "max_rows_per_person >= 2")
})

# ---- exact budget composition ----------------------------------------------

test_that("the length + initial + transition release composes exactly (Laplace)", {
  d <- long_data()
  cap <- 4L; nV <- 3L
  dp <- dp_control(epsilon = 5, mechanism = "laplace",
                   max_rows_per_person = cap, bounds = pub)   # public -> no domain spend
  res <- synth(d, structure = ~ id / visit, privacy = dp, seed = 1)
  n_init <- nV + nV * (nV - 1L) / 2L            # tree: one-way + all pairwise
  total_l1 <- 1 + n_init + nV * (cap - 1L)      # length + initial + transitions
  expect_equal(res$privacy$noise, total_l1 / 5) # scale = total_l1 / eps, marg_frac 1
  expect_equal(res$privacy$longitudinal$n_transitions, nV)
  expect_equal(res$privacy$longitudinal$n_init_marg, n_init)
  expect_equal(res$privacy$n_marginals, 1L + n_init + nV)  # total histograms
})

test_that("the release composes exactly under Gaussian zCDP", {
  d <- long_data()
  cap <- 4L; nV <- 3L
  dp <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                   max_rows_per_person = cap, bounds = pub)
  res <- synth(d, structure = ~ id / visit, privacy = dp, seed = 1)
  n_init <- nV + nV * (nV - 1L) / 2L
  sum_sq <- 1 + n_init + nV * (cap - 1L)^2      # summed squared L2 sensitivities
  rho <- flexsynth:::zcdp_rho_for(5, 1e-6)
  expect_equal(res$privacy$noise, sqrt(sum_sq / (2 * rho)))
  expect_equal(res$privacy$epsilon, 5)          # exact end-to-end guarantee
})

test_that("DP-estimated edges take an accounted slice under longitudinal DP", {
  d <- long_data()
  dp <- dp_control(epsilon = 6, mechanism = "laplace", max_rows_per_person = 4)
  expect_silent(res <- synth(d, structure = ~ id / visit, privacy = dp, seed = 1))
  expect_setequal(res$privacy$domain$vars, c("age", "sbp"))
  expect_equal(res$privacy$domain$frac, 0.1)
  expect_equal(res$privacy$epsilon, 6)          # still exactly eps
})

test_that("dp_calibrate is the flat special case of dp_make_noise", {
  dp <- dp_control(epsilon = 2, mechanism = "laplace")
  a <- flexsynth:::dp_calibrate(dp, n_marginals = 5L, cap = 3)
  b <- flexsynth:::dp_make_noise(dp, total_l1 = 5 * 3, sum_sq = 5 * 3^2)
  expect_equal(a$scale, b$scale)
  dpg <- dp_control(epsilon = 2, delta = 1e-6, mechanism = "gaussian")
  ag <- flexsynth:::dp_calibrate(dpg, n_marginals = 5L, cap = 3)
  bg <- flexsynth:::dp_make_noise(dpg, total_l1 = 5 * 3, sum_sq = 5 * 3^2)
  expect_equal(ag$sigma, bg$sigma)
})

# ---- utility: temporal structure is preserved ------------------------------

test_that("the Markov model keeps within-person autocorrelation", {
  d <- long_data(np = 500)
  expect_gt(lag1_cor(d), 0.5)                   # strong real AR(1)
  dp <- dp_control(epsilon = 40, mechanism = "laplace",  # generous budget
                   max_rows_per_person = 4, bounds = pub)
  syn <- as.data.frame(synth(d, structure = ~ id / visit, privacy = dp, seed = 2))
  # DP attenuates but a clear positive lag-1 correlation survives.
  expect_gt(lag1_cor(syn), 0.2)
})

# ---- housekeeping ----------------------------------------------------------

test_that("longitudinal DP is reproducible and supports m > 1", {
  d <- long_data()
  dp <- dp_control(epsilon = 5, max_rows_per_person = 4, bounds = pub)
  a <- as.data.frame(synth(d, ~ id / visit, privacy = dp, seed = 7))
  b <- as.data.frame(synth(d, ~ id / visit, privacy = dp, seed = 7))
  expect_equal(a, b)

  res <- synth(d, ~ id / visit, privacy = dp, m = 3, seed = 7)
  expect_equal(res$m, 3L)
  expect_length(res$syn, 3L)
  expect_equal(res$privacy$epsilon, 5)          # one spend, three draws
})

test_that("domain='public' still requires bounds for every numeric column", {
  d <- long_data()
  expect_error(
    synth(d, ~ id / visit,
          privacy = dp_control(epsilon = 5, domain = "public",
                               max_rows_per_person = 4),
          seed = 1),
    "requires .bounds. for every numeric")
})

test_that("longitudinal accounting prints the Markov composition", {
  d <- long_data()
  dp <- dp_control(epsilon = 5, max_rows_per_person = 4, bounds = pub)
  res <- synth(d, ~ id / visit, privacy = dp, seed = 1)
  expect_output(print(res$privacy), "DP Markov")
  expect_output(print(res$privacy), "transition")
})
