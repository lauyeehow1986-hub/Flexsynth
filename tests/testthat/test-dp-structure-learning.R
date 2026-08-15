# Budget-efficient structure learning for the flat DP tree release:
# dp_control(structure_frac = f) selects the Chow-Liu tree from a cheap all-pairs
# scan (fraction f of the marginal budget), then concentrates the remaining
# 1 - f on re-measuring only the chosen edges. Both passes must compose into the
# same exact (epsilon, delta) budget.

# A chained-correlation numeric table: v1 -> v2 -> ... so a tree is meaningful.
mk_chain <- function(n = 500L, d = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(0, n, d)
  m[, 1L] <- stats::rnorm(n)
  for (j in 2:d) m[, j] <- 0.9 * m[, j - 1L] + stats::rnorm(n, sd = 0.4)
  df <- as.data.frame(round(m, 2))
  names(df) <- paste0("v", seq_len(d))
  df$id <- seq_len(n)
  df
}
chain_bounds <- function(d) stats::setNames(rep(list(c(-6, 6)), d),
                                            paste0("v", seq_len(d)))

test_that("structure_frac splits a Laplace release into an exact composition", {
  df <- mk_chain(d = 5L)
  bd <- chain_bounds(5L)
  eps <- 2

  naive <- synth(df, structure = ~ id, seed = 7,
                 privacy = dp_control(epsilon = eps, mechanism = "laplace",
                                      bounds = bd, domain = "public"))
  learn <- synth(df, structure = ~ id, seed = 7,
                 privacy = dp_control(epsilon = eps, mechanism = "laplace",
                                      bounds = bd, domain = "public",
                                      structure_frac = 0.25))
  an <- naive$privacy; al <- learn$privacy
  d <- 5L; cap <- 1
  n_pairs <- d * (d - 1L) / 2L          # 10
  n_param <- 2L * d - 1L                # 9

  # Naive: d + C(d,2) = 15 marginals, single scale.
  expect_equal(an$n_marginals, d + n_pairs)
  expect_null(an$learn)
  expect_equal(an$noise, (d + n_pairs) * cap / eps)

  # Learn: reports the two-phase count; parameter noise is the headline.
  expect_equal(al$n_marginals, n_pairs + n_param)     # 19
  expect_false(is.null(al$learn))
  expect_equal(al$learn$n_struct, n_pairs)
  expect_equal(al$learn$n_param, n_param)
  expect_equal(al$noise, n_param * cap / ((1 - 0.25) * eps))            # 6
  expect_equal(al$learn$struct_noise, n_pairs * cap / (0.25 * eps))     # 20
})

test_that("structure_frac splits a Gaussian release into an exact zCDP composition", {
  df <- mk_chain(d = 5L)
  bd <- chain_bounds(5L)
  eps <- 2; delta <- 1e-6
  rho <- flexsynth:::zcdp_rho_for(eps, delta)
  cap <- 1; d <- 5L
  n_pairs <- d * (d - 1L) / 2L; n_param <- 2L * d - 1L

  learn <- synth(df, structure = ~ id, seed = 7,
                 privacy = dp_control(epsilon = eps, delta = delta,
                                      mechanism = "gaussian", bounds = bd,
                                      domain = "public", structure_frac = 0.25))
  al <- learn$privacy
  # Parameter phase gets (1 - f) of rho over n_param queries; scan f over n_pairs.
  expect_equal(al$noise, cap * sqrt(n_param / (2 * (1 - 0.25) * rho)))
  expect_equal(al$learn$struct_noise, cap * sqrt(n_pairs / (2 * 0.25 * rho)))
})

test_that("the domain slice, structure scan and parameters compose to the full budget", {
  # No public bounds -> DP domain estimation takes domain_frac; the remaining
  # marg_frac is split by structure_frac. Scales must reflect marg_frac exactly.
  df <- mk_chain(d = 5L)
  eps <- 3; f <- 0.3; dfrac <- 0.1
  marg <- 1 - dfrac
  d <- 5L; cap <- 1
  n_param <- 2L * d - 1L; n_pairs <- d * (d - 1L) / 2L

  learn <- synth(df, structure = ~ id, seed = 3,
                 privacy = dp_control(epsilon = eps, mechanism = "laplace",
                                      domain = "dp", domain_frac = dfrac,
                                      structure_frac = f))
  al <- learn$privacy
  expect_equal(al$noise, n_param * cap / (marg * (1 - f) * eps))
  expect_equal(al$learn$struct_noise, n_pairs * cap / (marg * f * eps))
  expect_true(length(al$domain$vars) == d)             # all estimated under DP
  expect_equal(al$domain$frac, dfrac)
})

test_that("structure_frac is inert for < 3 variables and for the independent model", {
  # d = 2: the tree has a single edge = the only pair; nothing to learn.
  df2 <- mk_chain(d = 2L)
  r2 <- synth(df2, structure = ~ id, seed = 1,
              privacy = dp_control(epsilon = 2, mechanism = "laplace",
                                   bounds = chain_bounds(2L), domain = "public",
                                   structure_frac = 0.25))
  expect_null(r2$privacy$learn)
  expect_equal(r2$privacy$n_marginals, 2L + 1L)        # 2 one-way + 1 pair

  # independent: no structure to learn.
  df5 <- mk_chain(d = 5L)
  r5 <- synth(df5, structure = ~ id, seed = 1,
              privacy = dp_control(epsilon = 2, mechanism = "laplace",
                                   bounds = chain_bounds(5L), domain = "public",
                                   dependence = "independent",
                                   structure_frac = 0.25))
  expect_null(r5$privacy$learn)
  expect_equal(r5$privacy$n_marginals, 5L)             # one-ways only
})

test_that("concentrating budget on the tree lowers the parameter noise (d >= 5)", {
  df <- mk_chain(d = 6L)
  bd <- chain_bounds(6L)
  naive <- synth(df, structure = ~ id, seed = 2,
                 privacy = dp_control(epsilon = 2, mechanism = "laplace",
                                      bounds = bd, domain = "public"))
  learn <- synth(df, structure = ~ id, seed = 2,
                 privacy = dp_control(epsilon = 2, mechanism = "laplace",
                                      bounds = bd, domain = "public",
                                      structure_frac = 0.3))
  expect_lt(learn$privacy$noise, naive$privacy$noise)
})

test_that("structure learning is reproducible and preserves shape", {
  df <- mk_chain(d = 5L)
  dp <- dp_control(epsilon = 2, mechanism = "laplace",
                   bounds = chain_bounds(5L), domain = "public",
                   structure_frac = 0.25)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 42))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 42))
  expect_identical(a, b)
  expect_identical(names(a), names(df))
  expect_equal(anyDuplicated(a$id), 0L)
})

test_that("budget-efficient learning improves correlation fidelity on average", {
  df <- mk_chain(n = 500L, d = 8L, seed = 1L)
  vs <- paste0("v", 1:8)
  bd <- chain_bounds(8L)
  truec <- stats::cor(df[vs])
  cor_err <- function(rr) {
    s <- as.data.frame(rr)[vs]
    mean(abs(stats::cor(s) - truec)[upper.tri(truec)])
  }
  dn <- dp_control(epsilon = 2, mechanism = "laplace", bounds = bd,
                   domain = "public", bins = 8L)
  dl <- dp_control(epsilon = 2, mechanism = "laplace", bounds = bd,
                   domain = "public", bins = 8L, structure_frac = 0.25)
  seeds <- 1:20
  en <- vapply(seeds, function(s) cor_err(synth(df, ~ id, privacy = dn, seed = s)),
               numeric(1))
  el <- vapply(seeds, function(s) cor_err(synth(df, ~ id, privacy = dl, seed = s)),
               numeric(1))
  expect_lt(mean(el), mean(en))
})

test_that("dp_control validates structure_frac", {
  expect_silent(dp_control(epsilon = 1, structure_frac = NULL))
  expect_silent(dp_control(epsilon = 1, structure_frac = 0.3))
  expect_error(dp_control(epsilon = 1, structure_frac = 0),   "structure_frac")
  expect_error(dp_control(epsilon = 1, structure_frac = 1),   "structure_frac")
  expect_error(dp_control(epsilon = 1, structure_frac = -0.1), "structure_frac")
  expect_error(dp_control(epsilon = 1, structure_frac = 1.5), "structure_frac")
  expect_error(dp_control(epsilon = 1, structure_frac = c(0.2, 0.3)),
               "structure_frac")
})

test_that("structure learning is announced in the printed records", {
  dp <- dp_control(epsilon = 2, mechanism = "laplace",
                   bounds = chain_bounds(5L), domain = "public",
                   structure_frac = 0.25)
  expect_output(print(dp), "budget-efficient")
  acct <- synth(mk_chain(d = 5L), structure = ~ id, privacy = dp, seed = 1)$privacy
  expect_output(print(acct), "structure : budget-efficient")
  expect_output(print(acct), "pairwise scans")
})
