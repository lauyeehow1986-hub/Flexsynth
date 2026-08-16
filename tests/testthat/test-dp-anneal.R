# AIM-style budget annealing (data-adaptive round schedule) on top of the
# adaptive junction-tree selector, plus treewidth >= 3. Covers backward-compat,
# exact composition over a *variable* round count (both mechanisms), the
# spanning invariant, that refinement actually fires when budget is plentiful,
# that re-measurement sharpens a clique, treewidth >= 3 structure + behaviour +
# the cell-blowup warning, gating and prints.

# A three-way generator (cc depends on the JOINT a == b, pairwise-independent of
# both) reused from the adaptive tests; here mainly to exercise annealing on a
# small variable count where surplus budget drives refinement rounds.
mk_three_way <- function(n = 4000, seed = 1) {
  set.seed(seed)
  a <- sample(1:3, n, TRUE)
  b <- ifelse(runif(n) < 0.8, a, sample(1:3, n, TRUE))
  cc <- ifelse(a == b, sample(1:2, n, TRUE, c(.85, .15)),
                       sample(1:2, n, TRUE, c(.15, .85)))
  data.frame(id = seq_len(n), a = factor(a), b = factor(b), cc = factor(cc))
}

# A four-way generator: dd depends on the 3-bit XOR parity of independent fair
# bits a, b, e. Any two of {a, b, e} leave the third fair, so the parity - and
# hence dd's dependence - is invisible to every pairwise and every three-way
# marginal. Only the full four-way clique {a, b, e, dd} (treewidth 3) can hold it.
mk_four_way <- function(n = 8000, seed = 1) {
  set.seed(seed)
  a <- sample(0:1, n, TRUE)
  b <- sample(0:1, n, TRUE)
  e <- sample(0:1, n, TRUE)
  par <- (a + b + e) %% 2
  dd <- ifelse(par == 0, sample(1:2, n, TRUE, c(.9, .1)),
                         sample(1:2, n, TRUE, c(.1, .9)))
  data.frame(id = seq_len(n), a = factor(a), b = factor(b),
             e = factor(e), dd = factor(dd))
}
four_way_gap <- function(d) {
  bit <- function(x) as.integer(as.character(x))
  par <- (bit(d$a) + bit(d$b) + bit(d$e)) %% 2
  lv <- levels(d$dd)[1]
  mean(d$dd[par == 0] == lv) - mean(d$dd[par == 1] == lv)
}


test_that("anneal defaults to FALSE and leaves the fixed-schedule adaptive path untouched", {
  expect_false(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          select = "adaptive")$anneal)
  df <- mk_three_way(1500, 3)
  dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 2)
  r1 <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 5))
  r2 <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 5))
  expect_identical(r1, r2)
  # No annealing record on the non-annealed adaptive path.
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 5)$privacy
  expect_false(isTRUE(ac$adaptive$anneal))
})

test_that("annealed release covers every variable (spanning junction tree)", {
  df <- mk_three_way(2000, 8)
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 2, anneal = TRUE)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  expect_named(a, names(df))
  expect_false(anyNA(a))
  expect_equal(a$id, seq_len(nrow(a)))
})

test_that("annealed gaussian budget composes exactly over the variable round count", {
  df <- mk_three_way(2000, 4)
  sf <- 0.3
  dp <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 1, select_frac = sf,
                   anneal = TRUE)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 2)$privacy
  rho_total <- flexsynth:::zcdp_rho_for(5, 1e-6)
  ad <- ac$adaptive
  # The engine reports the exact per-round zCDP it spent on measurement and on
  # selection; both must sum to the pre-committed marginal split.
  expect_equal(sum(ad$rho_meas_rounds), (1 - sf) * rho_total, tolerance = 1e-6)
  expect_equal(sum(ad$rho_sel_rounds), sf * rho_total, tolerance = 1e-6)
  expect_equal(sum(ad$rho_meas_rounds) + sum(ad$rho_sel_rounds),
               rho_total, tolerance = 1e-6)
})

test_that("annealed laplace budget composes exactly over the variable round count", {
  df <- mk_three_way(2000, 6)
  sf <- 0.25; eps <- 6
  dp <- dp_control(epsilon = eps, mechanism = "laplace",
                   select = "adaptive", treewidth = 1, select_frac = sf,
                   anneal = TRUE)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 3)$privacy
  ad <- ac$adaptive
  expect_equal(sum(ad$eps_meas_rounds), (1 - sf) * eps, tolerance = 1e-6)
  expect_equal(sum(ad$eps_sel_rounds), sf * eps, tolerance = 1e-6)
  expect_equal(sum(ad$eps_meas_rounds) + sum(ad$eps_sel_rounds),
               eps, tolerance = 1e-6)
})

test_that("with surplus budget the schedule adds refinement rounds beyond the spanning tree", {
  # Small variable count + generous budget => spanning tree is cheap, so the
  # annealer should keep spending on re-measurement rounds.
  df <- mk_three_way(3000, 9)
  d <- 3L; w <- 2L
  dp <- dp_control(epsilon = 12, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = w, anneal = TRUE)
  ad <- synth(df, structure = ~ id, privacy = dp, seed = 1)$privacy$adaptive
  expect_true(isTRUE(ad$anneal))
  expect_gte(ad$n_rounds, d - w)               # at least the spanning rounds
  expect_gte(ad$n_refine, 1L)                  # and at least one refinement
  expect_equal(ad$n_rounds, (d - w) + ad$n_refine)
})

test_that("re-measuring a clique lowers its effective per-cell noise (rho adds)", {
  # White-box: combine() two noisy measurements of the same true array.
  set.seed(1)
  truth <- matrix(sample(50:150, 6), 2, 3)
  s1 <- 5; s2 <- 3
  A1 <- truth + matrix(stats::rnorm(6, sd = s1), 2, 3)
  A2 <- truth + matrix(stats::rnorm(6, sd = s2), 2, 3)
  comb <- flexsynth:::dp_combine_gaussian(A1, s1, A2, s2)
  # Inverse-variance weighting => effective sd with 1/s^2 = 1/s1^2 + 1/s2^2.
  s_eff <- sqrt(1 / (1 / s1^2 + 1 / s2^2))
  expect_lt(s_eff, min(s1, s2))
  expect_equal(comb$sd, s_eff, tolerance = 1e-9)
  # The combined estimate is the inverse-variance weighted mean.
  w1 <- 1 / s1^2; w2 <- 1 / s2^2
  expect_equal(comb$value, (w1 * A1 + w2 * A2) / (w1 + w2), tolerance = 1e-9)
})

test_that("the sigma-halving signal test fires when a round's budget is too small", {
  # A high-noise regime (tiny epsilon, many rounds' worth of variables) should
  # trip the annealer's noise test at least once, recorded in the schedule.
  df <- data.frame(id = 1:4000,
                   v1 = factor(sample(1:4, 4000, TRUE)),
                   v2 = factor(sample(1:4, 4000, TRUE)),
                   v3 = factor(sample(1:4, 4000, TRUE)),
                   v4 = factor(sample(1:4, 4000, TRUE)))
  dp <- dp_control(epsilon = 0.3, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 1, anneal = TRUE)
  ad <- synth(df, structure = ~ id, privacy = dp, seed = 7)$privacy$adaptive
  expect_true(isTRUE(ad$anneal))
  expect_gte(ad$n_anneal_steps, 1L)            # sigma halved at least once
})

test_that("annealed synthesis is reproducible", {
  df <- mk_three_way(1500, 13)
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 2, anneal = TRUE)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 99))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 99))
  expect_identical(a, b)
})

# ---- treewidth >= 3 --------------------------------------------------------

test_that("dp_control now accepts treewidth >= 3 and warns on cell blow-up", {
  expect_silent(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                           select = "adaptive", treewidth = 3, bins = 3))
  # A fat clique (bins^(w+1) huge) warns the user that noise will dominate.
  expect_warning(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                            select = "adaptive", treewidth = 3, bins = 40),
                 "cell")
})

test_that("treewidth 3 builds 4-variable cliques covering every variable", {
  df <- mk_four_way(2000, 5)
  vars <- setdiff(names(df), "id")
  d <- length(vars)                            # 4
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 3)
  # Reach into the fitted model via the internal fitter for structure checks.
  w_eff <- min(3L, d - 1L)
  cap <- 1L
  n_marg <- d + (d - w_eff)
  calib <- flexsynth:::dp_calibrate(dp, n_marg, cap, budget_frac = 1 - dp$select_frac)
  sel_eps <- flexsynth:::dp_select_eps(dp, d - w_eff, dp$select_frac)
  dom <- flexsynth:::dp_build_domain(df[vars], vars, dp, NULL)
  nbins <- vapply(vars, function(v) dom[[v]]$nbin, integer(1))
  codes <- stats::setNames(
    lapply(vars, function(v) flexsynth:::dp_encode(dom[[v]], df[[v]])), vars)
  m <- flexsynth:::dp_fit_model_adaptive(codes, nbins, dp, calib, w_eff, sel_eps, cap)
  expect_equal(m$treewidth, 3L)
  expect_equal(length(m$cliques), d - w_eff)   # one 4-clique for d = 4
  expect_setequal(m$cliques[[1L]]$vars, seq_len(d))
  expect_true(all(vapply(m$cliques, function(cl) length(cl$vars), integer(1)) == 4L))
})

test_that("treewidth 3 recovers a four-way interaction that treewidth 2 misses", {
  df <- mk_four_way(8000, 11)
  real <- four_way_gap(df)
  mk <- function(w) {
    dp <- dp_control(epsilon = 12, delta = 1e-6, mechanism = "gaussian",
                     select = "adaptive", treewidth = w)
    as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 21))
  }
  g2 <- four_way_gap(mk(2))
  g3 <- four_way_gap(mk(3))
  expect_gt(real, 0.5)
  expect_lt(abs(g2), 0.3)                       # pairwise/triangle miss it
  expect_gt(g3, 0.4)                            # the 4-clique recovers it
  expect_gt(g3 - g2, 0.25)
})

# ---- prints ----------------------------------------------------------------

test_that("prints surface the annealing controls and treewidth 3", {
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 3, anneal = TRUE, bins = 3)
  expect_output(print(dp), "annealed")
  expect_output(print(dp), "treewidth 3")

  df <- mk_three_way(1000, 15)
  dp2 <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                    select = "adaptive", treewidth = 2, anneal = TRUE)
  ac <- synth(df, structure = ~ id, privacy = dp2, seed = 1)$privacy
  expect_output(print(ac), "annealed")
  expect_output(print(ac), "round\\(s\\)")
})

test_that("anneal is rejected without adaptive selection", {
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          anneal = TRUE),
               "anneal")
})
