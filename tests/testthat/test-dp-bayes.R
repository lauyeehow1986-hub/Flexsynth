# PrivBayes degree-k Bayesian networks (GreedyBayes) on top of the flat Track B
# engine. Covers: degree defaults to 1 and leaves the Chow-Liu path untouched;
# validation of the new argument; a degree-2 network recovering a two-parent
# interaction a degree-1 tree structurally cannot; exact budget composition
# (both mechanisms); the bounded-fan-in DAG structure (<= degree parents, valid
# ancestral order, full coverage); that the exponential-mechanism selection
# spends its slice; reproducibility; gates on longitudinal / linked releases;
# and the accounting print.

# A two-parent interaction that is (near-)invisible to any single pairwise
# marginal but strong in the joint: a driver H pulls A and B in first (high MI),
# while C is a noisy XOR of A and B, so C carries almost no one-parent signal and
# is added last -- its only good explanation is the {A, B} pair. A degree-1 tree
# can give C one parent and misses it; a degree-2 network takes {A, B} and
# recovers it.
mk_xor <- function(n = 6000, seed = 1) {
  set.seed(seed)
  h <- sample(0:1, n, TRUE)
  a <- ifelse(runif(n) < 0.95, h, 1 - h)
  b <- ifelse(runif(n) < 0.85, h, 1 - h)
  x <- as.integer(xor(a == 1, b == 1))          # A XOR B
  cc <- ifelse(runif(n) < 0.9, x, 1 - x)        # noisy parity of A, B
  data.frame(id = seq_len(n), h = factor(h), a = factor(a),
             b = factor(b), c = factor(cc))
}
# Gap between P(C = 1 | A xor B = 1) and P(C = 1 | A xor B = 0). ~0.8 in the real
# data; near 0 for a model that cannot see the {A, B} interaction.
xor_gap <- function(d) {
  bit <- function(z) as.integer(as.character(z))
  x <- as.integer(xor(bit(d$a) == 1, bit(d$b) == 1))
  cc <- bit(d$c)
  mean(cc[x == 1]) - mean(cc[x == 0])
}

# ---- argument + back-compat ------------------------------------------------

test_that("degree defaults to 1 and leaves the Chow-Liu tree path byte-identical", {
  expect_identical(dp_control(epsilon = 1, delta = 1e-6,
                              mechanism = "gaussian")$degree, 1L)
  df <- mk_xor(1500, 3)
  dp0 <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian")
  dp1 <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian", degree = 1)
  r0 <- as.data.frame(synth(df, structure = ~ id, privacy = dp0, seed = 5))
  r1 <- as.data.frame(synth(df, structure = ~ id, privacy = dp1, seed = 5))
  expect_identical(r0, r1)
})

test_that("dp_control validates degree", {
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          degree = 0), "degree")
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          degree = 2.5), "degree")
  # degree > 1 is a Bayesian-network structure over the tree model.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          dependence = "independent", degree = 2), "degree")
  # It is an alternative to the two structure-learning knobs.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          degree = 2, structure_frac = 0.3), "degree|structure")
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          degree = 2, select = "adaptive"), "degree|adaptive")
  # A fat family warns about cell blow-up (bins^(degree + 1)).
  expect_warning(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                            degree = 3, bins = 40), "cell")
})

# ---- the interaction a tree cannot reach -----------------------------------

test_that("a degree-2 network recovers a two-parent interaction a degree-1 tree misses", {
  df <- mk_xor(8000, 4)
  real <- xor_gap(df)
  mk <- function(deg) {
    dp <- dp_control(epsilon = 12, delta = 1e-6, mechanism = "gaussian",
                     degree = deg)
    as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 7))
  }
  g1 <- xor_gap(mk(1))       # Chow-Liu tree: one parent for C, misses it
  g2 <- xor_gap(mk(2))       # degree-2 network: takes {A, B}, recovers it
  expect_gt(real, 0.6)
  expect_lt(abs(g1), 0.3)
  expect_gt(g2, 0.45)
  expect_gt(g2 - g1, 0.25)
})

# ---- exact budget composition ----------------------------------------------

test_that("degree-k gaussian budget composes exactly (measurement + selection)", {
  df <- mk_xor(2000, 6)
  d <- 4L; sf <- 0.3
  dp <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                   degree = 2, select_frac = sf)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 2)$privacy
  by <- ac$bayes
  rho_total <- flexsynth:::zcdp_rho_for(5, 1e-6)
  expect_equal(by$rho_meas, (1 - sf) * rho_total, tolerance = 1e-8)
  expect_equal(by$rho_sel, sf * rho_total, tolerance = 1e-8)
  expect_equal(by$rho_meas + by$rho_sel, rho_total, tolerance = 1e-8)
  # The reported per-cell noise is exactly the sd calibrated from the measurement
  # slice over the 2d - 1 measured marginals at sensitivity cap = 1.
  n_meas <- 2L * d - 1L
  sigma <- 1 * sqrt(n_meas / (2 * by$rho_meas))
  expect_equal(by$meas_noise, sigma, tolerance = 1e-8)
})

test_that("degree-k laplace budget composes exactly (measurement + selection)", {
  df <- mk_xor(2000, 8)
  d <- 4L; sf <- 0.25; eps <- 6
  dp <- dp_control(epsilon = eps, mechanism = "laplace",
                   degree = 2, select_frac = sf)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 3)$privacy
  by <- ac$bayes
  expect_equal(by$eps_meas, (1 - sf) * eps, tolerance = 1e-8)
  expect_equal(by$eps_sel, sf * eps, tolerance = 1e-8)
  expect_equal(by$eps_meas + by$eps_sel, eps, tolerance = 1e-8)
  n_meas <- 2L * d - 1L
  scale <- n_meas * 1 / by$eps_meas
  expect_equal(by$meas_noise, scale, tolerance = 1e-8)
})

# ---- DAG structure ---------------------------------------------------------

test_that("the fitted network is a bounded-fan-in DAG covering every variable", {
  df <- mk_xor(3000, 5)
  vars <- setdiff(names(df), "id")
  d <- length(vars)                              # 4
  degree <- 2L; cap <- 1L
  dp <- dp_control(epsilon = 10, delta = 1e-6, mechanism = "gaussian",
                   degree = degree)
  n_meas <- 2L * d - 1L
  calib <- flexsynth:::dp_calibrate(dp, n_meas, cap,
                                    budget_frac = 1 - dp$select_frac)
  sel_eps <- flexsynth:::dp_select_eps(dp, d - 1L, dp$select_frac)
  dom <- flexsynth:::dp_build_domain(df[vars], vars, dp, NULL)
  nbins <- vapply(vars, function(v) dom[[v]]$nbin, integer(1))
  codes <- stats::setNames(
    lapply(vars, function(v) flexsynth:::dp_encode(dom[[v]], df[[v]])), vars)
  m <- flexsynth:::dp_fit_model_bayes(codes, nbins, dp, calib, degree,
                                      sel_eps, cap)
  expect_identical(m$kind, "bayes")
  expect_equal(m$degree, 2L)
  # One root (no parents) + d - 1 families.
  expect_equal(length(m$cliques), d)
  roots <- vapply(m$cliques, function(cl) length(cl$sep) == 0L, logical(1))
  expect_equal(sum(roots), 1L)
  # Every non-root node has between 1 and `degree` parents.
  fam <- m$cliques[-which(roots)]
  np <- vapply(fam, function(cl) length(cl$sep), integer(1))
  expect_true(all(np >= 1L & np <= degree))
  # Valid ancestral order: each node's parents appear strictly earlier.
  order_added <- vapply(m$cliques, function(cl) cl$new, integer(1))
  seen <- integer(0)
  for (cl in m$cliques) {
    expect_true(all(cl$sep %in% seen))           # parents already generated
    seen <- c(seen, cl$new)
  }
  expect_setequal(order_added, seq_len(d))        # covers all variables once
})

test_that("degree-k selection actually spends its exponential-mechanism slice", {
  # With select_frac -> 0 the EM picks are near-uniform; with a healthy slice the
  # picks are guided, so the recovered interaction is stronger. This shows the
  # selection budget is doing real work.
  df <- mk_xor(8000, 4)
  weak <- dp_control(epsilon = 12, delta = 1e-6, mechanism = "gaussian",
                     degree = 2, select_frac = 0.02)
  strong <- dp_control(epsilon = 12, delta = 1e-6, mechanism = "gaussian",
                       degree = 2, select_frac = 0.3)
  gw <- xor_gap(as.data.frame(synth(df, structure = ~ id, privacy = weak, seed = 7)))
  gs <- xor_gap(as.data.frame(synth(df, structure = ~ id, privacy = strong, seed = 7)))
  expect_gt(gs, gw)
})

test_that("degree-k synthesis is reproducible", {
  df <- mk_xor(1500, 13)
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian", degree = 2)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 99))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 99))
  expect_identical(a, b)
})

# ---- gates -----------------------------------------------------------------

test_that("degree > 1 is refused on longitudinal and linked DP releases", {
  df <- data.frame(id = rep(1:200, each = 2),
                   visit = rep(1:2, times = 200),
                   x = factor(sample(0:1, 400, TRUE)),
                   y = factor(sample(0:1, 400, TRUE)))
  dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                   degree = 2, max_rows_per_person = 2)
  expect_error(synth(df, structure = ~ id / visit, privacy = dp, seed = 1),
               "degree")
})

# ---- print -----------------------------------------------------------------

test_that("prints surface the degree-k network", {
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian", degree = 3)
  expect_output(print(dp), "degree")

  df <- mk_xor(1000, 15)
  dp2 <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian", degree = 2)
  ac <- synth(df, structure = ~ id, privacy = dp2, seed = 1)$privacy
  expect_output(print(ac), "[Bb]ayesian network")
  expect_output(print(ac), "degree 2")
})
