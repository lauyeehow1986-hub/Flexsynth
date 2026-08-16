# AIM-style adaptive marginal selection (flat Track B). Covers backward-compat,
# exact budget composition (both mechanisms), junction-tree / treewidth
# invariants, the treewidth-2-beats-treewidth-1 behavioural win, gating, and
# prints.

# A distribution whose third variable depends on the JOINT (a == b) but is
# pairwise-independent of both a and b: a tree cannot capture it, a triangle can.
mk_three_way <- function(n = 4000, seed = 1) {
  set.seed(seed)
  a <- sample(1:3, n, TRUE)
  b <- ifelse(runif(n) < 0.8, a, sample(1:3, n, TRUE))          # b ~ a (pairwise)
  cc <- ifelse(a == b, sample(1:2, n, TRUE, c(.85, .15)),
                       sample(1:2, n, TRUE, c(.15, .85)))         # depends on (a==b)
  data.frame(id = seq_len(n), a = factor(a), b = factor(b), cc = factor(cc))
}
three_way_bounds <- NULL  # all-factor, no numeric bounds needed

# Direct fit of the internal adaptive model, for structural assertions.
fit_adaptive <- function(df, w, eps = 6, select_frac = 0.25) {
  dp <- dp_control(epsilon = eps, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = w, select_frac = select_frac)
  vars <- setdiff(names(df), "id")
  d <- length(vars)
  w_eff <- min(w, d - 1L)
  n_cliques <- d - w_eff
  n_marg <- d + n_cliques
  cap <- 1L
  calib <- flexsynth:::dp_calibrate(dp, n_marg, cap, budget_frac = 1 - select_frac)
  sel_eps <- flexsynth:::dp_select_eps(dp, n_cliques, select_frac)
  dom <- flexsynth:::dp_build_domain(df[vars], vars, dp, NULL)
  nbins <- vapply(vars, function(v) dom[[v]]$nbin, integer(1))
  codes <- stats::setNames(lapply(vars, function(v) flexsynth:::dp_encode(dom[[v]], df[[v]])), vars)
  flexsynth:::dp_fit_model_adaptive(codes, nbins, dp, calib, w_eff, sel_eps, cap)
}

three_way_gap <- function(d) {
  eq <- d$a == d$b
  lv <- levels(d$cc)[1]
  mean(d$cc[eq] == lv) - mean(d$cc[!eq] == lv)
}


test_that("select defaults to fixed and leaves the fixed path untouched", {
  expect_equal(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian")$select,
               "fixed")
  df <- mk_three_way(1500, 3)
  dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian")
  r <- synth(df, structure = ~ id, privacy = dp, seed = 5)
  # Fixed tree release: no adaptive record, and the classic d + C(d,2) count.
  expect_null(r$privacy$adaptive)
  d <- 3L
  expect_equal(r$privacy$n_marginals, d + d * (d - 1L) / 2L)
  # Reproducible and unchanged by the presence of the new code path.
  r2 <- synth(df, structure = ~ id, privacy = dp, seed = 5)
  expect_identical(as.data.frame(r), as.data.frame(r2))
})

test_that("adaptive n_marginals = d one-ways + (d - treewidth) cliques", {
  df <- mk_three_way(1200, 2)
  d <- 3L
  for (w in 1:2) {
    dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                     select = "adaptive", treewidth = w)
    r <- synth(df, structure = ~ id, privacy = dp, seed = 1)
    w_eff <- min(w, d - 1L)
    expect_equal(r$privacy$adaptive$treewidth, w_eff)
    expect_equal(r$privacy$adaptive$n_cliques, d - w_eff)
    expect_equal(r$privacy$n_marginals, d + (d - w_eff))
  }
})

test_that("gaussian budget composes exactly: measurement + selection = total rho", {
  df <- mk_three_way(1500, 4)
  d <- 3L; cap <- 1L; sf <- 0.3
  dp <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 2, select_frac = sf)
  # all-factor data => no numeric domain estimation => marg_frac = 1
  r <- synth(df, structure = ~ id, privacy = dp, seed = 2)
  ac <- r$privacy
  rho_total <- flexsynth:::zcdp_rho_for(5, 1e-6)
  # Measurement rho recovered from the reported per-cell Gaussian sd.
  rho_meas <- ac$n_marginals * cap^2 / (2 * ac$noise^2)
  # Selection rho: each of n_cliques rounds is (sel_eps^2 / 2)-zCDP.
  rho_sel <- ac$adaptive$n_cliques * ac$adaptive$select_eps^2 / 2
  expect_equal(rho_meas, (1 - sf) * rho_total, tolerance = 1e-6)
  expect_equal(rho_sel, sf * rho_total, tolerance = 1e-6)
  expect_equal(rho_meas + rho_sel, rho_total, tolerance = 1e-6)
})

test_that("laplace budget composes exactly: measurement eps + selection eps = total", {
  df <- mk_three_way(1500, 6)
  cap <- 1L; sf <- 0.25; eps <- 6
  dp <- dp_control(epsilon = eps, mechanism = "laplace",
                   select = "adaptive", treewidth = 2, select_frac = sf)
  r <- synth(df, structure = ~ id, privacy = dp, seed = 3)
  ac <- r$privacy
  eps_meas <- ac$n_marginals * cap / ac$noise      # from scale = n*cap/eps_meas
  eps_sel  <- ac$adaptive$n_cliques * ac$adaptive$select_eps
  expect_equal(eps_meas, (1 - sf) * eps, tolerance = 1e-6)
  expect_equal(eps_sel, sf * eps, tolerance = 1e-6)
  expect_equal(eps_meas + eps_sel, eps, tolerance = 1e-6)
})

test_that("treewidth 1 builds a spanning tree; treewidth 2 a partial 2-tree", {
  df <- mk_three_way(1500, 7)
  d <- 3L
  m1 <- fit_adaptive(df, 1)
  expect_equal(length(m1$cliques), d - 1L)             # spanning tree edges
  expect_true(all(vapply(m1$cliques, function(cl) length(cl$vars), integer(1)) == 2L))
  expect_setequal(unlist(lapply(m1$cliques, `[[`, "vars")), seq_len(d))
  # root has empty separator; every other clique's separator has 1 var.
  expect_length(m1$cliques[[1]]$sep, 0L)
  expect_true(all(vapply(m1$cliques[-1], function(cl) length(cl$sep), integer(1)) == 1L))

  m2 <- fit_adaptive(df, 2)
  expect_equal(length(m2$cliques), d - 2L)             # one triangle for d = 3
  expect_equal(m2$treewidth, 2L)
  expect_true(all(vapply(m2$cliques, function(cl) length(cl$vars), integer(1)) == 3L))
  expect_setequal(m2$cliques[[1]]$vars, seq_len(d))
})

test_that("treewidth is capped at d - 1", {
  df <- data.frame(id = 1:800, x = factor(sample(1:3, 800, TRUE)),
                   y = factor(sample(1:2, 800, TRUE)))
  m <- fit_adaptive(df, 2)                              # d = 2 -> w capped to 1
  expect_equal(m$treewidth, 1L)
  expect_equal(length(m$cliques), 1L)                  # the single pair
})

test_that("treewidth 2 recovers a three-way interaction that treewidth 1 misses", {
  df <- mk_three_way(4000, 11)
  real <- three_way_gap(df)
  mk <- function(w) {
    dp <- dp_control(epsilon = 9, delta = 1e-6, mechanism = "gaussian",
                     select = "adaptive", treewidth = w)
    as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 21))
  }
  g1 <- three_way_gap(mk(1))
  g2 <- three_way_gap(mk(2))
  expect_gt(real, 0.5)
  expect_lt(abs(g1), 0.25)                 # tree sees cc as ~independent
  expect_gt(g2, 0.45)                      # triangle recovers the interaction
  expect_gt(g2 - g1, 0.3)
})

test_that("adaptive synthesis is reproducible and shape-correct", {
  df <- mk_three_way(1000, 13)
  dp <- dp_control(epsilon = 5, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 2)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 99))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 99))
  expect_identical(a, b)
  expect_named(a, names(df))
  expect_equal(a$id, seq_len(nrow(a)))
})

test_that("adaptive is refused on longitudinal and linked releases", {
  long <- data.frame(id = rep(1:40, each = 3), visit = rep(1:3, 40),
                     x = rnorm(120))
  dp <- dp_control(epsilon = 2, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive")
  expect_error(synth(long, structure = ~ id / visit, privacy = dp, seed = 1),
               "flat-table only")

  patients <- data.frame(id = 1:30, sex = factor(sample(c("F", "M"), 30, TRUE)))
  visits <- do.call(rbind, lapply(patients$id, function(p)
    data.frame(id = p, v = 1:2, val = factor(sample(1:3, 2, TRUE)))))
  expect_error(
    synth_linked(tables = list(patients = patients, visits = visits),
                 structures = list(patients = ~ id, visits = ~ id / v),
                 keys = list(patients = "id", visits = c("id", "v")),
                 privacy = dp, seed = 1),
    "flat-table only")
})

test_that("dp_control rejects structure_frac + adaptive and bad treewidth / select_frac", {
  # treewidth > 2 is now supported (see test-dp-anneal.R); the mutual-exclusion
  # and range validations still hold.
  expect_error(dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                          select = "adaptive", structure_frac = 0.3),
               "structure_frac")
  expect_error(dp_control(epsilon = 1, treewidth = 0), "treewidth")
  expect_error(dp_control(epsilon = 1, select_frac = 1.2), "select_frac")
})

test_that("prints surface the adaptive controls", {
  dp <- dp_control(epsilon = 3, delta = 1e-6, mechanism = "gaussian",
                   select = "adaptive", treewidth = 2, select_frac = 0.3)
  expect_output(print(dp), "adaptive AIM-style, treewidth 2")

  df <- mk_three_way(800, 15)
  r <- synth(df, structure = ~ id, privacy = dp, seed = 1)
  expect_output(print(r$privacy), "adaptive junction tree \\(treewidth 2\\)")
  expect_output(print(r$privacy), "selection : adaptive")
  expect_output(print(r$privacy), "adaptively selected cliques")
})
