# Track B: higher-order (own-lag) and cross-variable DP Markov transitions in the
# flat DP longitudinal engine (dp_control(transition_order=, transition_cross=)).

# A longitudinal table with two time-varying, correlated vitals (hr, sbp; sbp
# tracks hr) plus a subject-invariant sex. 3..5 visits per person.
mk_traj <- function(n = 90, seed = 1) {
  set.seed(seed)
  parts <- lapply(seq_len(n), function(i) {
    v   <- sample(3:5, 1L)
    sex <- sample(c("F", "M"), 1L)
    hr  <- round(72 + cumsum(stats::rnorm(v, 0, 6)))
    sbp <- round(0.8 * hr + 60 + stats::rnorm(v, 0, 5))
    data.frame(id = i, visit = seq_len(v),
               sex = factor(sex, levels = c("F", "M")),
               hr = hr, sbp = sbp)
  })
  do.call(rbind, parts)
}

traj_bounds <- list(hr = c(40, 140), sbp = c(60, 240))

dp_traj <- function(order = 1L, cross = 0L, dependence = "tree",
                    mechanism = "laplace", epsilon = 4,
                    delta = if (mechanism == "gaussian") 1e-6 else 0,
                    cap = 5L) {
  dp_control(epsilon = epsilon, delta = delta, mechanism = mechanism,
             dependence = dependence, transition_order = order,
             transition_cross = cross, max_rows_per_person = cap,
             bounds = traj_bounds, domain = "public")
}

# ---- exact budget composition -------------------------------------------------

test_that("higher order gives exact Laplace composition and lower noise", {
  df <- mk_traj()
  cap <- 5; nV <- 3; nT <- 3
  # independent -> n_init_marg = nV exactly (no pairwise term).
  o1 <- synth(df, structure = ~ id / visit,
              privacy = dp_traj(1L, 0L, dependence = "independent"),
              seed = 1)$privacy
  o2 <- synth(df, structure = ~ id / visit,
              privacy = dp_traj(2L, 0L, dependence = "independent"),
              seed = 1)$privacy
  o3 <- synth(df, structure = ~ id / visit,
              privacy = dp_traj(3L, 0L, dependence = "independent"),
              seed = 1)$privacy
  expect_equal(o1$noise, (1 + nV + nT * (cap - 1)) / 4, tolerance = 1e-9)
  expect_equal(o2$noise, (1 + nV + nT * (cap - 2)) / 4, tolerance = 1e-9)
  expect_equal(o3$noise, (1 + nV + nT * (cap - 3)) / 4, tolerance = 1e-9)
  # Splitting a person's transitions into fewer, deeper tuples lowers sensitivity.
  expect_lt(o2$noise, o1$noise)
  expect_lt(o3$noise, o2$noise)
})

test_that("higher order gives exact Gaussian (zCDP) composition", {
  df <- mk_traj()
  cap <- 5; nV <- 3; nT <- 3; ord <- 2
  g <- synth(df, structure = ~ id / visit,
             privacy = dp_traj(ord, 0L, dependence = "independent",
                               mechanism = "gaussian"), seed = 1)$privacy
  sum_sq <- 1 + nV + nT * (cap - ord)^2
  rho <- zcdp_rho_for(4, 1e-6)
  expect_equal(g$noise, sqrt(sum_sq / (2 * rho)), tolerance = 1e-9)
})

test_that("cross-parents are budget-neutral (same sensitivity as no cross)", {
  df <- mk_traj()
  c0 <- synth(df, structure = ~ id / visit,
              privacy = dp_traj(1L, 0L, dependence = "tree"), seed = 1)$privacy
  c1 <- synth(df, structure = ~ id / visit,
              privacy = dp_traj(1L, 1L, dependence = "tree"), seed = 1)$privacy
  c2 <- synth(df, structure = ~ id / visit,
              privacy = dp_traj(1L, 2L, dependence = "tree"), seed = 1)$privacy
  # Adding conditioning columns must not move the noise scale at all.
  expect_equal(c1$noise, c0$noise, tolerance = 1e-12)
  expect_equal(c2$noise, c0$noise, tolerance = 1e-12)
  expect_equal(c1$longitudinal$n_transitions, c0$longitudinal$n_transitions)
  # Each variable records up to `cross` selected parents (never itself).
  cp <- c1$longitudinal$cross_parents
  expect_type(cp, "list")
  expect_true(all(vapply(names(cp),
    function(v) length(cp[[v]]) <= 1L && !(v %in% cp[[v]]), logical(1))))
})

# ---- backward compatibility ---------------------------------------------------

test_that("default (order 1, cross 0) is unchanged and prints no transition line", {
  df <- mk_traj()
  a <- synth(df, structure = ~ id / visit, privacy = dp_traj(1L, 0L), seed = 3)
  # first-order own-lag composition, tree n_init_marg = nV + C(nV,2).
  cap <- 5; nV <- 3
  n_init <- nV + nV * (nV - 1) / 2
  expect_equal(a$privacy$noise, (1 + n_init + nV * (cap - 1)) / 4, tolerance = 1e-9)
  expect_equal(a$privacy$longitudinal$order, 1L)
  expect_equal(a$privacy$longitudinal$cross, 0L)
  txt <- capture.output(print(a$privacy))
  expect_false(any(grepl("transitions: order", txt)))
})

test_that("default release is reproducible; order/cross release is reproducible", {
  df <- mk_traj()
  d1 <- as.data.frame(synth(df, structure = ~ id / visit,
                            privacy = dp_traj(2L, 1L), seed = 5))
  d2 <- as.data.frame(synth(df, structure = ~ id / visit,
                            privacy = dp_traj(2L, 1L), seed = 5))
  expect_identical(d1, d2)
  expect_identical(names(d1), names(df))
  expect_true(all(d1$hr >= 40 & d1$hr <= 140))
  expect_true(all(d1$sbp >= 60 & d1$sbp <= 240))
})

# ---- early-row marginalisation (deterministic, no-noise tensor) ---------------

test_that("depth tables are exact and depth-1 marginalises out the deeper own-lag", {
  # Two units of length 3, one variable x with 2 bins; identity noise so counts
  # are exact. cur rows (pos > order=2) are rows 3 and 6.
  codes <- list(x = c(1L, 2L, 1L, 2L, 1L, 2L))
  pos   <- c(1L, 2L, 3L, 1L, 2L, 3L)
  nbins <- c(x = 2L)
  te <- dp_fit_transition_tensors(codes, nbins, pos, order = 2L,
                                  cross_parents = list(x = character(0)),
                                  tv_vars = "x", add_noise = function(z) z)[["x"]]
  d2 <- te$depth_tables[[2]]     # rows = combo(own1, own2) col-major, cols = current
  d1 <- te$depth_tables[[1]]     # rows = own1, cols = current
  # tuple (cur=1, own1=2, own2=1) at row3; tuple (cur=2, own1=1, own2=2) at row6.
  # combo(own1,own2), own1 fastest: (1,1)=1 (2,1)=2 (1,2)=3 (2,2)=4.
  expect_equal(d2[2, ], c(1, 0))          # own1=2,own2=1 -> cur 1
  expect_equal(d2[3, ], c(0, 1))          # own1=1,own2=2 -> cur 2
  expect_equal(d2[1, ], c(0.5, 0.5))      # unseen -> uniform
  expect_equal(d2[4, ], c(0.5, 0.5))
  # depth-1: sum over own2. own1=1 saw cur2 once; own1=2 saw cur1 once.
  expect_equal(d1[1, ], c(0, 1))          # own1=1 -> cur 2
  expect_equal(d1[2, ], c(1, 0))          # own1=2 -> cur 1
})

test_that("order 1, cross 0 tensor path reproduces the simple transition matrix", {
  # depth-1 table with order=1 equals table(prev, cur) row-normalised.
  codes <- list(x = c(1L, 1L, 2L, 2L, 1L, 2L))
  pos   <- c(1L, 2L, 3L, 1L, 2L, 3L)
  nbins <- c(x = 2L)
  te <- dp_fit_transition_tensors(codes, nbins, pos, order = 1L,
                                  cross_parents = list(x = character(0)),
                                  tv_vars = "x", add_noise = function(z) z)[["x"]]
  # cur rows (pos>1): 2,3,5,6 -> pairs (prev,cur): (1,1)(1,2)(1,1)(2,2).
  simple <- dp_fit_transitions(list(x = codes$x[c(1L, 2L, 4L, 5L)]),
                               list(x = codes$x[c(2L, 3L, 5L, 6L)]),
                               nbins["x"], function(z) z)[["x"]]
  expect_equal(te$depth_tables[[1]], simple)
})

# ---- cross-parent selection ---------------------------------------------------

test_that("cross-parents are the highest mutual-information companions", {
  W <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  W["a", "b"] <- W["b", "a"] <- 0.9
  W["a", "c"] <- W["c", "a"] <- 0.1
  W["b", "c"] <- W["c", "b"] <- 0.4
  sel <- dp_select_cross_parents(W, tv_vars = c("a", "b", "c"),
                                 all_vars = c("a", "b", "c"), cross = 1L)
  expect_identical(sel$a, "b")            # a's strongest companion is b
  expect_identical(sel$b, "a")
  expect_identical(sel$c, "b")            # c prefers b (0.4) over a (0.1)
  # cross = 0 or NULL W -> no parents.
  expect_identical(dp_select_cross_parents(W, c("a"), c("a", "b"), 0L)$a,
                   character(0))
  expect_identical(dp_select_cross_parents(NULL, c("a"), c("a", "b"), 1L)$a,
                   character(0))
})

# ---- degenerate / interaction cases ------------------------------------------

test_that("all-baseline plus high order is a valid constant release", {
  df <- mk_traj()
  res <- synth(df, structure = ~ id / visit,
               privacy = dp_control(epsilon = 4, dependence = "tree",
                 transition_order = 2L, baseline = c("sex", "hr", "sbp"),
                 max_rows_per_person = 5L, bounds = traj_bounds,
                 domain = "public"), seed = 2)
  expect_equal(res$privacy$longitudinal$n_transitions, 0L)
  out <- as.data.frame(res)
  by_unit <- split(out, out$id)
  all_const <- vapply(by_unit, function(u)
    all(vapply(c("sex", "hr", "sbp"),
               function(v) length(unique(u[[v]])) == 1L, logical(1))),
    logical(1))
  expect_true(all(all_const))
})

test_that("baseline coexists with higher-order transitions on the rest", {
  df <- mk_traj()
  res <- synth(df, structure = ~ id / visit,
               privacy = dp_control(epsilon = 5, dependence = "tree",
                 transition_order = 2L, transition_cross = 1L,
                 baseline = "sex", max_rows_per_person = 5L,
                 bounds = traj_bounds, domain = "public"), seed = 4)
  # sex held constant; hr, sbp are the 2 time-varying transition histograms.
  expect_equal(res$privacy$longitudinal$n_transitions, 2L)
  expect_identical(res$privacy$longitudinal$baseline, "sex")
  out <- as.data.frame(res)
  by_unit <- split(out, out$id)
  const_sex <- vapply(by_unit, function(u) length(unique(u$sex)) == 1L, logical(1))
  expect_true(all(const_sex))
})

# ---- validation & prints ------------------------------------------------------

test_that("dp_control validates transition_order and transition_cross", {
  expect_silent(dp_control(epsilon = 1, transition_order = 2L))
  expect_silent(dp_control(epsilon = 1, transition_cross = 2L))
  expect_error(dp_control(epsilon = 1, transition_order = 0L), "transition_order")
  expect_error(dp_control(epsilon = 1, transition_order = 1.5), "transition_order")
  expect_error(dp_control(epsilon = 1, transition_cross = -1L), "transition_cross")
  expect_error(dp_control(epsilon = 1, transition_cross = 1.5), "transition_cross")
  expect_error(
    dp_control(epsilon = 1, dependence = "independent", transition_cross = 1L),
    "dependence")
})

test_that("order beyond the row cap is refused with a clear error", {
  df <- mk_traj()
  expect_error(
    synth(df, structure = ~ id / visit, privacy = dp_traj(order = 5L, cap = 5L),
          seed = 1),
    "transition_order")
})

test_that("order/cross are announced in the control and accounting prints", {
  ctrl <- capture.output(print(
    dp_control(epsilon = 1, transition_order = 2L, transition_cross = 1L)))
  expect_true(any(grepl("transitions: order 2", ctrl)))

  df <- mk_traj()
  acct <- synth(df, structure = ~ id / visit,
                privacy = dp_traj(2L, 1L), seed = 1)$privacy
  txt <- capture.output(print(acct))
  expect_true(any(grepl("transitions: order 2", txt)))
  expect_true(any(grepl("cross-parents", txt)))
})
