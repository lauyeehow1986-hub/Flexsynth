# Full AIM (select = "aim"): loopy adaptive marginal selection over a bounded-
# treewidth graphical model, reconciled and sampled with Private-PGM. Covers the
# triangulation primitive (measured marginals -> maximal cliques + junction tree),
# the general PGM sampler (set-valued separators / new-variable blocks), the
# end-to-end AIM fit that captures a loop of pairwise dependence a tree cannot,
# the control validation and gates, budget-neutral accounting, reproducibility,
# and the prints.

# ---- triangulation ---------------------------------------------------------

# For a set of clique variable sets and their junction-tree edges, the running
# intersection property: for every variable, the cliques containing it induce a
# connected subtree of the junction tree.
jt_rip_holds <- function(cliques, edges, d) {
  K <- length(cliques)
  if (K <= 1L) return(TRUE)
  adj <- vector("list", K)
  for (e in edges) { adj[[e[1L]]] <- c(adj[[e[1L]]], e[2L])
                     adj[[e[2L]]] <- c(adj[[e[2L]]], e[1L]) }
  for (v in seq_len(d)) {
    has <- which(vapply(cliques, function(cl) v %in% cl, logical(1)))
    if (length(has) <= 1L) next
    # BFS over the sub-adjacency restricted to cliques containing v.
    seen <- has[1L]; queue <- has[1L]
    while (length(queue)) {
      a <- queue[1L]; queue <- queue[-1L]
      for (b in adj[[a]]) if (b %in% has && !(b %in% seen)) {
        seen <- c(seen, b); queue <- c(queue, b) }
    }
    if (!setequal(seen, has)) return(FALSE)
  }
  TRUE
}

test_that("dp_triangulate collapses a 3-cycle into one triangle clique", {
  tr <- flexsynth:::dp_triangulate(3L, list(c(1L, 2L), c(2L, 3L), c(1L, 3L)))
  expect_length(tr$cliques, 1L)
  expect_setequal(tr$cliques[[1L]], 1:3)
  expect_equal(tr$width, 3L)
})

test_that("dp_triangulate keeps a path as two pairwise cliques with a valid JT", {
  tr <- flexsynth:::dp_triangulate(3L, list(c(1L, 2L), c(2L, 3L)))
  expect_length(tr$cliques, 2L)
  expect_equal(tr$width, 2L)
  # both input edges are covered by a clique
  covered <- function(e) any(vapply(tr$cliques, function(cl) all(e %in% cl), logical(1)))
  expect_true(covered(c(1L, 2L))); expect_true(covered(c(2L, 3L)))
  expect_true(jt_rip_holds(tr$cliques, tr$edges, 3L))
})

test_that("dp_triangulate emits a singleton clique for an isolated variable", {
  tr <- flexsynth:::dp_triangulate(3L, list(c(1L, 2L)))
  # variable 3 is in no edge: it must still appear in some clique (singleton).
  inclq <- function(v) any(vapply(tr$cliques, function(cl) v %in% cl, logical(1)))
  expect_true(inclq(1L)); expect_true(inclq(2L)); expect_true(inclq(3L))
  expect_equal(tr$width, 2L)
  expect_true(any(vapply(tr$cliques, function(cl) identical(sort(cl), 3L), logical(1))))
  expect_true(jt_rip_holds(tr$cliques, tr$edges, 3L))
})

test_that("dp_triangulate triangulates a 4-cycle into width-3 cliques", {
  tr <- flexsynth:::dp_triangulate(
    4L, list(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(1L, 4L)))
  expect_equal(tr$width, 3L)           # a chord is added -> two triangles
  expect_true(jt_rip_holds(tr$cliques, tr$edges, 4L))
  # every var present, every input edge covered
  for (v in 1:4)
    expect_true(any(vapply(tr$cliques, function(cl) v %in% cl, logical(1))))
})

# ---- general PGM sampler ---------------------------------------------------

# A factor over a sorted variable set, values column-major (first var fastest).
fac <- function(vars, dims, vals) list(vars = as.integer(vars),
                                       dims = as.integer(dims), vals = vals)

test_that("dp_sample_codes_pgm reproduces clique marginals over a chain", {
  # Chain A-B-C: clique {1,2} with joint P(A,B), clique {2,3} with joint P(B,C),
  # sharing the same P(B). Built from P(A), P(B|A), P(C|B).
  pa <- c(0.6, 0.4); pbA <- rbind(c(0.8, 0.2), c(0.2, 0.8))
  pab <- matrix(0, 2, 2); for (a in 1:2) for (b in 1:2) pab[a, b] <- pa[a] * pbA[a, b]
  pb <- colSums(pab); pcB <- rbind(c(0.7, 0.3), c(0.3, 0.7))
  pbc <- matrix(0, 2, 2); for (b in 1:2) for (cc in 1:2) pbc[b, cc] <- pb[b] * pcB[b, cc]
  model <- list(kind = "pgm", vars = c("A", "B", "C"), nbins = c(2L, 2L, 2L),
                cliques = list(c(1L, 2L), c(2L, 3L)), edges = list(c(1L, 2L)),
                beliefs = list(fac(c(1, 2), c(2, 2), as.vector(pab)),
                               fac(c(2, 3), c(2, 2), as.vector(pbc))))
  set.seed(1)
  codes <- flexsynth:::dp_sample_codes_pgm(model, 40000L)
  eab <- prop.table(table(factor(codes[, 1], 1:2), factor(codes[, 2], 1:2)))
  ebc <- prop.table(table(factor(codes[, 2], 1:2), factor(codes[, 3], 1:2)))
  expect_equal(as.numeric(eab), as.vector(pab), tolerance = 0.02)
  expect_equal(as.numeric(ebc), as.vector(pbc), tolerance = 0.02)
})

test_that("dp_sample_codes_pgm samples a single-triangle three-way belief", {
  # A strong three-way structure (near-XOR): only even-parity cells carry mass.
  T <- array(0, dim = c(2, 2, 2))
  for (a in 1:2) for (b in 1:2) for (cc in 1:2)
    T[a, b, cc] <- if (((a + b + cc) %% 2L) == 1L) 0.24 else 0.01
  T <- T / sum(T)
  model <- list(kind = "pgm", vars = c("A", "B", "C"), nbins = c(2L, 2L, 2L),
                cliques = list(c(1L, 2L, 3L)), edges = list(),
                beliefs = list(fac(1:3, c(2, 2, 2), as.vector(T))))
  set.seed(2)
  codes <- flexsynth:::dp_sample_codes_pgm(model, 60000L)
  E <- table(factor(codes[, 1], 1:2), factor(codes[, 2], 1:2), factor(codes[, 3], 1:2))
  E <- E / sum(E)
  expect_equal(as.numeric(E), as.numeric(T), tolerance = 0.02)
})

test_that("dp_sample_codes_pgm handles an isolated variable via an empty separator", {
  pab <- matrix(c(0.3, 0.2, 0.1, 0.4), 2, 2)     # P(A,B)
  pc <- c(0.75, 0.25)                            # P(C), independent
  model <- list(kind = "pgm", vars = c("A", "B", "C"), nbins = c(2L, 2L, 2L),
                cliques = list(c(1L, 2L), 3L), edges = list(c(1L, 2L)),
                beliefs = list(fac(c(1, 2), c(2, 2), as.vector(pab)),
                               fac(3L, 2L, pc)))
  set.seed(3)
  codes <- flexsynth:::dp_sample_codes_pgm(model, 40000L)
  eab <- prop.table(table(factor(codes[, 1], 1:2), factor(codes[, 2], 1:2)))
  ec  <- prop.table(table(factor(codes[, 3], 1:2)))
  expect_equal(as.numeric(eab), as.vector(pab), tolerance = 0.02)
  expect_equal(as.numeric(ec), pc, tolerance = 0.02)
})

# ---- end-to-end AIM fit ----------------------------------------------------

# Three binary variables with a *loop* of pairwise dependence: every pair is
# correlated, but no single tree edge-set can hold all three correlations. AIM at
# treewidth 2 measures all three pairs and reconciles them into one triangle.
loop_data <- function(n = 4000L, seed = 11L) {
  set.seed(seed)
  a <- sample(0:1, n, TRUE)
  b <- ifelse(runif(n) < 0.85, a, 1L - a)              # B tracks A
  cc <- ifelse(runif(n) < 0.85, b, 1L - b)             # C tracks B
  a2 <- ifelse(runif(n) < 0.7, cc, a)                  # and A tracks C (closes loop)
  data.frame(id = seq_len(n),
             A = factor(a2), B = factor(b), C = factor(cc))
}

test_that("AIM at treewidth 2 reconciles a three-variable loop into one triangle", {
  df <- loop_data()
  codes <- lapply(c("A", "B", "C"), function(v) as.integer(df[[v]]))
  names(codes) <- c("A", "B", "C")
  nbins <- c(A = 2L, B = 2L, C = 2L)
  dp <- dp_control(epsilon = 40, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2)
  calib <- flexsynth:::dp_calibrate(dp, 6L, 1, budget_frac = 0.75)
  set.seed(5)
  model <- flexsynth:::dp_fit_model_aim(codes, nbins, dp, calib, w = 2L,
                                        sel_eps = 5, cap = 1)
  expect_identical(model$kind, "pgm")
  # all three pairs -> complete graph on 3 vars -> a single triangle clique.
  expect_length(model$cliques, 1L)
  expect_setequal(model$cliques[[1L]], 1:3)
  # the reconciled model reproduces every pairwise marginal (a tree cannot).
  bel <- model$beliefs[[1L]]
  A <- array(bel$vals, dim = bel$dims)
  pab <- apply(A, c(1, 2), sum); pac <- apply(A, c(1, 3), sum)
  tab <- function(i, j) prop.table(table(factor(codes[[i]], 1:2), factor(codes[[j]], 1:2)))
  expect_equal(as.numeric(pab), as.numeric(tab(1, 2)), tolerance = 0.03)
  expect_equal(as.numeric(pac), as.numeric(tab(1, 3)), tolerance = 0.03)
})

test_that("AIM at treewidth 1 stays a tree (no clique larger than a pair)", {
  df <- loop_data()
  codes <- lapply(c("A", "B", "C"), function(v) as.integer(df[[v]]))
  names(codes) <- c("A", "B", "C")
  nbins <- c(A = 2L, B = 2L, C = 2L)
  dp <- dp_control(epsilon = 40, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 1)
  calib <- flexsynth:::dp_calibrate(dp, 5L, 1, budget_frac = 0.75)
  set.seed(6)
  model <- flexsynth:::dp_fit_model_aim(codes, nbins, dp, calib, w = 1L,
                                        sel_eps = 5, cap = 1)
  expect_true(max(vapply(model$cliques, length, integer(1))) <= 2L)
})

# ---- control validation ----------------------------------------------------

test_that("dp_control accepts select = aim and stores it", {
  dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2)
  expect_identical(dp$select, "aim")
})

test_that("dp_control rejects aim combined with incompatible knobs", {
  base <- list(epsilon = 1, delta = 1e-6, mechanism = "gaussian", select = "aim")
  # annealed AIM is not offered yet.
  expect_error(do.call(dp_control, c(base, list(anneal = TRUE))),
               "aim|anneal|AIM")
  # structure_frac is the fixed tree's knob.
  expect_error(do.call(dp_control, c(base, list(structure_frac = 0.3))),
               "aim|structure|AIM")
  # degree > 1 is a Bayesian network, an alternative structure search.
  expect_error(do.call(dp_control, c(base, list(degree = 2))),
               "aim|degree|network|AIM")
})

# ---- end-to-end release + accounting ---------------------------------------

test_that("select = aim produces a well-formed, reproducible release", {
  df <- data.frame(id = 1:600,
                   x = factor(sample(0:1, 600, TRUE)),
                   y = factor(sample(0:2, 600, TRUE)),
                   z = factor(sample(0:1, 600, TRUE)))
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2)
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 7))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 7))
  expect_identical(sort(names(a)), sort(names(df)))
  expect_identical(a, b)                          # reproducible
  expect_true(all(as.character(a$y) %in% as.character(0:2)))
})

test_that("aim accounting reports d one-way + n_rounds selected marginals", {
  df <- data.frame(id = 1:600,
                   x = factor(sample(0:1, 600, TRUE)),
                   y = factor(sample(0:2, 600, TRUE)),
                   z = factor(sample(0:1, 600, TRUE)))
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2)
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 8)$privacy
  d <- 3L
  n_rounds <- min(d * (d - 1L) / 2L, 2L * (d - 1L))   # min(3, 4) = 3
  expect_equal(ac$n_marginals, d + n_rounds)
  expect_equal(ac$epsilon, 8)
  expect_false(is.null(ac$aim))
})

# ---- gates -----------------------------------------------------------------

test_that("select = aim is refused on longitudinal and linked DP releases", {
  dfl <- data.frame(id = rep(1:200, each = 2),
                    visit = rep(1:2, times = 200),
                    x = factor(sample(0:1, 400, TRUE)),
                    y = factor(sample(0:1, 400, TRUE)))
  dp <- dp_control(epsilon = 4, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2, max_rows_per_person = 2)
  expect_error(synth(dfl, structure = ~ id / visit, privacy = dp, seed = 1),
               "aim|AIM|flat")

  patients <- data.frame(id = 1:60, sex = factor(sample(c("F", "M"), 60, TRUE)))
  adm <- do.call(rbind, lapply(patients$id, function(pid)
    data.frame(id = pid, aid = seq_len(1 + rpois(1, 0.6)))))
  adm$los <- factor(sample(0:2, nrow(adm), TRUE))
  expect_error(
    synth_linked(tables = list(patients = patients, admissions = adm),
                 structures = list(patients = ~ id, admissions = ~ id / aid),
                 keys = list(patients = "id", admissions = c("id", "aid")),
                 privacy = dp, seed = 1),
    "aim|AIM|flat")
})

# ---- prints ----------------------------------------------------------------

test_that("prints surface the AIM selector", {
  dp <- dp_control(epsilon = 8, delta = 1e-6, mechanism = "gaussian",
                   select = "aim", treewidth = 2)
  expect_output(print(dp), "AIM|aim|loop")

  df <- data.frame(id = 1:400,
                   x = factor(sample(0:1, 400, TRUE)),
                   y = factor(sample(0:1, 400, TRUE)),
                   z = factor(sample(0:2, 400, TRUE)))
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 1)$privacy
  expect_output(print(ac), "AIM|aim|loop")
})
