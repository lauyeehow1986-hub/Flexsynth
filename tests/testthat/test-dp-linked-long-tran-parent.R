# Track B: parent-conditioned transitions on a longitudinally-modelled linked
# child - dp_control(transition_parent = p). Today cross_table = TRUE
# cross-conditions only a longitudinal child's INITIAL state; the parent's
# influence then rides the own-lag chain and decays. transition_parent re-injects
# the (static) immediate-parent attributes into the transition tensor at every
# step, so parent -> child dependence stays anchored across the whole trajectory.
# It reuses the parent-by-child joints already measured (and paid for) by the
# cross-conditioned initial state, so it is BUDGET-NEUTRAL and REQUIRES
# cross_table = TRUE for that child.
#
# patients (root, `sex`) -> visits (longitudinal child):
#   * status - a 2-state chain whose stickiness AND marginal both depend on the
#     parent's sex (males start "worse" and are sticky; females start "stable"
#     and flip almost every visit), so sex predicts both the value (MI > 0, so it
#     is selected) and the DYNAMICS (the thing the transitions must carry);
#   * grade  - a sex-independent time-varying factor (a decoy companion).
lk_tp <- function(np = 500, seed = 1, stick_m = 0.9, stick_f = 0.1,
                  worse_m = 0.7, worse_f = 0.3) {
  set.seed(seed)
  patients <- data.frame(id = seq_len(np),
                         sex = factor(sample(c("F", "M"), np, TRUE)),
                         stringsAsFactors = FALSE)
  visits <- do.call(rbind, lapply(seq_len(np), function(i) {
    pid <- patients$id[i]; male <- patients$sex[i] == "M"
    k <- 3L + stats::rpois(1, 1.5)
    stick <- if (male) stick_m else stick_f
    p0    <- if (male) worse_m else worse_f
    st <- character(k)
    st[1L] <- if (stats::runif(1) < p0) "worse" else "stable"
    for (j in 2:k)
      st[j] <- if (stats::runif(1) < stick) st[j - 1L]
               else setdiff(c("stable", "worse"), st[j - 1L])
    data.frame(id = pid, visit_num = seq_len(k),
               status = factor(st, levels = c("stable", "worse")),
               grade  = factor(sample(c("A", "B", "C"), k, TRUE),
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

# Within-unit lag-1 agreement of `status`, split by the (synthetic) parent sex.
agree_by_sex <- function(res) {
  syn <- as.list(res)
  sex <- as.character(syn$patients$sex)[match(syn$visits$id, syn$patients$id)]
  ag <- function(status, id) {
    s <- split(as.character(status), id)
    num <- 0L; den <- 0L
    for (x in s) if (length(x) >= 2L) {
      num <- num + sum(x[-1L] == x[-length(x)]); den <- den + length(x) - 1L
    }
    if (den == 0L) NA_real_ else num / den
  }
  v <- syn$visits
  c(M = ag(v$status[sex == "M"], v$id[sex == "M"]),
    F = ag(v$status[sex == "F"], v$id[sex == "F"]))
}

mk_tp <- function(p, eps = 10, cap = 6L, cross = TRUE, extra = list()) {
  args <- c(list(epsilon = eps, mechanism = "laplace", dependence = "tree",
                 max_rows_per_person = c(visits = cap), longitudinal = "visits",
                 cross_table = cross, transition_parent = p, domain = "public"),
            extra)
  do.call(dp_control, args)
}

test_that("transition_parent is budget-neutral (same noise and histogram count)", {
  d <- lk_tp(np = 200)
  r0 <- synth_linked(d$tables, d$structures, d$keys, privacy = mk_tp(0L), seed = 1)
  r1 <- synth_linked(d$tables, d$structures, d$keys, privacy = mk_tp(1L), seed = 1)
  expect_equal(r1$privacy$noise, r0$privacy$noise)          # identical calibration
  expect_equal(r1$privacy$n_marginals, r0$privacy$n_marginals)
  expect_equal(by_tbl(r1)$visits$tran_parent, 1L)
  expect_equal(by_tbl(r0)$visits$tran_parent, 0L)
})

test_that("the exact composition matches the combined cross + longitudinal release", {
  d <- lk_tp(np = 200)
  cap <- 6L; eps <- 10
  res <- synth_linked(d$tables, d$structures, d$keys,
                      privacy = mk_tp(1L, eps = eps, cap = cap), seed = 1)
  # patients: sex, 1 one-way at cap 1                                          = 1
  # visits (cross + longi, nC = 2, nP = 1): init (cross) =
  #   dp_child_nvarmarg(2, 1, tree) = (2 + 1) + 2*1 = 5, at pcp 1              = 5
  #   2 transitions at pcp*(cap - 1) = 5 each                                  = 10
  #   count at cs = 1                                                          = 1
  # transition_parent adds NO histogram and NO sensitivity.
  total_l1 <- 1 + 5 + 10 + 1
  expect_equal(res$privacy$noise, total_l1 / eps)
  expect_equal(res$privacy$n_marginals, 1L + 5L + 2L + 1L)
})

test_that("transition_parent > 0 requires cross_table for that child", {
  d <- lk_tp(np = 80)
  expect_error(
    synth_linked(d$tables, d$structures, d$keys,
                 privacy = mk_tp(1L, cross = FALSE), seed = 1),
    "transition_parent")
})

test_that("parent dependence is carried into the transitions, not just the start", {
  d <- lk_tp(np = 600, seed = 3)
  # ground truth: males sticky, females flip -> a large per-sex agreement gap.
  r0 <- synth_linked(d$tables, d$structures, d$keys,
                     privacy = mk_tp(0L, eps = 12), seed = 11)   # init-cross only
  r1 <- synth_linked(d$tables, d$structures, d$keys,
                     privacy = mk_tp(1L, eps = 12), seed = 11)   # + parent-cond. steps
  g0 <- agree_by_sex(r0)
  g1 <- agree_by_sex(r1)
  gap0 <- abs(g0[["M"]] - g0[["F"]])
  gap1 <- abs(g1[["M"]] - g1[["F"]])
  # Conditioning the transitions on the parent restores the sex-specific
  # stickiness that pooled (init-cross-only) transitions wash out.
  expect_gt(gap1, gap0 + 0.2)
  expect_gt(gap1, 0.4)
})

test_that("a parent-conditioned longitudinal child reproduces and keeps integrity", {
  d <- lk_tp(np = 300, seed = 4)
  dp <- mk_tp(1L, eps = 10, cap = 6L)
  a <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 7)
  b <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 7)
  expect_equal(as.list(a)$visits, as.list(b)$visits)         # deterministic
  cl <- check_linkage(a)
  expect_true(all(cl$orphan_rows[!is.na(cl$orphan_rows)] == 0))
  # within-unit positions are contiguous 1..k
  syn <- as.list(a)
  pos <- split(syn$visits$visit_num, syn$visits$id)
  expect_true(all(vapply(pos, function(v) identical(sort(v), seq_along(v)),
                         logical(1))))
})

test_that("transition_parent composes with baseline, order and cross", {
  d <- lk_tp(np = 250, seed = 5)
  # add a baseline column to the child so all four controls are exercised at once.
  d$tables$visits$egfr_band <- factor(
    ave(as.character(d$tables$visits$id), d$tables$visits$id,
        FUN = function(x) rep(sample(c("G1", "G2", "G3"), 1L), length(x))),
    levels = c("G1", "G2", "G3"))
  dp <- dp_control(epsilon = 10, mechanism = "laplace", dependence = "tree",
                   max_rows_per_person = c(visits = 6L), longitudinal = "visits",
                   cross_table = TRUE, baseline = "egfr_band",
                   transition_order = 2L, transition_cross = 1L,
                   transition_parent = 1L, domain = "public")
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  bt <- by_tbl(res)$visits
  expect_true(bt$cross_init)
  expect_equal(bt$baseline, "egfr_band")
  expect_equal(bt$tran_order, 2L)
  expect_equal(bt$tran_cross, 1L)
  expect_equal(bt$tran_parent, 1L)
  # baseline still byte-constant within a unit
  syn <- as.list(res)
  nuniq <- tapply(as.character(syn$visits$egfr_band), syn$visits$id,
                  function(x) length(unique(x)))
  expect_true(all(nuniq == 1L))
})

test_that("the accounting print flags parent-conditioned transitions", {
  d <- lk_tp(np = 80, seed = 5)
  dp <- mk_tp(1L, eps = 8, cap = 5L)
  res <- synth_linked(d$tables, d$structures, d$keys, privacy = dp, seed = 1)
  out <- paste(utils::capture.output(print(res$privacy)), collapse = "\n")
  expect_match(out, "parent-attr")
  expect_match(out, "parent-attrs: (status|grade) ~ sex")
})
