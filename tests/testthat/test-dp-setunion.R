# DP set-union for character category discovery under domain = "dp". A bare
# character column no longer has to be pre-converted to a factor: its present
# category set is discovered privately with a stability histogram (add Laplace
# noise to each present category's count, keep those over a threshold that hides
# any category a single person could have created). Rare/unique categories are
# dropped into an appended "(other)" catch-all. Needs delta > 0 (hiding a
# category's very presence needs the threshold's delta slack). Budget reuses the
# domain_frac slice; composition stays exact.

setunion_data <- function(n = 3000L, seed = 1L) {
  set.seed(seed)
  grp <- sample(c("aa", "bb", "cc"), n, TRUE, prob = c(0.5, 0.3, 0.2))
  n_rare <- round(n * 0.02)
  grp[seq_len(n_rare)] <- sprintf("rare%05d", seq_len(n_rare))   # ~unique tail
  data.frame(id = seq_len(n), grp = grp,
             age = round(stats::rnorm(n, 60, 8)),
             stringsAsFactors = FALSE)
}

# ---- the stability-histogram core ------------------------------------------

test_that("dp_discover_categories keeps frequent, drops rare categories", {
  x <- c(rep("aa", 400), rep("bb", 300), rep("cc", 200), sprintf("u%d", 1:80))
  set.seed(1)
  kept <- flexsynth:::dp_discover_categories(x, eps_op = 1, delta_cat = 1e-6,
                                             C = 1L, cap = 1L)
  expect_true(all(c("aa", "bb", "cc") %in% kept))
  expect_false(any(sprintf("u%d", 1:80) %in% kept))       # singletons dropped
})

test_that("dp_discover_categories keeps (weakly) more as eps_op grows", {
  x <- c(rep("aa", 400), rep("bb", 300), rep("cc", 200), sprintf("u%d", 1:80))
  set.seed(2); lo <- flexsynth:::dp_discover_categories(x, 0.05, 1e-6, 1L, 1L)
  set.seed(2); hi <- flexsynth:::dp_discover_categories(x, 5, 1e-6, 1L, 1L)
  expect_gte(length(hi), length(lo))
  expect_true(all(hi %in% c("aa", "bb", "cc")))           # never invents a level
})

test_that("dp_discover_categories tolerates NA and empty input", {
  expect_identical(flexsynth:::dp_discover_categories(character(0), 1, 1e-6, 1L, 1L),
                   character(0))
  x <- c(rep("aa", 300), NA, NA)
  set.seed(1)
  expect_true("aa" %in% flexsynth:::dp_discover_categories(x, 1, 1e-6, 1L, 1L))
})

# ---- end-to-end release ----------------------------------------------------

test_that("a character column is discovered under domain = 'dp' + gaussian", {
  df <- setunion_data(3000L, 3L)
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   bounds = list(age = c(20, 100)))
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  b <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 4))
  expect_named(a, names(df))
  expect_false(anyNA(a))
  expect_identical(a, b)                                   # reproducible
  # Synthetic categories are a subset of {discovered frequent} u {"(other)"}.
  expect_true(all(unique(a$grp) %in% c("aa", "bb", "cc", "(other)")))
  # The frequent categories survived; the rare tail collapsed to "(other)".
  expect_true(all(c("aa", "bb", "cc") %in% unique(a$grp)))
  expect_false(any(grepl("^rare", a$grp)))
})

test_that("discovery is recorded in the accounting and reported budget is exact", {
  df <- setunion_data(2500L, 5L)
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   bounds = list(age = c(20, 100)))
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 6)$privacy
  expect_equal(ac$epsilon, 6)
  expect_equal(ac$delta, 1e-6)
  expect_false(is.null(ac$domain$categorical))
  expect_true("grp" %in% ac$domain$categorical$vars)
})

# ---- delta > 0 requirement -------------------------------------------------

test_that("pure-eps laplace refuses a character column with a pointed message", {
  df <- setunion_data(1500L, 7L)
  dp <- dp_control(epsilon = 4, mechanism = "laplace",   # delta = 0
                   bounds = list(age = c(20, 100)))
  expect_error(synth(df, structure = ~ id, privacy = dp, seed = 1),
               "delta")
})

test_that("delta > 0 laplace discovers a character column", {
  df <- setunion_data(2000L, 8L)
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "laplace",
                   bounds = list(age = c(20, 100)))
  a <- as.data.frame(synth(df, structure = ~ id, privacy = dp, seed = 2))
  expect_false(anyNA(a))
  expect_true(all(unique(a$grp) %in% c("aa", "bb", "cc", "(other)")))
})

# ---- no effect on non-character releases -----------------------------------

test_that("a numeric-only release is unaffected (no categorical block engaged)", {
  df <- data.frame(id = 1:2000, age = round(stats::rnorm(2000, 60, 8)),
                   sbp = round(stats::rnorm(2000, 130, 15)))
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   bounds = list(age = c(20, 100), sbp = c(60, 240)))
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 3)$privacy
  expect_null(ac$domain$categorical)
  expect_equal(ac$epsilon, 6)
})

# ---- print -----------------------------------------------------------------

test_that("accounting print surfaces the discovered categories", {
  df <- setunion_data(1500L, 9L)
  dp <- dp_control(epsilon = 6, delta = 1e-6, mechanism = "gaussian",
                   bounds = list(age = c(20, 100)))
  ac <- synth(df, structure = ~ id, privacy = dp, seed = 1)$privacy
  expect_output(print(ac), "categor")
})
