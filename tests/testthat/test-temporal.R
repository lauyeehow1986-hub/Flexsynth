# Phase 2: learned structural (count) model + within-unit temporal dependence.

# A panel with strong within-subject autocorrelation. A latent per-subject level
# drives both visits, but is NOT a column, so it can only be recovered through
# the lag-1 predictor. n is large to keep the correlation check stable.
make_autocorr <- function(n_subj = 300, seed = 101) {
  set.seed(seed)
  level <- rnorm(n_subj, 130, 16)          # latent subject level (not a column)
  data.frame(
    id    = rep(seq_len(n_subj), each = 2),
    visit = rep(1:2, times = n_subj),
    sbp   = round(c(rbind(level + rnorm(n_subj, 0, 3),
                          level + rnorm(n_subj, 0, 3)))),
    stringsAsFactors = FALSE
  )
}

wide_by_visit <- function(df) {
  v1 <- df$sbp[df$visit == 1][order(df$id[df$visit == 1])]
  v2 <- df$sbp[df$visit == 2][order(df$id[df$visit == 2])]
  data.frame(v1 = v1, v2 = v2)
}

test_that("within-unit autocorrelation is carried into the synthetic data", {
  df <- make_autocorr()
  real_cor <- cor(wide_by_visit(df)$v1, wide_by_visit(df)$v2)
  expect_gt(real_cor, 0.7)                  # the signal is genuinely there

  syn <- as.data.frame(synth(df, structure = ~ id / visit, seed = 4))
  syn_cor <- cor(wide_by_visit(syn)$v1, wide_by_visit(syn)$v2)
  expect_gt(syn_cor, 0.4)                   # the lag model recovers most of it
})

test_that("the learned count model reproduces the rows-per-unit distribution", {
  set.seed(202)
  n_subj <- 250
  sizes <- sample(1:4, n_subj, replace = TRUE, prob = c(0.2, 0.4, 0.25, 0.15))
  df <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(id = i, visit = seq_len(sizes[i]),
               sbp = round(rnorm(sizes[i], 130, 15)))
  }))

  res <- synth(df, structure = ~ id / visit, seed = 6)
  syn <- as.data.frame(res)
  syn_sizes <- as.integer(table(syn$id))

  expect_true(all(syn_sizes %in% sizes))                 # only real sizes appear
  expect_lt(abs(mean(syn_sizes) - mean(sizes)), 0.5)     # mean count preserved
  # within a synthetic unit the visit index is ordered 1..size
  ok <- tapply(syn$visit, syn$id, function(v) identical(v, seq_along(v)))
  expect_true(all(ok))
})

test_that("variable-length units synthesise without error and keep schema", {
  set.seed(303)
  sizes <- sample(2:3, 60, replace = TRUE)
  df <- do.call(rbind, lapply(seq_along(sizes), function(i) {
    data.frame(id = i, visit = seq_len(sizes[i]),
               age = 60L, sbp = round(rnorm(sizes[i], 130, 15)),
               drug = sample(c("A", "B"), sizes[i], replace = TRUE),
               stringsAsFactors = FALSE)
  }))
  syn <- as.data.frame(synth(df, structure = ~ id / visit, seed = 8))
  expect_named(syn, names(df))
  expect_identical(vapply(syn, class, character(1)),
                   vapply(df, class, character(1)))
  expect_true(all(syn$drug %in% c("A", "B")))
})
