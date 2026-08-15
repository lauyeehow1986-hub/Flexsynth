# Mixed-type (Gower-style) distance helpers used by the disclosure-risk
# diagnostics. Not exported.

# Learn per-column scaling from a reference frame: numeric columns are scaled by
# their observed range, categorical columns compared by (in)equality.
distance_scales <- function(ref, vars) {
  stats::setNames(lapply(vars, function(v) {
    x <- ref[[v]]
    if (is.numeric(x)) {
      rng <- suppressWarnings(diff(range(x, na.rm = TRUE)))
      list(type = "numeric", range = if (is.finite(rng) && rng > 0) rng else 1)
    } else {
      list(type = "categorical")
    }
  }), vars)
}

# For each row of `a`, the minimum Gower distance to any row of `b`.
# `exclude_self = TRUE` (used only when `a` and `b` are the same frame) ignores
# the trivial zero-distance match of a row to itself.
nn_distance <- function(a, b, vars, scales, exclude_self = FALSE) {
  nb <- nrow(b)
  p  <- length(vars)
  out <- numeric(nrow(a))
  bcols <- lapply(vars, function(v) b[[v]])
  acols <- lapply(vars, function(v) a[[v]])
  for (i in seq_len(nrow(a))) {
    di <- numeric(nb)
    for (j in seq_len(p)) {
      s  <- scales[[j]]
      av <- acols[[j]][i]
      bv <- bcols[[j]]
      if (identical(s$type, "numeric")) {
        dj <- abs(bv - av) / s$range
      } else {
        dj <- as.numeric(as.character(bv) != as.character(av))
      }
      dj[is.na(dj)] <- 1                 # NA / mismatch counts as max distance
      di <- di + dj
    }
    di <- di / p
    if (exclude_self) di[i] <- Inf
    out[i] <- min(di)
  }
  out
}

# Rank-based (Mann-Whitney) AUC of `score` separating positives from negatives.
# Larger score should indicate a positive. Returns 0.5 when uninformative.
auc_mw <- function(score, positive) {
  positive <- as.logical(positive)
  np <- sum(positive)
  nn <- sum(!positive)
  if (np == 0L || nn == 0L) return(NA_real_)
  r <- rank(score)
  (sum(r[positive]) - np * (np + 1) / 2) / (np * nn)
}
