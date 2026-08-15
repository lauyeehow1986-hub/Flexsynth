# Utility diagnostics: how close is the synthetic data to the real data?
# Track A carries no formal guarantee, so honest utility (and, in disclosure.R,
# risk) diagnostics are how quality is judged.

# --- input coercion -------------------------------------------------------

# Coerce a synth object or data.frame to a plain data.frame; leave lists alone.
as_diag_frame <- function(x) {
  if (inherits(x, "synth_result")) return(as.data.frame(x))
  if (inherits(x, "synth_linked_result")) return(as.list(x))
  x
}

# The columns present in both frames (order follows `real`).
common_columns <- function(real, syn) {
  intersect(names(real), names(syn))
}

# --- per-variable (marginal) comparison -----------------------------------

# Kolmogorov-Smirnov style max ECDF gap between two numeric samples.
ecdf_gap <- function(xr, xs) {
  xr <- xr[!is.na(xr)]; xs <- xs[!is.na(xs)]
  if (!length(xr) || !length(xs)) return(NA_real_)
  grid <- sort(unique(c(xr, xs)))
  max(abs(stats::ecdf(xr)(grid) - stats::ecdf(xs)(grid)))
}

# Total variation distance between two categorical distributions.
tv_distance <- function(xr, xs) {
  xr <- as.character(xr); xs <- as.character(xs)
  lev <- union(unique(xr), unique(xs))
  pr <- tabulate(match(xr, lev), length(lev)) / length(xr)
  ps <- tabulate(match(xs, lev), length(lev)) / length(xs)
  0.5 * sum(abs(pr - ps))
}

univariate_diagnostics <- function(real, syn, vars) {
  rows <- lapply(vars, function(v) {
    xr <- real[[v]]; xs <- syn[[v]]
    if (is.numeric(xr)) {
      data.frame(
        variable  = v,
        type      = "numeric",
        real_mean = mean(xr, na.rm = TRUE),
        syn_mean  = mean(xs, na.rm = TRUE),
        real_sd   = stats::sd(xr, na.rm = TRUE),
        syn_sd    = stats::sd(xs, na.rm = TRUE),
        distance  = ecdf_gap(xr, xs),          # KS statistic in [0, 1]
        metric    = "ks",
        stringsAsFactors = FALSE
      )
    } else {
      novel <- setdiff(unique(as.character(xs)), unique(as.character(xr)))
      data.frame(
        variable  = v,
        type      = if (is.logical(xr)) "logical" else "categorical",
        real_mean = NA_real_, syn_mean = NA_real_,
        real_sd   = NA_real_, syn_sd   = NA_real_,
        distance  = tv_distance(xr, xs),       # total variation in [0, 1]
        metric    = "tvd",
        stringsAsFactors = FALSE
      )
    }
  })
  do.call(rbind, rows)
}

# --- correlation structure ------------------------------------------------

correlation_diagnostics <- function(real, syn, vars) {
  num <- vars[vapply(vars, function(v) is.numeric(real[[v]]), logical(1))]
  if (length(num) < 2L) return(NULL)
  cr <- stats::cor(real[num], use = "pairwise.complete.obs")
  cs <- stats::cor(syn[num],  use = "pairwise.complete.obs")
  ut <- upper.tri(cr)
  diff <- (cr - cs)[ut]
  diff <- diff[is.finite(diff)]
  list(
    vars          = num,
    real          = cr,
    syn           = cs,
    frobenius     = sqrt(sum(diff^2)),
    mean_abs_diff = mean(abs(diff)),
    max_abs_diff  = if (length(diff)) max(abs(diff)) else NA_real_
  )
}

# --- propensity score utility (pMSE) --------------------------------------

# General utility: fit a model discriminating real from synthetic; if the two
# are indistinguishable the fitted propensities all sit at c = n_syn / N and the
# pMSE is near its null expectation (ratio ~ 1).
pmse_diagnostics <- function(real, syn, vars, propensity = "logistic") {
  keep <- vars[vapply(vars, function(v)
    length(unique(c(as.character(real[[v]]), as.character(syn[[v]])))) > 1L,
    logical(1))]
  if (!length(keep)) return(NULL)

  stacked <- rbind(real[keep], syn[keep])
  z <- c(rep(0L, nrow(real)), rep(1L, nrow(syn)))
  N <- length(z)
  cc <- nrow(syn) / N

  fit_ok <- FALSE
  if (identical(propensity, "cart") && requireNamespace("rpart", quietly = TRUE)) {
    m <- tryCatch(
      rpart::rpart(z ~ ., data = cbind(z = z, stacked), method = "class"),
      error = function(e) NULL)
    if (!is.null(m)) {
      ph <- stats::predict(m)[, 2L]
      npar <- NA_integer_; fit_ok <- TRUE
    }
  }
  if (!fit_ok) {                            # default: main-effects logistic
    m <- tryCatch(
      suppressWarnings(stats::glm(z ~ ., data = cbind(z = z, stacked),
                                  family = stats::binomial())),
      error = function(e) NULL)
    if (is.null(m)) return(NULL)
    ph   <- stats::fitted(m)
    npar <- length(stats::coef(m))
    propensity <- "logistic"
  }

  pmse <- mean((ph - cc)^2)
  # Snoke et al. (2018) null expectation for a logistic model with npar terms.
  expected <- if (!is.na(npar))
    (npar - 1) * (1 - cc)^2 * cc / N else NA_real_
  list(
    method     = propensity,
    pmse       = pmse,
    expected   = expected,
    ratio      = if (!is.na(expected) && expected > 0) pmse / expected else NA_real_,
    c          = cc,
    npar       = npar,
    n_vars     = length(keep)
  )
}

# --- public entry point ---------------------------------------------------

#' Utility diagnostics comparing synthetic data to the real data
#'
#' Quantifies how faithfully synthetic data reproduces the real data's
#' distributions and dependence structure. Because Track A output carries no
#' formal privacy guarantee, these utility diagnostics (and the risk
#' diagnostics in [disclosure_risk()]) are how synthesis quality is judged.
#'
#' Three views are reported:
#' \itemize{
#'   \item **Univariate** — per-variable marginal fit: a Kolmogorov-Smirnov
#'     statistic for numeric variables and a total-variation distance for
#'     categorical ones (both in \eqn{[0, 1]}; smaller is better).
#'   \item **Correlation** — the Frobenius and mean absolute difference between
#'     the real and synthetic numeric correlation matrices.
#'   \item **Propensity (pMSE)** — a general utility score: a model is fitted to
#'     tell real from synthetic records; when they are indistinguishable the
#'     propensity mean squared error sits near its null expectation, so the
#'     `ratio` is near 1 (larger means more distinguishable).
#' }
#'
#' If `real` and `syn` are named lists (or `syn` is a `synth_linked_result`),
#' each table is diagnosed and a per-table result is returned.
#'
#' @param real The real `data.frame` (or a named list of them for linked data).
#' @param syn The synthetic data: a `data.frame`, a [synth()] `synth_result`, a
#'   named list of tables, or a [synth_linked()] `synth_linked_result`.
#' @param vars Optional character vector restricting the variables compared;
#'   defaults to all columns present in both.
#' @param propensity Propensity model for the pMSE, `"logistic"` (default) or
#'   `"cart"` (needs the `rpart` package).
#' @param ... Unused.
#'
#' @return A `flexsynth_diagnostics` object (or a `flexsynth_diagnostics_list`
#'   for linked / list input). Has `print()` and `plot()` methods.
#' @seealso [disclosure_risk()] for privacy risk.
#' @export
#' @examples
#' df <- data.frame(
#'   id  = 1:200,
#'   age = round(rnorm(200, 60, 10)),
#'   sbp = round(rnorm(200, 130, 15))
#' )
#' res <- synth(df, ~ id, seed = 1)
#' diagnose(df, res)
diagnose <- function(real, syn, vars = NULL, propensity = "logistic", ...) {
  syn <- as_diag_frame(syn)

  if (is.list(real) && !is.data.frame(real) &&
      is.list(syn) && !is.data.frame(syn)) {
    tbls <- intersect(names(real), names(syn))
    if (!length(tbls)) stop("`real` and `syn` share no table names.", call. = FALSE)
    out <- stats::setNames(lapply(tbls, function(t)
      diagnose(real[[t]], syn[[t]], vars = vars, propensity = propensity)), tbls)
    return(structure(out, class = "flexsynth_diagnostics_list"))
  }

  if (!is.data.frame(real)) stop("`real` must be a data.frame.", call. = FALSE)
  if (!is.data.frame(syn))  stop("`syn` must be a data.frame or synth result.",
                                 call. = FALSE)

  cols <- common_columns(real, syn)
  if (!is.null(vars)) {
    miss <- setdiff(vars, cols)
    if (length(miss))
      stop(sprintf("`vars` not present in both frames: %s",
                   paste(miss, collapse = ", ")), call. = FALSE)
    cols <- intersect(vars, cols)
  }
  if (!length(cols)) stop("`real` and `syn` share no columns.", call. = FALSE)

  structure(
    list(
      univariate  = univariate_diagnostics(real, syn, cols),
      correlation = correlation_diagnostics(real, syn, cols),
      pmse        = pmse_diagnostics(real, syn, cols, propensity),
      vars        = cols,
      n_real      = nrow(real),
      n_syn       = nrow(syn),
      real        = real[cols],
      syn         = syn[cols]
    ),
    class = "flexsynth_diagnostics"
  )
}

#' @export
print.flexsynth_diagnostics <- function(x, ...) {
  cat("<flexsynth_diagnostics>\n")
  cat("  rows        : real", x$n_real, " synthetic", x$n_syn, "\n")
  cat("  variables   :", length(x$vars), "\n\n")

  u <- x$univariate
  cat("Univariate fit (smaller = closer):\n")
  disp <- data.frame(
    variable = u$variable,
    type     = u$type,
    metric   = u$metric,
    distance = round(u$distance, 4),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  cat(sprintf("  mean distance: %.4f   worst: %s (%.4f)\n",
              mean(u$distance, na.rm = TRUE),
              u$variable[which.max(u$distance)], max(u$distance, na.rm = TRUE)))

  if (!is.null(x$correlation)) {
    cr <- x$correlation
    cat(sprintf("\nCorrelation structure (%d numeric vars):\n", length(cr$vars)))
    cat(sprintf("  Frobenius diff: %.4f   mean |diff|: %.4f   max |diff|: %.4f\n",
                cr$frobenius, cr$mean_abs_diff, cr$max_abs_diff))
  }

  if (!is.null(x$pmse)) {
    p <- x$pmse
    cat(sprintf("\nPropensity utility (pMSE, %s):\n", p$method))
    cat(sprintf("  pMSE: %.5f", p$pmse))
    if (!is.na(p$ratio))
      cat(sprintf("   expected: %.5f   ratio: %.2f (1 = indistinguishable)",
                  p$expected, p$ratio))
    cat("\n")
  }
  invisible(x)
}

#' @export
print.flexsynth_diagnostics_list <- function(x, ...) {
  cat("<flexsynth_diagnostics_list>", length(x), "tables\n\n")
  for (nm in names(x)) {
    cat("== ", nm, " ==\n", sep = "")
    print(x[[nm]])
    cat("\n")
  }
  invisible(x)
}

#' Plot utility diagnostics
#'
#' Overlays the real and synthetic marginal distribution of each variable
#' (density for numeric, side-by-side bars for categorical) using base graphics.
#'
#' @param x A `flexsynth_diagnostics` object from [diagnose()].
#' @param vars Optional subset of variables to plot.
#' @param max_panels Maximum number of variables to draw (default 12).
#' @param ... Passed to the underlying plotting calls.
#' @return `x`, invisibly.
#' @export
plot.flexsynth_diagnostics <- function(x, vars = NULL, max_panels = 12L, ...) {
  vv <- vars %||% x$vars
  vv <- intersect(vv, x$vars)
  if (length(vv) > max_panels) vv <- vv[seq_len(max_panels)]

  op <- graphics::par(mfrow = grDevices::n2mfrow(length(vv)),
                      mar = c(3.2, 3.2, 2, 0.8), mgp = c(2, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)

  cr <- grDevices::adjustcolor("#1f77b4", alpha.f = 0.5)   # real
  cs <- grDevices::adjustcolor("#d62728", alpha.f = 0.5)   # synthetic
  for (v in vv) {
    xr <- x$real[[v]]; xs <- x$syn[[v]]
    if (is.numeric(xr)) {
      dr <- stats::density(xr[!is.na(xr)]); ds <- stats::density(xs[!is.na(xs)])
      graphics::plot(range(dr$x, ds$x), range(0, dr$y, ds$y), type = "n",
                     xlab = v, ylab = "density", main = v)
      graphics::polygon(dr$x, dr$y, col = cr, border = NA)
      graphics::polygon(ds$x, ds$y, col = cs, border = NA)
    } else {
      lev <- union(unique(as.character(xr)), unique(as.character(xs)))
      pr <- tabulate(match(as.character(xr), lev), length(lev)) / length(xr)
      ps <- tabulate(match(as.character(xs), lev), length(lev)) / length(xs)
      graphics::barplot(rbind(pr, ps), beside = TRUE, names.arg = lev,
                        col = c(cr, cs), border = NA, main = v, ylab = "prop",
                        las = 2, cex.names = 0.7)
    }
    graphics::legend("topright", c("real", "syn"), fill = c(cr, cs),
                     bty = "n", cex = 0.7)
  }
  invisible(x)
}
