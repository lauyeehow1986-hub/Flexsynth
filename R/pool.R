# Valid inference from synthetic data: combine an analysis fitted on each of the
# m synthetic datasets into one estimate with a standard error that reflects the
# extra variability synthesis introduces. Without this, naive standard errors
# from a single synthetic dataset are (usually) too small.
#
# flexsynth regenerates every column from models fitted on the real data with a
# freshly drawn skeleton, so its output is *fully synthetic* in the Reiter (2003)
# sense. Two combining rules are offered:
#
#  * "synthpop" (default) - the large-sample estimators of Raab, Nowok & Dibben
#    (2016), matching synthpop's summary.fit.synds. The within-synthesis variance
#    (mean of the per-fit variances) is rescaled by k / n (synthetic / original
#    size) and inflated for synthesis: for simple synthesis
#    Tf = vbar (1 + n / (k m)); for proper synthesis Tf = vbar (1 + (n/k + 1)/m).
#    Stable at small m and never negative.
#  * "reiter" - the classic fully-synthetic estimator Tf = (1 + 1/m) b_m - vbar
#    (Reiter 2003), using the between-synthesis variance b_m. Theoretically
#    standard but can be negative (then the standard error is NA) and needs
#    m >= 2, ideally >= 5.
#
# Both report large-sample normal (z) intervals, as synthpop does.

# Combine per-synthesis point estimates and variances into a pooled table.
# `ests` / `vars` are length-m lists of named numeric vectors (point estimates
# and their variances, i.e. squared standard errors) that must share names.
combine_estimates <- function(ests, vars, n, k, proper, rule,
                              population_inference, conf.level) {
  m     <- length(ests)
  terms <- names(ests[[1L]])
  E <- do.call(rbind, lapply(ests, function(e) e[terms]))   # m x p estimates
  U <- do.call(rbind, lapply(vars, function(v) v[terms]))   # m x p variances

  qbar        <- colMeans(E)                 # pooled point estimate
  ubar        <- colMeans(U)                 # mean within-synthesis variance
  vars_scaled <- ubar * (k / n)              # synthpop's `vars`

  if (rule == "synthpop") {
    if (!population_inference) {
      Tf <- vars_scaled
    } else if (isTRUE(proper)) {
      Tf <- vars_scaled * (1 + (n / k + 1) / m)
    } else {
      Tf <- vars_scaled * (1 + (n / k) / m)
    }
  } else {                                   # reiter
    if (m < 2L)
      stop("`rule = \"reiter\"` needs m >= 2 synthetic datasets.", call. = FALSE)
    if (!population_inference) {
      Tf <- vars_scaled
    } else {
      bm <- apply(E, 2L, stats::var)         # between-synthesis variance
      Tf <- (1 + 1 / m) * bm - ubar
    }
  }

  neg <- is.finite(Tf) & Tf < 0
  if (any(neg)) {
    warning(sprintf(paste0("pool_synth(): negative total variance for %s under ",
                           "rule = \"reiter\"; standard error set to NA. Use more ",
                           "synthetic datasets (m) or rule = \"synthpop\"."),
                    paste(terms[neg], collapse = ", ")), call. = FALSE)
    Tf[neg] <- NA_real_
  }

  se <- sqrt(Tf)
  zc <- stats::qnorm(1 - (1 - conf.level) / 2)
  data.frame(
    term      = terms,
    estimate  = qbar,
    std.error = se,
    statistic = qbar / se,
    p.value   = 2 * stats::pnorm(-abs(qbar / se)),
    conf.low  = qbar - zc * se,
    conf.high = qbar + zc * se,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Pull a named estimate vector and a matching variance vector out of whatever the
# analysis returned: a fitted model with coef()/vcov() methods, or an explicit
# list(estimate = , variance = ).
extract_est_var <- function(fit) {
  if (is.list(fit) && all(c("estimate", "variance") %in% names(fit))) {
    e <- fit$estimate; v <- fit$variance
    if (length(e) != length(v))
      stop("analysis: `estimate` and `variance` must have equal length.",
           call. = FALSE)
    nm <- names(e) %||% paste0("term", seq_along(e))
    return(list(est = stats::setNames(as.numeric(e), nm),
                var = stats::setNames(as.numeric(v), nm)))
  }
  e <- tryCatch(stats::coef(fit), error = function(...) NULL)
  V <- tryCatch(stats::vcov(fit), error = function(...) NULL)
  if (is.null(e) || is.null(V))
    stop(paste0("analysis must return a fitted model supporting coef()/vcov() ",
                "(e.g. lm/glm), or a list(estimate = , variance = )."),
         call. = FALSE)
  keep <- !is.na(e)                          # drop aliased (NA) coefficients
  list(est = e[keep], var = diag(V)[names(e)[keep]])
}

#' Valid inference from synthetic data (pooled analysis)
#'
#' Fit an analysis on each of the `m` synthetic datasets in a [synth()] result
#' and combine the results into one estimate whose standard error reflects the
#' extra variability that synthesis introduces. A naive analysis of a single
#' synthetic dataset generally reports standard errors that are too small; this
#' is the correct way to obtain confidence intervals and tests.
#'
#' `analysis` is run once per synthetic dataset and must return either a fitted
#' model supporting `coef()` and `vcov()` (such as [lm()] or [glm()]) or a list
#' with numeric `estimate` and `variance` (squared standard error) vectors. The
#' per-dataset estimates are then combined with the chosen `rule`: `"synthpop"`
#' (Raab, Nowok & Dibben 2016, matching synthpop's `summary.fit.synds`) or
#' `"reiter"` (Reiter 2003 fully-synthetic, `Tf = (1 + 1/m) b_m - vbar`, which
#' needs `m >= 2` and can be negative). Intervals are large-sample normal, as in
#' synthpop.
#'
#' Only Track A ([synth()] without `privacy`) results are supported: differentially
#' private (Track B) inference must additionally account for the DP noise and is
#' out of scope here. Pooling of a [synth_linked()] result is not yet supported;
#' analyse a single table's synthetic frame instead.
#'
#' @param object A `synth_result` from [synth()] (ideally with `m >= 5`).
#' @param analysis A function of one synthetic `data.frame` returning a fitted
#'   model (with `coef()`/`vcov()`) or a `list(estimate =, variance =)`.
#' @param rule Combining rule: `"synthpop"` (default, Raab/Nowok/Dibben 2016) or
#'   `"reiter"` (Reiter 2003 fully-synthetic; needs `m >= 2`).
#' @param population_inference If `TRUE` (default) the standard error targets the
#'   population the real data was drawn from, inflating for synthesis. If `FALSE`
#'   it reports only the (size-rescaled) within-synthesis variance — the standard
#'   error an analyst of the original data would have seen.
#' @param conf.level Confidence level for the intervals (default `0.95`).
#' @param ... Reserved for future use.
#'
#' @return A `flexsynth_pool` object; `$estimates` is a data frame with `term`,
#'   `estimate`, `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`.
#' @seealso [synth_glm()] for the common linear / generalised-linear case.
#' @export
#' @examples
#' d <- data.frame(id = 1:400, x = rnorm(400))
#' d$y <- 1 + 2 * d$x + rnorm(400)
#' res <- synth(d, ~ id, m = 5, seed = 1)
#' pool_synth(res, function(dat) lm(y ~ x, dat))
pool_synth <- function(object, analysis,
                       rule = c("synthpop", "reiter"),
                       population_inference = TRUE,
                       conf.level = 0.95, ...) {
  if (inherits(object, "synth_linked_result"))
    stop(paste0("pool_synth() does not (yet) support linked results; analyse a ",
                "single table's synthetic frame."), call. = FALSE)
  if (!inherits(object, "synth_result"))
    stop("`object` must be a synth_result from synth().", call. = FALSE)
  if (!is.null(object$privacy))
    stop(paste0("pool_synth() is for Track A results. Valid inference from a ",
                "differentially private (Track B) release must account for the ",
                "DP noise and is out of scope."), call. = FALSE)
  if (!is.function(analysis))
    stop("`analysis` must be a function of one synthetic data.frame.", call. = FALSE)
  rule <- match.arg(rule)
  if (!is.logical(population_inference) || length(population_inference) != 1L ||
      is.na(population_inference))
    stop("`population_inference` must be a single TRUE/FALSE.", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L || is.na(conf.level) ||
      conf.level <= 0 || conf.level >= 1)
    stop("`conf.level` must be a single number in (0, 1).", call. = FALSE)

  m <- object$m
  if (rule == "reiter" && m < 2L)
    stop("`rule = \"reiter\"` needs m >= 2 synthetic datasets.", call. = FALSE)

  dsets <- if (m == 1L) list(object$syn) else object$syn

  ev <- lapply(dsets, function(d) extract_est_var(analysis(d)))
  terms <- names(ev[[1L]]$est)
  for (i in seq_along(ev)) {
    if (!identical(names(ev[[i]]$est), terms))
      stop("analysis returned different terms across synthetic datasets.",
           call. = FALSE)
  }
  ests <- lapply(ev, `[[`, "est")
  vars <- lapply(ev, `[[`, "var")

  n <- object$n
  k <- mean(vapply(dsets, nrow, integer(1L)))
  est <- combine_estimates(ests, vars, n = n, k = k,
                           proper = isTRUE(object$proper), rule = rule,
                           population_inference = population_inference,
                           conf.level = conf.level)

  structure(
    list(estimates = est, m = m, n = n, k = k, rule = rule,
         proper = isTRUE(object$proper),
         population_inference = population_inference, conf.level = conf.level,
         call = match.call()),
    class = "flexsynth_pool"
  )
}

#' Pooled generalised-linear analysis of synthetic data
#'
#' Convenience wrapper over [pool_synth()] for the common case: fit the same
#' [glm()] on each synthetic dataset and combine. With the default
#' `family = gaussian()` this is an ordinary linear model.
#'
#' @param object A `synth_result` from [synth()].
#' @param formula A model formula (e.g. `y ~ x1 + x2`).
#' @param family A [family()] for [glm()]; default [gaussian()] (linear model).
#' @param rule,population_inference,conf.level Passed to [pool_synth()].
#' @param ... Further arguments to [glm()].
#'
#' @return A `flexsynth_pool` object (see [pool_synth()]).
#' @export
#' @examples
#' d <- data.frame(id = 1:400, x = rnorm(400))
#' d$y <- rbinom(400, 1, plogis(-0.5 + d$x))
#' res <- synth(d, ~ id, m = 5, seed = 1)
#' synth_glm(res, y ~ x, family = binomial())
synth_glm <- function(object, formula, family = stats::gaussian(),
                      rule = c("synthpop", "reiter"),
                      population_inference = TRUE, conf.level = 0.95, ...) {
  rule <- match.arg(rule)
  pool_synth(object,
             analysis = function(d)
               stats::glm(formula, family = family, data = d, ...),
             rule = rule, population_inference = population_inference,
             conf.level = conf.level)
}

#' @export
print.flexsynth_pool <- function(x, ...) {
  cat("<flexsynth_pool> pooled inference from", x$m, "synthetic dataset(s)\n")
  cat(sprintf("  rule: %s   %s   (n = %d, k = %.0f%s)\n",
              x$rule,
              if (x$population_inference) "population inference" else "sample inference",
              x$n, x$k, if (x$proper) ", proper" else ""))
  est <- x$estimates
  disp <- data.frame(
    term      = est$term,
    estimate  = round(est$estimate, 5),
    std.error = round(est$std.error, 5),
    statistic = round(est$statistic, 3),
    p.value   = signif(est$p.value, 3),
    conf.low  = round(est$conf.low, 5),
    conf.high = round(est$conf.high, 5),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  invisible(x)
}
