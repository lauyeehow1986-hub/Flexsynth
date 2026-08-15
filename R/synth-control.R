#' Tuning controls for synthesis
#'
#' Collects the fine-grained knobs for the sequential synthesis engine. Strong
#' defaults mean beginners can ignore this; power users can tune the
#' utility / privacy / speed trade-off. Returns a validated `synth_control`
#' object consumed by [synth()] and [synth_linked()].
#'
#' @param visit_sequence Optional character vector giving the order in which
#'   variables are synthesised. `NULL` uses the natural (left-to-right) order.
#' @param predictor_matrix Optional 0/1 matrix marking, for each variable, which
#'   variables may predict it (including cross-table predictors). `NULL` lets
#'   the engine derive a sensible default.
#' @param method Optional per-variable method override; a single string applies
#'   globally. `NULL` uses the call-level `method`.
#' @param smoothing Optional smoothing specification for numeric variables.
#' @param proper Logical; use proper (`TRUE`) vs improper (`FALSE`) synthesis.
#' @param k Optional size of each synthetic dataset. `NULL` matches the input.
#' @param cart,forest Named lists of hyperparameters passed to the CART / random
#'   forest backends. For `forest`, `ntree` (number of trees, default 10) and
#'   `mtry` (predictors tried per tree, default all) are recognised.
#' @param smoothing Numeric-variable smoothing: `NULL` (none), `TRUE` /
#'   `"density"` (kernel-smooth every numeric draw), or a character vector of
#'   variable names to smooth.
#' @param constraint_max_tries Integer; how many times [synth()] may regenerate
#'   while rejection-sampling to satisfy `constraints`. Defaults to 50.
#' @param parallel Logical; enable parallel synthesis where supported.
#'
#' @return An object of class `synth_control` (a validated list).
#' @export
#' @examples
#' ctrl <- synth_control(proper = TRUE, cart = list(minbucket = 5))
#' ctrl
synth_control <- function(visit_sequence = NULL,
                          predictor_matrix = NULL,
                          method = NULL,
                          smoothing = NULL,
                          proper = FALSE,
                          k = NULL,
                          cart = list(),
                          forest = list(),
                          constraint_max_tries = 50L,
                          parallel = FALSE) {
  if (!is.null(visit_sequence) && !is.character(visit_sequence)) {
    stop("`visit_sequence` must be NULL or a character vector.", call. = FALSE)
  }
  if (!is.null(predictor_matrix) && !is.matrix(predictor_matrix)) {
    stop("`predictor_matrix` must be NULL or a matrix.", call. = FALSE)
  }
  if (!is.logical(proper) || length(proper) != 1L || is.na(proper)) {
    stop("`proper` must be a single TRUE/FALSE.", call. = FALSE)
  }
  if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || k < 1)) {
    stop("`k` must be NULL or a single positive number.", call. = FALSE)
  }
  if (!is.list(cart) || !is.list(forest)) {
    stop("`cart` and `forest` must be lists of hyperparameters.", call. = FALSE)
  }
  if (!is.null(smoothing) &&
      !(isTRUE(smoothing) || is.character(smoothing))) {
    stop("`smoothing` must be NULL, TRUE, \"density\", or a character vector.",
         call. = FALSE)
  }
  if (!is.numeric(constraint_max_tries) || length(constraint_max_tries) != 1L ||
      is.na(constraint_max_tries) || constraint_max_tries < 1) {
    stop("`constraint_max_tries` must be a single positive integer.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be a single TRUE/FALSE.", call. = FALSE)
  }

  structure(
    list(
      visit_sequence = visit_sequence,
      predictor_matrix = predictor_matrix,
      method = method,
      smoothing = smoothing,
      proper = proper,
      k = k,
      cart = cart,
      forest = forest,
      constraint_max_tries = as.integer(constraint_max_tries),
      parallel = parallel
    ),
    class = "synth_control"
  )
}

#' @export
print.synth_control <- function(x, ...) {
  fmt <- function(v) if (is.null(v)) "none" else paste(v, collapse = ", ")
  cat("<synth_control>\n")
  cat("  proper           :", x$proper, "\n")
  cat("  k                :", if (is.null(x$k)) "match input" else x$k, "\n")
  cat("  method           :", if (is.null(x$method)) "call-level" else fmt(x$method), "\n")
  cat("  smoothing        :", if (isTRUE(x$smoothing)) "all numeric" else fmt(x$smoothing), "\n")
  cat("  predictor_matrix :", if (is.null(x$predictor_matrix)) "default" else "custom", "\n")
  cat("  cart params      :", length(x$cart), "\n")
  cat("  forest params    :", length(x$forest), "\n")
  cat("  constraint tries :", x$constraint_max_tries, "\n")
  cat("  parallel         :", x$parallel, "\n")
  invisible(x)
}
