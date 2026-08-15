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
#'   forest backends.
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
      parallel = parallel
    ),
    class = "synth_control"
  )
}

#' @export
print.synth_control <- function(x, ...) {
  cat("<synth_control>\n")
  cat("  proper       :", x$proper, "\n")
  cat("  k            :", if (is.null(x$k)) "match input" else x$k, "\n")
  cat("  parallel     :", x$parallel, "\n")
  cat("  method       :", if (is.null(x$method)) "call-level" else x$method, "\n")
  cat("  cart params  :", length(x$cart), "\n")
  cat("  forest params:", length(x$forest), "\n")
  invisible(x)
}
