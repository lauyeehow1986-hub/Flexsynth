#' Differential-privacy controls (Track B)
#'
#' Opt into differentially private synthesis. Passing the result as
#' `synth(..., privacy = dp_control(...))` selects the DP-capable synthesiser
#' and yields a formal (\eqn{\epsilon}, \eqn{\delta}) guarantee. The privacy
#' unit is **person-level** by default: the guarantee protects a whole
#' individual (all of their linked rows across all tables), enforced by bounding
#' each person's contribution before the budget is spent.
#'
#' Note: the DP engine is a later roadmap phase (Track B). This constructor is
#' stable so the interface and validation are settled first.
#'
#' @param epsilon Positive privacy-loss budget. Smaller = more private, less
#'   utility.
#' @param delta Failure probability for approximate DP; `0` requests pure
#'   \eqn{\epsilon}-DP.
#' @param unit Privacy unit. `"person"` (default) protects an individual across
#'   all their rows/tables; `"row"` is weaker and cheaper.
#' @param max_rows_per_person Optional cap on how many rows one person may
#'   contribute (per table); required to bound sensitivity at `unit = "person"`.
#'   `NULL` derives a data-independent default at synthesis time.
#' @param mechanism Noise mechanism, e.g. `"laplace"` or `"gaussian"`.
#'
#' @return An object of class `dp_control` (a validated list).
#' @export
#' @examples
#' dp <- dp_control(epsilon = 1, delta = 1e-6, unit = "person")
#' dp
dp_control <- function(epsilon,
                       delta = 0,
                       unit = c("person", "row"),
                       max_rows_per_person = NULL,
                       mechanism = c("laplace", "gaussian")) {
  unit <- match.arg(unit)
  mechanism <- match.arg(mechanism)

  if (missing(epsilon) || !is.numeric(epsilon) || length(epsilon) != 1L ||
      is.na(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
      delta < 0 || delta >= 1) {
    stop("`delta` must be a single number in [0, 1).", call. = FALSE)
  }
  if (mechanism == "gaussian" && delta <= 0) {
    stop("The gaussian mechanism requires delta > 0.", call. = FALSE)
  }
  if (!is.null(max_rows_per_person) &&
      (!is.numeric(max_rows_per_person) || length(max_rows_per_person) != 1L ||
       max_rows_per_person < 1)) {
    stop("`max_rows_per_person` must be NULL or a single positive number.",
         call. = FALSE)
  }

  structure(
    list(
      epsilon = epsilon,
      delta = delta,
      unit = unit,
      max_rows_per_person = max_rows_per_person,
      mechanism = mechanism
    ),
    class = "dp_control"
  )
}

#' @export
print.dp_control <- function(x, ...) {
  cat("<dp_control> differentially private synthesis (Track B)\n")
  cat("  epsilon  :", x$epsilon, "\n")
  cat("  delta    :", x$delta, "\n")
  cat("  unit     :", x$unit, "\n")
  cat("  mechanism:", x$mechanism, "\n")
  cat("  max rows/person:",
      if (is.null(x$max_rows_per_person)) "auto" else x$max_rows_per_person, "\n")
  invisible(x)
}
