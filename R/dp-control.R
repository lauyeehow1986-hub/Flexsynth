#' Differential-privacy controls (Track B)
#'
#' Opt into differentially private synthesis. Passing the result as
#' `synth(..., privacy = dp_control(...))` selects the DP synthesiser and yields
#' a formal (\eqn{\epsilon}, \eqn{\delta}) guarantee. The privacy unit is
#' **person-level** by default: the guarantee protects a whole individual (all of
#' their rows), enforced by bounding each person's contribution before any budget
#' is spent.
#'
#' The DP engine (Track B) is a *marginal-based* synthesiser in the
#' PrivBayes / MST lineage. Continuous variables are discretised into a public
#' grid; low-order marginals are measured under the chosen noise mechanism with
#' correct budget composition; and synthetic records are drawn from the resulting
#' model. With `dependence = "tree"` a Chow-Liu tree of pairwise dependencies is
#' learned from the same noisy marginals (no extra budget) so second-order
#' structure is retained; `dependence = "independent"` keeps only one-way
#' marginals. All noise calibration and composition is reported back on the
#' result (see the accounting printed by [synth()]).
#'
#' Note on bin edges and category sets: unless `bounds` is supplied, numeric bin
#' edges are taken from the data range and categorical domains from the observed
#' levels. That is a small data-dependent step which is **not** included in the
#' (\eqn{\epsilon}, \eqn{\delta}) accounting; supply `bounds` (and factor columns
#' with their full `levels`) from public knowledge to make the release fully
#' data-independent. A warning is emitted when edges are derived from the data.
#'
#' @param epsilon Positive privacy-loss budget. Smaller = more private, less
#'   utility.
#' @param delta Failure probability for approximate DP; `0` requests pure
#'   \eqn{\epsilon}-DP (only the Laplace mechanism supports `delta = 0`).
#' @param unit Privacy unit. `"person"` (default) protects an individual across
#'   all their rows; `"row"` is weaker and cheaper (each row independent).
#' @param max_rows_per_person Cap on how many rows one person may contribute,
#'   used to bound sensitivity at `unit = "person"`. `NULL` (default) uses `1`
#'   (each person contributes at most one row); rows beyond the cap are dropped
#'   by per-person subsampling before measuring. Set this from public domain
#'   knowledge when a person legitimately has several rows.
#' @param mechanism Noise mechanism: `"laplace"` (pure \eqn{\epsilon}-DP) or
#'   `"gaussian"` (approximate DP with zCDP composition; needs `delta > 0`).
#' @param dependence Dependence structure of the generative model: `"tree"`
#'   (default; a Chow-Liu tree over pairwise marginals) or `"independent"`
#'   (one-way marginals only).
#' @param bins Number of equal-width bins used to discretise each numeric
#'   variable (default 12). Finer grids sharpen one-way marginals but make the
#'   noisy two-way marginals used by `dependence = "tree"` weaker per cell, so a
#'   moderate value usually gives the best correlation fidelity. Means and sums
#'   are largely unaffected (bin contents are decoded uniformly within the bin).
#' @param bounds Optional named list giving `c(lower, upper)` for numeric
#'   variables, used as public, data-independent bin edges. Variables not named
#'   fall back to the observed data range (with a warning).
#'
#' @return An object of class `dp_control` (a validated list).
#' @export
#' @examples
#' dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian")
#' dp
dp_control <- function(epsilon,
                       delta = 0,
                       unit = c("person", "row"),
                       max_rows_per_person = NULL,
                       mechanism = c("laplace", "gaussian"),
                       dependence = c("tree", "independent"),
                       bins = 12L,
                       bounds = NULL) {
  unit <- match.arg(unit)
  mechanism <- match.arg(mechanism)
  dependence <- match.arg(dependence)

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
       is.na(max_rows_per_person) || max_rows_per_person < 1 ||
       max_rows_per_person != as.integer(max_rows_per_person))) {
    stop("`max_rows_per_person` must be NULL or a single positive integer.",
         call. = FALSE)
  }
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 2 ||
      bins != as.integer(bins)) {
    stop("`bins` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.null(bounds)) {
    if (!is.list(bounds) || is.null(names(bounds)) || any(names(bounds) == "")) {
      stop("`bounds` must be a named list of numeric c(lower, upper) pairs.",
           call. = FALSE)
    }
    ok <- vapply(bounds, function(b)
      is.numeric(b) && length(b) == 2L && all(is.finite(b)) && b[1L] < b[2L],
      logical(1))
    if (!all(ok)) {
      stop("each entry of `bounds` must be numeric c(lower, upper) with lower < upper.",
           call. = FALSE)
    }
  }

  structure(
    list(
      epsilon = epsilon,
      delta = delta,
      unit = unit,
      max_rows_per_person =
        if (is.null(max_rows_per_person)) NULL else as.integer(max_rows_per_person),
      mechanism = mechanism,
      dependence = dependence,
      bins = as.integer(bins),
      bounds = bounds
    ),
    class = "dp_control"
  )
}

#' @export
print.dp_control <- function(x, ...) {
  cat("<dp_control> differentially private synthesis (Track B)\n")
  cat("  epsilon   :", x$epsilon, "\n")
  cat("  delta     :", x$delta, "\n")
  cat("  unit      :", x$unit, "\n")
  cat("  mechanism :", x$mechanism, "\n")
  cat("  dependence:", x$dependence, "\n")
  cat("  bins      :", x$bins, "\n")
  cat("  max rows/person:",
      if (is.null(x$max_rows_per_person)) "1 (auto)" else x$max_rows_per_person, "\n")
  invisible(x)
}
