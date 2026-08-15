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
#' Where the bin edges come from matters, because reading them from the data (its
#' min / max) is itself a data-dependent step that can leak an individual's
#' presence. `domain` controls this:
#' \describe{
#'   \item{`"dp"` (default)}{Rigorous and automatic. Numeric variables named in
#'     `bounds` use those public edges at no cost; any other numeric variable has
#'     its working range **estimated under differential privacy** (a clamp-free
#'     exponential-mechanism quantile at each end), and that estimation spends an
#'     accounted slice `domain_frac` of the budget. The reported
#'     (\eqn{\epsilon}, \eqn{\delta}) is therefore exact end to end — discretisation
#'     adds no unaccounted leakage.}
#'   \item{`"public"`}{Fully data-independent and free: every numeric variable
#'     **must** be given a public range in `bounds` (an error is raised otherwise)
#'     and no budget is spent on the domain. The recommended mode when public
#'     ranges (codebook / physiological limits) are available.}
#'   \item{`"data"`}{The non-rigorous legacy behaviour: bin edges are taken from
#'     the data range with a warning, and that step is **excluded** from the
#'     accounting. Kept only for benchmarking; do not use for a governed release.}
#' }
#' Categorical variables carry their domain in their type: `factor` and `logical`
#' columns use their declared levels (public metadata, no leakage). In the `"dp"`
#' and `"public"` modes a bare `character` column is refused — convert it to a
#' `factor` with its full `levels` so the category set is public.
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
#'   knowledge when a person legitimately has several rows. For a longitudinal DP
#'   release (a `structure` with a nesting index), this must be `>= 2` — it caps
#'   the transition sensitivity — and the per-person prefix of that many rows is
#'   kept in temporal order. For a linked DP release ([synth_linked()]) it is the
#'   maximum children kept per parent — a single integer applied to every child
#'   table, or a **named** vector/list keyed by child-table name (e.g.
#'   `list(admissions = 8, labs = 20)`) — and the root cap is always 1.
#' @param mechanism Noise mechanism: `"laplace"` (pure \eqn{\epsilon}-DP) or
#'   `"gaussian"` (approximate DP with zCDP composition; needs `delta > 0`).
#' @param dependence Dependence structure of the generative model: `"tree"`
#'   (default; a Chow-Liu tree over pairwise marginals) or `"independent"`
#'   (one-way marginals only).
#' @param cross_table Linked DP only ([synth_linked()]). When `TRUE`, a child
#'   table's variables are conditioned on the **synthetic parent's** attributes:
#'   for each child table with a modellable immediate parent, parent-by-child
#'   joint marginals are measured (at the child's per-entity path-cap sensitivity,
#'   the same as a child one-way marginal) and folded into the child's Chow-Liu
#'   structure as fixed context nodes, so each child variable's single strongest
#'   predictor may be a parent variable or another child variable. The synthetic
#'   parent's already-drawn value then conditions the child draw, so cross-table
#'   statistical dependence — not just referential integrity — survives the noise.
#'   Costs extra budget (the added joints compose into the same
#'   (\eqn{\epsilon}, \eqn{\delta})); `FALSE` (default) keeps child variables on
#'   their own within-table marginals. Ignored by flat / longitudinal `synth()`
#'   releases, which have no parent table.
#' @param bins Number of equal-width bins used to discretise each numeric
#'   variable (default 12). Finer grids sharpen one-way marginals but make the
#'   noisy two-way marginals used by `dependence = "tree"` weaker per cell, so a
#'   moderate value usually gives the best correlation fidelity. Means and sums
#'   are largely unaffected (bin contents are decoded uniformly within the bin).
#' @param bounds Optional named list giving `c(lower, upper)` for numeric
#'   variables, used as public, data-independent bin edges. How variables *not*
#'   named here are handled depends on `domain`.
#' @param domain How numeric bin edges are chosen for variables without a public
#'   range in `bounds`: `"dp"` (default) estimates them under differential privacy
#'   and accounts for the cost; `"public"` requires `bounds` for every numeric
#'   variable (error otherwise) and spends no budget on the domain; `"data"` reads
#'   them from the data range with a warning and *excludes* that step from the
#'   accounting (non-rigorous; benchmarking only).
#' @param domain_frac Fraction of the privacy budget spent estimating bin edges
#'   under `domain = "dp"` (default `0.1`). Split evenly across the two
#'   quantile queries per estimated variable and composed with the marginal
#'   measurements, so the total spend is exactly (\eqn{\epsilon}, \eqn{\delta}).
#'   Ignored unless some numeric variable actually needs estimating.
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
                       cross_table = FALSE,
                       bins = 12L,
                       bounds = NULL,
                       domain = c("dp", "public", "data"),
                       domain_frac = 0.1) {
  unit <- match.arg(unit)
  mechanism <- match.arg(mechanism)
  dependence <- match.arg(dependence)
  domain <- match.arg(domain)

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
  if (!is.null(max_rows_per_person)) {
    mrp <- max_rows_per_person
    is_named <- !is.null(names(mrp)) && all(names(mrp) != "") &&
      (is.list(mrp) || is.numeric(mrp)) && length(mrp) >= 1L
    if (is_named) {
      vals <- suppressWarnings(as.numeric(unlist(mrp, use.names = FALSE)))
      ok <- length(vals) == length(mrp) && all(is.finite(vals)) &&
        all(vals >= 1) && all(vals == as.integer(vals))
      if (!ok) {
        stop(paste0("each per-table `max_rows_per_person` must be a positive ",
                    "integer."), call. = FALSE)
      }
      max_rows_per_person <- stats::setNames(as.integer(vals), names(mrp))
    } else if (is.numeric(mrp) && length(mrp) == 1L && !is.na(mrp) &&
               mrp >= 1 && mrp == as.integer(mrp)) {
      max_rows_per_person <- as.integer(mrp)
    } else {
      stop(paste0("`max_rows_per_person` must be NULL, a single positive ",
                  "integer, or a named vector/list of positive integers (one ",
                  "per table, for linked DP)."), call. = FALSE)
    }
  }
  if (!is.logical(cross_table) || length(cross_table) != 1L ||
      is.na(cross_table)) {
    stop("`cross_table` must be a single TRUE or FALSE.", call. = FALSE)
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
  if (!is.numeric(domain_frac) || length(domain_frac) != 1L ||
      is.na(domain_frac) || domain_frac <= 0 || domain_frac >= 1) {
    stop("`domain_frac` must be a single number in (0, 1).", call. = FALSE)
  }

  structure(
    list(
      epsilon = epsilon,
      delta = delta,
      unit = unit,
      max_rows_per_person = max_rows_per_person,   # NULL, integer, or named ints
      mechanism = mechanism,
      dependence = dependence,
      cross_table = cross_table,
      bins = as.integer(bins),
      bounds = bounds,
      domain = domain,
      domain_frac = domain_frac
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
  if (isTRUE(x$cross_table))
    cat("  cross-table: parent-conditioned child vars (linked DP)\n")
  cat("  bins      :", x$bins, "\n")
  cat("  domain    :", x$domain,
      if (x$domain == "dp") paste0("(", signif(x$domain_frac, 3),
                                   " of budget for edges)") else "", "\n")
  mrp <- x$max_rows_per_person
  cat("  max rows/person:",
      if (is.null(mrp)) "1 (auto)"
      else if (!is.null(names(mrp)))
        paste(paste0(names(mrp), "=", mrp), collapse = ", ")
      else mrp, "\n")
  invisible(x)
}
