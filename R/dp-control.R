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
#' @param structure_frac Budget-efficient structure learning for the flat
#'   `dependence = "tree"` release ([synth()] with three or more variables).
#'   `NULL` (default) measures all \eqn{\binom{d}{2}} pairwise marginals at full
#'   fidelity and reuses them as both the Chow-Liu weights and the tree's
#'   conditionals. A number in `(0, 1)` instead spends only that fraction of the
#'   marginal budget on a cheap all-pairs scan used **solely** to select the tree,
#'   then concentrates the remaining `1 - structure_frac` on re-measuring just the
#'   \eqn{d - 1} chosen edges (plus the `d` one-way marginals) — so the budget
#'   lands on the parameters that survive into the model instead of on pairs that
#'   are discarded. Structure selection tolerates noise well, so a small slice
#'   (e.g. `0.2`–`0.3`) usually sharpens the retained edges, more so as the number
#'   of variables grows. Both passes compose into the same exact
#'   (\eqn{\epsilon}, \eqn{\delta}) budget. Inert for `dependence = "independent"`,
#'   for fewer than three variables (the tree is then trivial), and for
#'   longitudinal / linked DP releases (which keep the single-pass fitter).
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
#' @param longitudinal Linked DP only ([synth_linked()]). Model a child table's
#'   repeated rows as a within-unit time series under DP, instead of exchangeable
#'   records: for such a table the children-per-parent count model doubles as a
#'   trajectory-length model, an initial-state model is measured over each unit's
#'   first (temporally earliest) child row, and a first-order Markov transition
#'   matrix \eqn{P(v_t \mid v_{t-1})} is measured per variable over consecutive
#'   within-unit rows (ordered by the child's own key index). Sensitivities differ
#'   by histogram — initial marginals at the parent's path cap, transitions at
#'   `path_cap[parent] * (branching_cap - 1)` — and fold into the same exact
#'   (\eqn{\epsilon}, \eqn{\delta}) budget. Pass `TRUE` to model every eligible
#'   child table (a child with a branching cap \eqn{\ge} 2 and at least one
#'   variable), or a character vector of child-table names to model only those.
#'   The child's branching cap in `max_rows_per_person` bounds the transition
#'   sensitivity, so it must be \eqn{\ge} 2 for a named table, and the prefix of
#'   that many rows is kept in temporal order. Setting `cross_table = TRUE`
#'   together with a longitudinal model on the same table **combines** them: the
#'   table's initial-state model is cross-conditioned on the synthetic parent
#'   (adding the parent-by-child joints at the same first-row sensitivity as the
#'   initial marginals), and the transition chain then carries that parent
#'   dependence across the trajectory — the parent shapes where a trajectory
#'   starts, while the transitions stay parent-free, so the only extra cost is the
#'   initial-state joints. `FALSE` (default) treats child rows as exchangeable.
#'   Ignored by flat / longitudinal `synth()` releases (which pick the DP Markov
#'   engine straight from the `structure` formula).
#' @param baseline Longitudinal `synth()` only (a `structure` with a nesting
#'   index). Names of **subject-invariant** columns — baseline covariates that do
#'   not change across a person's rows (e.g. birth sex, a baseline measurement).
#'   These are held **exactly constant** within each synthetic unit: they are
#'   modelled once in the initial-state model (so their joint distribution and
#'   their correlation with the first visit are preserved) and then broadcast to
#'   every row, rather than being stepped through a transition matrix that would
#'   let them drift. Declaring a column baseline is public schema knowledge, so it
#'   costs no budget; it also **removes** that column's transition histogram from
#'   the release, sharpening every remaining measurement at the same
#'   (\eqn{\epsilon}, \eqn{\delta}). `NULL` (default) treats every column as
#'   time-varying. Names that match no modelled column are ignored. Ignored by a
#'   flat `synth()` release (one row per unit) and by linked releases.
#' @param transition_order Longitudinal `synth()` only (a `structure` with a
#'   nesting index). Markov order of the within-unit transition model: how many of
#'   a variable's **own** immediately preceding values condition its next value.
#'   `1` (default) is the first-order model \eqn{P(v_t \mid v_{t-1})}; `2` gives
#'   \eqn{P(v_t \mid v_{t-1}, v_{t-2})}, and so on. Higher orders capture momentum
#'   but need a larger row cap: the order must be at most
#'   `max_rows_per_person - 1` (each modelled transition consumes that many prior
#'   rows). A person then contributes at most `cap - order` transition tuples per
#'   variable, so the transition sensitivity — and hence the noise — actually
#'   *drops* with order, at the price of a finer conditioning grid (sparser cells)
#'   and of not modelling the earliest rows' own dynamics separately: rows before
#'   position `order + 1` are generated by marginalising the same measured tensor,
#'   which is post-processing and costs no budget. Ignored by flat / linked
#'   releases.
#' @param transition_cross Longitudinal `synth()` only. Number of **other**
#'   variables (each at lag 1) that additionally condition each variable's
#'   transition — moving from \eqn{P(v_t \mid v_{t-1}, \dots)} to also conditioning
#'   on \eqn{u_{t-1}} for the `transition_cross` most strongly associated
#'   companions `u`. Cross-parents are selected automatically and budget-neutrally
#'   from the pairwise marginals the Chow-Liu tree already measures
#'   (contemporaneous mutual information, a leak-free proxy for lag-1 cross-
#'   predictiveness), so `transition_cross > 0` requires `dependence = "tree"`.
#'   Adding conditioning columns does **not** change the (\eqn{\epsilon},
#'   \eqn{\delta}) budget — a transition tuple still lands in exactly one cell — it
#'   trades budget-free structure for cell sparsity, so keep it small. `0`
#'   (default) conditions each variable on its own past only. Ignored by flat /
#'   linked releases.
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
                       structure_frac = NULL,
                       cross_table = FALSE,
                       longitudinal = FALSE,
                       baseline = NULL,
                       transition_order = 1L,
                       transition_cross = 0L,
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
  if (!is.null(structure_frac)) {
    if (!is.numeric(structure_frac) || length(structure_frac) != 1L ||
        is.na(structure_frac) || structure_frac <= 0 || structure_frac >= 1) {
      stop("`structure_frac` must be NULL or a single number in (0, 1).",
           call. = FALSE)
    }
    structure_frac <- as.numeric(structure_frac)
  }
  if (!is.logical(cross_table) || length(cross_table) != 1L ||
      is.na(cross_table)) {
    stop("`cross_table` must be a single TRUE or FALSE.", call. = FALSE)
  }
  long_ok <- (is.logical(longitudinal) && length(longitudinal) == 1L &&
                !is.na(longitudinal)) ||
    (is.character(longitudinal) && length(longitudinal) >= 1L &&
       !anyNA(longitudinal) && all(nzchar(longitudinal)))
  if (!long_ok) {
    stop(paste0("`longitudinal` must be a single TRUE/FALSE, or a character ",
                "vector of child-table names (for linked DP)."), call. = FALSE)
  }
  if (!is.null(baseline)) {
    if (!is.character(baseline) || anyNA(baseline) || !all(nzchar(baseline))) {
      stop(paste0("`baseline` must be NULL or a character vector of column ",
                  "names (no NA, no empty strings)."), call. = FALSE)
    }
    baseline <- unique(baseline)
  }
  if (!is.numeric(transition_order) || length(transition_order) != 1L ||
      is.na(transition_order) || transition_order < 1 ||
      transition_order != as.integer(transition_order)) {
    stop("`transition_order` must be a single integer >= 1.", call. = FALSE)
  }
  transition_order <- as.integer(transition_order)
  if (!is.numeric(transition_cross) || length(transition_cross) != 1L ||
      is.na(transition_cross) || transition_cross < 0 ||
      transition_cross != as.integer(transition_cross)) {
    stop("`transition_cross` must be a single non-negative integer.",
         call. = FALSE)
  }
  transition_cross <- as.integer(transition_cross)
  if (transition_cross > 0L && dependence == "independent") {
    stop(paste0("`transition_cross` > 0 needs `dependence = \"tree\"`: cross-",
                "variable temporal parents are selected from the pairwise ",
                "marginals that only the tree model measures."), call. = FALSE)
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
      structure_frac = structure_frac,          # NULL or number in (0, 1)
      cross_table = cross_table,
      longitudinal = longitudinal,               # FALSE, TRUE, or table names
      baseline = baseline,                        # NULL or character column names
      transition_order = transition_order,        # integer >= 1 (own-lag order)
      transition_cross = transition_cross,        # integer >= 0 (cross-parents)
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
  if (!is.null(x$structure_frac))
    cat("  structure : budget-efficient (", signif(x$structure_frac, 3),
        " of budget selects the tree)\n", sep = "")
  if (isTRUE(x$cross_table))
    cat("  cross-table: parent-conditioned child vars (linked DP)\n")
  if (isTRUE(x$longitudinal))
    cat("  longitudinal: DP Markov over child rows, all tables (linked DP)\n")
  else if (is.character(x$longitudinal))
    cat("  longitudinal: DP Markov over child rows of",
        paste(x$longitudinal, collapse = ", "), "(linked DP)\n")
  if (!is.null(x$baseline))
    cat("  baseline  : held constant within unit:",
        paste(x$baseline, collapse = ", "), "\n")
  ord <- if (is.null(x$transition_order)) 1L else x$transition_order
  crs <- if (is.null(x$transition_cross)) 0L else x$transition_cross
  if (ord > 1L || crs > 0L)
    cat("  transitions: order ", ord, " + ", crs,
        " cross-parent(s) (longitudinal DP)\n", sep = "")
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
