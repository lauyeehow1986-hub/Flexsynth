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
#' @param degree Fan-in of the Bayesian network fitted for a flat `synth()`
#'   `dependence = "tree"` release. `1` (default) keeps the Chow-Liu tree exactly
#'   as before (each variable conditions on at most one parent). `2` or more opts
#'   into **PrivBayes' GreedyBayes**: a degree-`k` Bayesian network in which each
#'   variable may condition on up to `degree` of the already-generated variables,
#'   its parent set chosen greedily with the **exponential mechanism** (the parents
#'   whose joint the variable depends on most, crediting the parents\eqn{\to}child
#'   association only). Because a tree gives a variable a single parent, it cannot
#'   represent a variable that genuinely depends on two otherwise-unrelated
#'   predecessors (a "v-structure", e.g. `C` a noisy XOR of independent `A` and
#'   `B`); a degree-2 network takes both parents and recovers it. Construction
#'   measures the `d` one-way marginals plus one `(parents, node)` family joint per
#'   non-root node (`2d - 1` marginals total) and spends a `select_frac` slice on
#'   the `d - 1` greedy picks; every slice composes into the same exact
#'   (\eqn{\epsilon}, \eqn{\delta}) (zCDP for Gaussian, pure \eqn{\epsilon} for
#'   Laplace). The network is forward-sampled ancestrally (no PGM inference), so
#'   like the tree it stays inference-free. `degree` is capped to `d - 1`, and a
#'   high fan-in makes each family a `(degree + 1)`-way histogram — the same
#'   `bins^(degree + 1)` cell-count warning as `treewidth` applies. Requires
#'   `dependence = "tree"`; it is an **alternative** to `structure_frac` and to
#'   `select = "adaptive"` (setting either alongside `degree > 1` is an error), and
#'   is refused on longitudinal / linked releases. Ignored when `degree = 1`.
#' @param select How the marginals that make up the model are chosen, for a flat
#'   [synth()] release. `"fixed"` (default) measures a predetermined set (all
#'   one-way marginals, plus — for `dependence = "tree"` — pairwise marginals),
#'   exactly as described above; every existing path is unchanged. `"adaptive"`
#'   opts into an **AIM-style** iterative selector: after measuring the one-way
#'   marginals it builds the dependency model one marginal at a time, at each step
#'   using the **exponential mechanism** to privately pick the marginal whose true
#'   value the model-so-far fits worst, then measuring that marginal under the main
#'   mechanism. Selection and measurement each spend budget and compose exactly
#'   into the same (\eqn{\epsilon}, \eqn{\delta}) (zCDP for Gaussian, pure
#'   \eqn{\epsilon} for Laplace); the number of rounds is fixed in advance from the
#'   variable count and `treewidth`, so the accounting is data-independent. Unlike
#'   `structure_frac` (which picks a tree from one free noisy scan), adaptive
#'   selection is model-error-guided *and* — at `treewidth >= 2` — can measure
#'   marginals that form loops, capturing higher-order interactions a Chow-Liu tree
#'   structurally cannot. Ignores `dependence` and `structure_frac` (supply
#'   `treewidth` / `select_frac` instead — setting `structure_frac` alongside is an
#'   error). `"aim"` is **Full AIM**: it lifts the adaptive selector's
#'   running-intersection constraint, so the exponential mechanism may pick *loopy*
#'   marginals — a pair between two variables already in the model, the cycle a
#'   junction-tree selector structurally cannot close. Because a loopy set has no
#'   forward junction-tree sampler, `"aim"` reconciles the whole measured set (the
#'   `d` one-way marginals plus the selected pairs) into one graphical model over a
#'   **triangulated** junction tree by **Private-PGM** (belief propagation + mirror
#'   descent; see `estimator`) and samples from that. It runs a data-independent
#'   `min(d (d - 1) / 2, treewidth * (d - 1))` selection rounds, each rejecting any
#'   new pair whose triangulated clique would exceed `treewidth + 1`; at
#'   `treewidth = 1` no loop can be closed, so it reduces to an adaptively-selected
#'   tree. Like `"adaptive"` it uses `treewidth` / `select_frac` (not
#'   `structure_frac` / `degree`), and the reconciliation is budget-neutral. All
#'   selectors are currently **flat-table only**: `"adaptive"` / `"aim"` on a
#'   longitudinal or linked structure is refused.
#' @param treewidth Adaptive selection only. Maximum clique size minus one in the
#'   junction-tree model — the ceiling on interaction order the model can hold.
#'   `1` (default) builds a tree (each measured marginal is a pair; equivalent
#'   *model class* to `dependence = "tree"`, but selected adaptively). `2` lets the
#'   selector measure three-way marginals whose variables form triangles, so
#'   pairwise-invisible three-way structure survives, and `3` measures four-way
#'   cliques (interactions no three-way marginal can see, e.g. a 3-bit parity).
#'   Higher clique orders cost more per cell — a `(w+1)`-way histogram over `bins`
#'   cells per axis is `bins^(w+1)` cells, so the per-cell DP noise grows fast; a
#'   warning is raised when that count is large (`bins` is the numeric proxy;
#'   low-cardinality factor cliques are smaller). Raise it only when interactions
#'   of that order genuinely matter. For `select = "aim"` it bounds the model's
#'   *triangulated* clique size the same way (so `1` forbids loops and `2` permits
#'   triangles). Ignored unless `select = "adaptive"` or `select = "aim"`.
#' @param anneal Adaptive selection only. `FALSE` (default) uses the fixed
#'   `d - treewidth`-round schedule with a uniform per-round budget. `TRUE` opts
#'   into **AIM-style budget annealing** — a data-adaptive round schedule: the
#'   `d` one-way marginals take a fair fixed share, then the clique measurements
#'   and their exponential-mechanism selections start at a small per-round quantum
#'   (large noise) and **double** whenever a round's measured signal fails to beat
#'   its noise floor (AIM's \eqn{\sigma}-halving rule). After the mandatory
#'   `d - treewidth` spanning cliques — which guarantee every variable is covered,
#'   keeping the sampler inference-free — any **surplus** budget is spent on extra
#'   rounds that re-measure the worst-fit existing clique (inverse-variance
#'   combined, so \eqn{\rho} / \eqn{\epsilon} adds exactly). The final round
#'   absorbs the exact remainder, so the total spend is still exactly
#'   (\eqn{\epsilon}, \eqn{\delta}) — now over a **variable** number of rounds
#'   chosen from the data. Because the model stays a spanning junction tree (no
#'   PGM inference), refinement can only sharpen the measured cliques, not add
#'   loopy marginals; it is most useful when there are few variables and the fixed
#'   schedule would otherwise leave budget on the table. Ignored unless
#'   `select = "adaptive"`; setting it without adaptive selection is an error.
#' @param select_frac Adaptive selection only. Fraction of the marginal budget
#'   spent on the private **selection** (the exponential-mechanism rounds); the
#'   remaining `1 - select_frac` is spent **measuring** the chosen marginals (and
#'   the one-way marginals). Default `0.25`. Selection tolerates noise well, so a
#'   modest slice usually gives the best fidelity (the same intuition as
#'   `structure_frac`); both slices compose into the same exact
#'   (\eqn{\epsilon}, \eqn{\delta}). Ignored unless `select = "adaptive"` or
#'   `select = "aim"`.
#' @param estimator How the measured marginals are turned into the generative
#'   model, for a flat [synth()] tree or adaptive release. `"local"` (default)
#'   is every existing path unchanged: each marginal is used **locally** — the
#'   tree takes the root's one-way and each edge's raw noisy 2-way as
#'   \eqn{P(\text{child} \mid \text{parent})}, the adaptive junction tree each
#'   clique's own array — so the other one-way marginals, and the disagreement
#'   between overlapping noisy marginals, are discarded. `"pgm"` opts into
#'   **Private-PGM inference** (McKenna et al.'s MST / graphical-model
#'   estimation): the *whole* set of measured marginals is reconciled into the
#'   single graphical-model distribution that best fits all of them at once (least
#'   squares), by belief propagation on a junction tree of the measured cliques
#'   plus entropic mirror descent, and the model is sampled from that. Because a
#'   tree (or the adaptive junction tree) has bounded treewidth, the inference is
#'   exact and cheap. Reconciliation is **pure post-processing** of the
#'   already-privatised marginals, so it spends **no extra budget** — the
#'   (\eqn{\epsilon}, \eqn{\delta}) is identical to the same release with
#'   `"local"`; only the fitted model changes. It denoises (overlapping marginals
#'   are made mutually consistent) and lets the otherwise-discarded one-way
#'   marginals actually constrain the model, usually sharpening both the marginals
#'   and the conditionals at the same budget. Available for the flat
#'   `dependence = "tree"` release (degree-1) and for `select = "adaptive"`; it is
#'   an **alternative** to `structure_frac`, to `degree > 1`, and to
#'   `anneal = TRUE` (setting any of those alongside `estimator = "pgm"` is an
#'   error), and there is nothing to reconcile in `dependence = "independent"`.
#'   Refused on longitudinal / linked releases (flat-table only). Ignored (no-op)
#'   for `"local"`. `select = "aim"` always reconciles with Private-PGM by
#'   construction (a loopy marginal set has no local sampler), so `estimator` is
#'   moot there.
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
#' @param baseline Longitudinal releases. Names of **subject-invariant** columns —
#'   baseline covariates that do not change across a unit's rows (e.g. birth sex, a
#'   baseline measurement). These are held **exactly constant** within each
#'   synthetic unit: they are modelled once in the initial-state model (so their
#'   joint distribution and their correlation with the first visit are preserved)
#'   and then broadcast to every row, rather than being stepped through a
#'   transition matrix that would let them drift. Declaring a column baseline is
#'   public schema knowledge, so it costs no budget; it also **removes** that
#'   column's transition histogram from the release, sharpening every remaining
#'   measurement at the same (\eqn{\epsilon}, \eqn{\delta}). Applies to a
#'   longitudinal `synth()` release (a `structure` with a nesting index) and, per
#'   table, to any **longitudinally-modelled linked child** (the names are matched
#'   against each such child's own columns). `NULL` (default) treats every column
#'   as time-varying. Names that match no modelled column are ignored. Ignored by a
#'   flat `synth()` release (one row per unit) and by non-longitudinal linked
#'   child tables.
#' @param transition_order Longitudinal releases. Markov order of the within-unit
#'   transition model: how many of a variable's **own** immediately preceding
#'   values condition its next value. `1` (default) is the first-order model
#'   \eqn{P(v_t \mid v_{t-1})}; `2` gives \eqn{P(v_t \mid v_{t-1}, v_{t-2})}, and so
#'   on. Higher orders capture momentum but need a larger row cap: the order must be
#'   at most `max_rows_per_person - 1` — for a linked release, one less than the
#'   **branching cap** of each longitudinally-modelled child (validated per table).
#'   A unit then contributes at most `cap - order` transition tuples per variable,
#'   so the transition sensitivity — and hence the noise — actually *drops* with
#'   order, at the price of a finer conditioning grid (sparser cells) and of not
#'   modelling the earliest rows' own dynamics separately: rows before position
#'   `order + 1` are generated by marginalising the same measured tensor, which is
#'   post-processing and costs no budget. Applies to a longitudinal `synth()`
#'   release and, uniformly, to every longitudinally-modelled linked child; ignored
#'   by flat releases and non-longitudinal linked children.
#' @param transition_cross Longitudinal releases. Number of **other** variables
#'   (each at lag 1) that additionally condition each variable's transition —
#'   moving from \eqn{P(v_t \mid v_{t-1}, \dots)} to also conditioning on
#'   \eqn{u_{t-1}} for the `transition_cross` most strongly associated companions
#'   `u`. Cross-parents are selected automatically and budget-neutrally from the
#'   pairwise marginals the Chow-Liu tree already measures (contemporaneous mutual
#'   information, a leak-free proxy for lag-1 cross-predictiveness) — for a linked
#'   child that means its own within-child pairwise marginals — so
#'   `transition_cross > 0` requires `dependence = "tree"`. Adding conditioning
#'   columns does **not** change the (\eqn{\epsilon}, \eqn{\delta}) budget — a
#'   transition tuple still lands in exactly one cell — it trades budget-free
#'   structure for cell sparsity, so keep it small. `0` (default) conditions each
#'   variable on its own past only. Applies to a longitudinal `synth()` release and
#'   to every longitudinally-modelled linked child; ignored by flat releases and
#'   non-longitudinal linked children.
#' @param transition_parent Linked DP only ([synth_linked()]), longitudinal
#'   children. Number of **immediate-parent** attributes that additionally
#'   condition each time-varying child variable's transition. Today
#'   `cross_table = TRUE` cross-conditions only a longitudinal child's
#'   *initial state* on the synthetic parent; the parent's influence then rides
#'   the own-lag chain and decays. `transition_parent = p` instead re-injects the
#'   (subject-invariant) parent attributes into the transition tensor at **every**
#'   step — moving from \eqn{P(v_t \mid v_{t-1}, \dots)} to also conditioning on
#'   the parent's \eqn{w} for the `p` parent attributes `w` most strongly
#'   associated with `v` — so parent \eqn{\to} child dependence stays anchored
#'   across the whole trajectory rather than washing out. The parents are selected
#'   automatically and **budget-neutrally** from the parent-by-child joints that
#'   the cross-conditioned initial state already measures, so it costs nothing in
#'   the (\eqn{\epsilon}, \eqn{\delta}) budget (a transition tuple still lands in
#'   exactly one cell — same sensitivity, same histogram count), trading
#'   budget-free structure for cell sparsity. Because it reuses those joints,
#'   `transition_parent > 0` **requires `cross_table = TRUE`** for the child (an
#'   error is raised otherwise). `0` (default) leaves the transitions parent-free
#'   (the initial-state cross-conditioning is unaffected). Composes with
#'   `baseline`, `transition_order` and `transition_cross`. Ignored by flat
#'   `synth()` releases and non-longitudinal linked children, which have no parent
#'   trajectory to condition.
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
                       degree = 1L,
                       select = c("fixed", "adaptive", "aim"),
                       treewidth = 1L,
                       select_frac = 0.25,
                       anneal = FALSE,
                       estimator = c("local", "pgm"),
                       cross_table = FALSE,
                       longitudinal = FALSE,
                       baseline = NULL,
                       transition_order = 1L,
                       transition_cross = 0L,
                       transition_parent = 0L,
                       bins = 12L,
                       bounds = NULL,
                       domain = c("dp", "public", "data"),
                       domain_frac = 0.1) {
  unit <- match.arg(unit)
  mechanism <- match.arg(mechanism)
  dependence <- match.arg(dependence)
  select <- match.arg(select)
  estimator <- match.arg(estimator)
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
  if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) ||
      degree < 1 || degree != as.integer(degree)) {
    stop("`degree` must be a single integer >= 1.", call. = FALSE)
  }
  degree <- as.integer(degree)
  if (degree > 1L) {
    if (dependence != "tree") {
      stop(paste0("`degree` > 1 builds a Bayesian network and needs ",
                  "`dependence = \"tree\"` (a tree is a degree-1 network)."),
           call. = FALSE)
    }
    if (!is.null(structure_frac)) {
      stop(paste0("`degree` > 1 and `structure_frac` are alternative structure ",
                  "learners for the tree model; set structure_frac = NULL."),
           call. = FALSE)
    }
    if (select %in% c("adaptive", "aim")) {
      stop(paste0("`degree` > 1 (a Bayesian network) and `select = \"", select,
                  "\"` are alternative structure searches; pick one."),
           call. = FALSE)
    }
  }
  if (!is.numeric(treewidth) || length(treewidth) != 1L || is.na(treewidth) ||
      treewidth < 1 || treewidth != as.integer(treewidth)) {
    stop("`treewidth` must be a single integer >= 1.", call. = FALSE)
  }
  treewidth <- as.integer(treewidth)
  if (!is.numeric(select_frac) || length(select_frac) != 1L ||
      is.na(select_frac) || select_frac <= 0 || select_frac >= 1) {
    stop("`select_frac` must be a single number in (0, 1).", call. = FALSE)
  }
  select_frac <- as.numeric(select_frac)
  if (!is.logical(anneal) || length(anneal) != 1L || is.na(anneal)) {
    stop("`anneal` must be a single TRUE or FALSE.", call. = FALSE)
  }
  if (select == "aim" && isTRUE(anneal)) {
    stop(paste0("annealed AIM (`select = \"aim\"` with `anneal = TRUE`) is not ",
                "yet supported; use anneal = FALSE (a fixed-round schedule)."),
         call. = FALSE)
  }
  if (isTRUE(anneal) && select != "adaptive") {
    stop(paste0("`anneal = TRUE` needs `select = \"adaptive\"`: it anneals the ",
                "adaptive selector's per-round budget over a data-adaptive ",
                "round schedule."), call. = FALSE)
  }
  if (select %in% c("adaptive", "aim")) {
    if (!is.null(structure_frac)) {
      stop(paste0("`structure_frac` applies to the fixed tree fitter; `select = \"",
                  select, "\"` uses `select_frac` (and `treewidth`) instead - set ",
                  "structure_frac = NULL."), call. = FALSE)
    }
  }
  # Private-PGM reconciliation is a post-processing estimator that reconciles the
  # measured marginals of a tree or an adaptive junction tree into one consistent
  # model. It is an alternative to the tree's other structure knobs and to the
  # annealed refiner, and is defined for a degree-1 tree (not a Bayesian network).
  if (estimator == "pgm") {
    if (isTRUE(anneal)) {
      stop(paste0("`estimator = \"pgm\"` and `anneal = TRUE` are alternative ",
                  "post-measurement refiners (PGM reconciles the measured ",
                  "marginals; annealing re-measures the worst-fit clique); pick ",
                  "one."), call. = FALSE)
    }
    if (degree > 1L) {
      stop(paste0("`estimator = \"pgm\"` reconciles a Chow-Liu tree (degree 1); ",
                  "the Bayesian-network `degree` > 1 is a different structure. ",
                  "Set degree = 1."), call. = FALSE)
    }
    if (!is.null(structure_frac)) {
      stop(paste0("`estimator = \"pgm\"` and `structure_frac` both restructure ",
                  "the tree release; they are not (yet) combined. Set ",
                  "structure_frac = NULL."), call. = FALSE)
    }
    if (select != "adaptive" && dependence != "tree") {
      stop(paste0("`estimator = \"pgm\"` reconciles a tree or an adaptive ",
                  "junction tree; it needs `dependence = \"tree\"` or ",
                  "`select = \"adaptive\"` (there is nothing to reconcile in the ",
                  "independent model)."), call. = FALSE)
    }
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
  if (!is.numeric(transition_parent) || length(transition_parent) != 1L ||
      is.na(transition_parent) || transition_parent < 0 ||
      transition_parent != as.integer(transition_parent)) {
    stop("`transition_parent` must be a single non-negative integer.",
         call. = FALSE)
  }
  transition_parent <- as.integer(transition_parent)
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 2 ||
      bins != as.integer(bins)) {
    stop("`bins` must be a single integer >= 2.", call. = FALSE)
  }
  # A high-order family makes each measured histogram a (order + 1)-way table;
  # over `bins` cells per axis that is bins^(order + 1) cells, and past a point
  # the per-cell noise swamps the signal. The order is the treewidth for adaptive
  # junction cliques, or the degree for a Bayesian network. Warn on that (bins is
  # the numeric-variable proxy; low-cardinality factor families are smaller).
  fam_order <- if (select %in% c("adaptive", "aim")) treewidth
               else if (degree > 1L) degree else 0L
  if (fam_order > 0L) {
    max_cells <- as.numeric(bins)^(fam_order + 1L)
    if (max_cells > 5e4) {
      what <- if (select %in% c("adaptive", "aim")) c("treewidth", "clique")
              else c("degree", "family")
      warning(sprintf(paste0("%s %d makes each %s up to bins^%d = %.0f cells; ",
                             "at that width the per-cell DP noise can dominate ",
                             "the signal (bins is the numeric proxy; ",
                             "low-cardinality factors are smaller)."),
                      what[1L], fam_order, what[2L], fam_order + 1L, max_cells),
              call. = FALSE)
    }
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
      degree = degree,                          # integer >= 1 (Bayesian network)
      select = select,                          # "fixed", "adaptive", or "aim"
      treewidth = treewidth,                    # integer >= 1 (adaptive only)
      select_frac = select_frac,                # number in (0, 1) (adaptive only)
      anneal = anneal,                          # logical (adaptive only)
      estimator = estimator,                    # "local" or "pgm" (reconcile)
      cross_table = cross_table,
      longitudinal = longitudinal,               # FALSE, TRUE, or table names
      baseline = baseline,                        # NULL or character column names
      transition_order = transition_order,        # integer >= 1 (own-lag order)
      transition_cross = transition_cross,        # integer >= 0 (cross-parents)
      transition_parent = transition_parent,      # integer >= 0 (parent attrs)
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
  if (!is.null(x$degree) && x$degree > 1L)
    cat("  network   : degree-", x$degree, " Bayesian network (GreedyBayes; ",
        signif(x$select_frac, 3), " of budget selects parents)\n", sep = "")
  if (identical(x$select, "adaptive"))
    cat("  select    : adaptive AIM-style, treewidth ", x$treewidth,
        if (isTRUE(x$anneal)) " (annealed, data-adaptive round schedule; "
        else " (", signif(x$select_frac, 3), " of budget selects marginals)\n",
        sep = "")
  if (identical(x$select, "aim"))
    cat("  select    : Full AIM (loopy marginals + Private-PGM), treewidth ",
        x$treewidth, " (", signif(x$select_frac, 3),
        " of budget selects marginals)\n", sep = "")
  if (identical(x$estimator, "pgm"))
    cat("  estimator : Private-PGM reconciliation",
        "(belief propagation + mirror descent; budget-neutral)\n")
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
  tpar <- if (is.null(x$transition_parent)) 0L else x$transition_parent
  if (ord > 1L || crs > 0L || tpar > 0L)
    cat("  transitions: order ", ord, " + ", crs, " cross-parent(s)",
        if (tpar > 0L) paste0(" + ", tpar, " parent-attr(s)") else "",
        " (longitudinal DP)\n", sep = "")
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
