# flexsynth 0.2.0

## New features

* **Valid inference from synthetic data.** `pool_synth()` fits an analysis on
  each of the `m` synthetic datasets and combines the results with a standard
  error that reflects the extra variability synthesis introduces — a single
  synthetic dataset analysed naively gives standard errors that are too small.
  Two fully-synthetic combining rules are provided (`rule=`): the large-sample
  synthpop estimators of Raab, Nowok & Dibben (2016) (default, matching
  `synthpop::summary.fit.synds`) and the classic Reiter (2003) estimator. The
  analysis may be any model with `coef()`/`vcov()` methods or a plain
  `list(estimate=, variance=)`. `synth_glm()` is a convenience wrapper for the
  linear / generalised-linear case. `population_inference` toggles between
  population- and sample-level standard errors. Track B (DP) results are refused
  (DP inference must additionally account for the noise). The `synth_result`
  now records whether `proper` synthesis was used, which selects the correct
  variance formula.

## Bug fixes

* **`method = "norm"` / `"normrank"` no longer crash on missing data.** The
  linear-model fit dropped `NA` rows via `model.matrix()` while keeping the
  full-length response, causing `"'qr' and 'y' must have the same number of
  rows"`. The fit now uses complete cases, and prediction propagates a missing
  value where a synthetic predictor is itself missing rather than erroring.
  (A dedicated missingness *model* remains future work.)

# flexsynth 0.1.2

* Documentation: `?synth_linked` previously stated, as flat limitations, that
  under differential privacy (Track B) cross-table statistical conditioning "is
  not preserved" and within-table longitudinal structure "is not modelled". That
  text predated the `dp_control(cross_table = ...)` and
  `dp_control(longitudinal = ...)` opt-ins, which do model both (parent-conditioned
  child variables and a within-unit DP Markov trajectory over a child's repeated
  rows), and so understated the package. The roxygen now documents the default
  (own within-table marginals; exchangeable rows) and both opt-ins accurately.
  Behaviour is unchanged.

# flexsynth 0.1.1

* **Bug fix: `method = "ctree"` on a bare `character` target.** The `partykit`
  conditional-inference-tree method fitted the tree with the raw response column,
  which tripped `partykit:::.y2infl()` with "unknown response class" whenever the
  synthesised column was a plain `character` vector (`rpart`-backed `cart`
  tolerated it, so this only surfaced under `ctree`). The response is now coerced
  to a factor for the fit while draws still return the original values. This was
  caught by win-builder R-devel, where `partykit` is available and the guarded
  test actually runs.
* Documentation: the README link to the roadmap now points at its GitHub URL
  (the `docs/` directory is not shipped in the package), and an `inst/WORDLIST`
  records the intentional British spellings used in the DESCRIPTION.

# flexsynth 0.1.0

First public release: a flexible sequential-synthesis engine for nested,
longitudinal and linked multi-table data (Track A, high-utility default) with an
opt-in differentially private track (Track B) providing a formal person-level
(\eqn{\epsilon}, \eqn{\delta}) guarantee. The notes below summarise what shipped.


* **Track B: DP set-union discovers `character` category sets under
  `domain = "dp"`.** A bare `character` column no longer has to be pre-converted to
  a `factor`: under the default `domain = "dp"`, its present category set is
  discovered privately with a **stability histogram** — each present category's
  count gets Laplace noise and is kept only if it clears a threshold that hides any
  category a single person could have created; rare/unique categories fold into an
  appended `"(other)"` catch-all. The cost is funded from the same `domain_frac`
  slice that numeric bin-edge estimation uses, and composes into the exact
  (\eqn{\epsilon}, \eqn{\delta}). Because a threshold cannot hide a lone category's
  presence at `delta = 0`, discovery **needs `delta > 0`**; a pure-\eqn{\epsilon}
  release still refuses `character` (pointing to `delta > 0` or public `factor`
  levels). Flat `synth()` only (linked / longitudinal still require public levels).
  New internal `dp_discover_categories()`; releases with no `character` column are
  byte-identical to before.

* **Track B: model-projection candidate scoring is now the Full AIM default.**
  `dp_control(select = "aim")` now scores each candidate marginal in AIM's
  exponential-mechanism rounds against the **current reconciled Private-PGM model's
  own marginal** over that pair (`scoring = "model"`), rather than the product of its
  one-way marginals — completing AIM's actual quality function as the faithful
  default. A loopy pair the model already explains through the marginals measured so
  far no longer looks surprising, so the budget is steered to the genuinely worst-fit
  interaction instead of re-selecting one the model already fits. The model reference
  is read from the **already-privatised** marginals (reconciled each round and
  projected onto the candidate — a not-yet-measured pair projects across cliques of
  the junction tree, via the internal `dp_pgm_project()`), so it is **pure
  post-processing**: the exponential mechanism's sensitivity and the exact
  (\eqn{\epsilon}, \eqn{\delta}) are unchanged; only which marginals get selected
  changes. It costs a reconciliation per selection round and composes with
  `anneal = TRUE`. **This can change the output of an existing `select = "aim"`
  release** (a different marginal set may be chosen); the previous one-way-product
  reference remains available as the cheaper opt-out via
  `dp_control(select = "aim", scoring = "independence")`. The new default is expressed
  through `scoring = "auto"` (model for `"aim"`, independence for every other
  selector). Flat-table only.

* **Track B: annealed Full AIM.** `dp_control(select = "aim", anneal = TRUE)` runs
  Full AIM on the same **data-adaptive \eqn{\sigma}-halving schedule** that
  `select = "adaptive"` already offers, instead of the fixed
  `min(d(d-1)/2, treewidth·(d-1))` rounds. The `d` one-way marginals take a fair
  fixed share; the pairwise measurements and their exponential-mechanism selections
  start at a small per-round quantum (large noise) and **double** whenever a round's
  measured signal fails to beat its noise floor. A baseline of treewidth-capped
  *new* loopy pairs is selected first, then any **surplus** budget re-measures the
  worst-fit already-measured pair (inverse-variance combined). The measurement and
  selection pools stay strictly proportional through every doubling and reserve, so
  they deplete in lockstep and the final round absorbs the exact remainder — the
  total spend is still exactly (\eqn{\epsilon}, \eqn{\delta}), now over a **variable**
  number of rounds chosen from the data. As in the fixed-round case the annealed set
  is triangulated and reconciled with **Private-PGM**, so refinement sharpens loopy
  marginals (not just tree cliques). Flat-table only.

* **Track B: Full AIM.** `dp_control(select = "aim")` lifts the running-intersection
  constraint on the `select = "adaptive"` junction-tree selector, so the
  exponential mechanism may pick *loopy* marginals — a pairwise marginal between two
  variables already in the model, the cycle a junction-tree selector structurally
  cannot close (and that no tree captures at any budget). Because a loopy set has no
  forward junction-tree sampler, AIM reconciles the whole measured set (the `d`
  one-way marginals plus the selected pairs) into one graphical model over a
  **triangulated** junction tree via **Private-PGM** (belief propagation + entropic
  mirror descent) and samples from that — following McKenna et al.'s *AIM*. It runs
  a data-independent `min(d(d-1)/2, treewidth·(d-1))` selection rounds, each
  rejecting any new pair whose triangulated clique would exceed `treewidth + 1` (so
  the model's treewidth is bounded and the inference stays exact); at
  `treewidth = 1` no loop can close, so it reduces to an adaptively-selected tree.
  Selection and measurement compose into exactly the same (\eqn{\epsilon},
  \eqn{\delta}) as any other flat slice, and the reconciliation is budget-neutral.
  Like `"adaptive"` it uses `treewidth` / `select_frac` (not `structure_frac` /
  `degree` / `anneal`, which are refused), and is flat-table only (refused on
  longitudinal / linked releases).

* **Track B: Private-PGM inference step.** `dp_control(estimator = "pgm")` adds the
  reconciliation the other (PGM-free) Track B samplers deliberately omit. Instead
  of using each measured marginal *locally* (the tree takes each edge's raw noisy
  2-way as \eqn{P(\text{child} \mid \text{parent})} and discards both the other
  one-way marginals and the disagreement between overlapping noisy marginals), it
  reconciles the *whole* measured set into the single graphical-model distribution
  that best fits all of it at once (least squares) — following McKenna et al.'s
  Private-PGM / MST — by belief propagation on a junction tree of the measured
  cliques plus **entropic mirror descent**, then samples from that. A tree (or the
  adaptive junction tree) has bounded treewidth, so the inference is exact and
  cheap. Reconciliation is **pure post-processing** of the already-privatised
  marginals, so it spends **no extra budget**: the (\eqn{\epsilon}, \eqn{\delta})
  is identical to the same release with `"local"`; only the fitted model changes.
  It denoises (overlapping marginals are made mutually consistent) and lets the
  otherwise-discarded one-way marginals constrain the model. Available for the flat
  `dependence = "tree"` release (degree-1) and for `select = "adaptive"`; an
  **alternative** to `structure_frac`, `degree > 1` and `anneal = TRUE`; refused on
  longitudinal / linked releases. `"local"` (default) leaves every existing path
  byte-identical.

* **Track B: PrivBayes degree>1 Bayesian networks.** `dp_control(dependence =
  "tree", degree = k)` generalises the Chow-Liu tree (a degree-1 Bayesian
  network) to **GreedyBayes**: a degree-`k` network in which each variable may
  condition on up to `k` of the already-generated variables. Each parent set is
  chosen greedily with the **exponential mechanism**, scored by the
  parents\eqn{\to}node association only (so a near-duplicate parent pair is not
  preferred over the parents a node truly depends on). Because a tree gives a
  variable a single parent, it cannot represent a variable that depends on two
  otherwise-unrelated predecessors (a v-structure, e.g. `C` a noisy XOR of
  independent `A` and `B`); a degree-2 network takes both parents and recovers it.
  Construction measures the `d` one-way marginals plus one `(parents, node)`
  family joint per non-root node (`2d - 1` marginals) and spends a `select_frac`
  slice on the `d - 1` greedy picks; every slice composes into the same exact
  (\eqn{\epsilon}, \eqn{\delta}) (zCDP for Gaussian, pure \eqn{\epsilon} for
  Laplace). The network is forward-sampled ancestrally, so like the tree it needs
  no PGM inference. `degree` is capped to `d - 1`, shares the
  `bins^(degree + 1)` cell-count warning, is an **alternative** to `structure_frac`
  and `select = "adaptive"`, and is refused on longitudinal / linked releases.
  `degree = 1` (default) leaves the Chow-Liu path byte-identical.
* **Track B: AIM-style budget annealing and higher treewidth for adaptive
  selection.** `dp_control(select = "adaptive")` gains `anneal = TRUE`, a
  data-adaptive round schedule in place of the fixed `d - treewidth` rounds. The
  one-way marginals take a fair fixed share; the clique measurements and their
  exponential-mechanism selections start at a small per-round quantum (large
  noise) and **double** whenever a round's measured signal fails to beat its noise
  floor (AIM's \eqn{\sigma}-halving rule). After the mandatory `d - treewidth`
  spanning cliques (which keep every variable covered, so the sampler stays
  PGM-free) any **surplus** budget is spent on extra rounds that re-measure the
  worst-fit existing clique — inverse-variance combined, so \eqn{\rho} /
  \eqn{\epsilon} adds exactly. The final round absorbs the exact remainder, so the
  release is still exactly (\eqn{\epsilon}, \eqn{\delta}) over a **variable**
  number of rounds. The accounting reports the realised schedule (rounds, spanning
  vs refinement, how often \eqn{\sigma} halved, and the noise range). The
  `treewidth` cap is also lifted: `treewidth = 3` now measures four-way cliques
  (interactions no three-way marginal can see, e.g. a 3-bit parity), with a
  warning when a clique's cell count (`bins^(treewidth + 1)`) grows large enough
  that per-cell noise would dominate. `anneal = FALSE` (default) leaves the
  fixed-schedule adaptive path byte-identical.
* **Track B: AIM-style adaptive marginal selection.** New
  `dp_control(select = "adaptive", treewidth = w)` opts a flat `synth()` DP
  release out of measuring a predetermined marginal set. After the one-way
  marginals it grows the dependency model one marginal at a time, at each round
  using the **exponential mechanism** to privately pick the marginal the
  model-so-far fits worst, then measuring it. The selected cliques form a
  bounded-treewidth junction tree by construction, so the model is a rooted
  junction-tree sampler (no PGM/IPF inference): `treewidth = 1` is a spanning tree
  (equivalent model class to `dependence = "tree"`), `treewidth = 2` measures
  three-way marginals that capture interactions a Chow-Liu tree structurally
  cannot. Selection and measurement each spend budget and compose exactly into the
  same (\eqn{\epsilon}, \eqn{\delta}) — zCDP for Gaussian, pure \eqn{\epsilon} for
  Laplace — with a data-independent round count (`d - treewidth`); `select_frac`
  splits the marginal budget between the two. Unlike `structure_frac` (which picks
  a tree from one free noisy scan), adaptive selection is model-error-guided and,
  at `treewidth >= 2`, escapes the tree's expressiveness ceiling. Flat-table only
  (refused on longitudinal / linked releases); `treewidth` shipped for 1 and 2.
  The default `select = "fixed"` leaves every existing path unchanged.
* **Track B: a linked longitudinal child's transitions can now condition on the
  parent.** New `dp_control(transition_parent = p)` re-injects the synthetic
  parent's (subject-invariant) attributes into a longitudinally-modelled child's
  **transition** tensor at every step, not just its initial state. Until now
  `cross_table = TRUE` cross-conditioned only the child's *initial state*; the
  parent's influence then rode the own-lag chain and decayed. With
  `transition_parent = p`, each time-varying child column additionally conditions
  its next value on the `p` immediate-parent attributes most strongly associated
  with it, so parent → child dependence stays anchored across the whole trajectory.
  The parents are chosen automatically and **budget-neutrally** from the
  parent-by-child joints the cross-conditioned initial state already measures — so
  it adds no histogram and no sensitivity (a transition tuple still lands in one
  cell), and the (\eqn{\epsilon}, \eqn{\delta}) accounting is unchanged. Because it
  reuses those joints, `transition_parent > 0` **requires `cross_table = TRUE`** for
  that child (an error is raised otherwise). It composes with `baseline`,
  `transition_order` and `transition_cross`. Linked-only (a flat longitudinal
  `synth()` has no parent); ignored by non-longitudinal children.
* **Track B: baseline columns and deeper transitions now apply to a linked
  longitudinal child.** The two within-unit transition controls of the flat DP
  Markov engine — `dp_control(baseline = ...)` (subject-invariant columns held
  exactly constant within a unit) and `dp_control(transition_order = k,
  transition_cross = m)` (each time-varying column steps on its own last `k` values
  plus the lag-1 values of its `m` most associated companions) — are now honoured
  by any **longitudinally-modelled linked child** (`synth_linked(..., privacy =
  dp_control(longitudinal = ...))`), applied per table: a baseline name is matched
  against each such child's own columns, and the order/cross settings apply
  uniformly to every longitudinal child. The budget stays exact — a baseline column
  drops its transition histogram (sharpening the rest), cross-conditioning is
  budget-neutral, and a higher order lowers the child's transition sensitivity to
  `path_cap[parent] * (branching_cap - order)` (so the order must be at most one
  less than the child's branching cap, validated per table). All three compose with
  the existing combined cross-table + longitudinal path (a cross-conditioned initial
  state can also hold baseline columns and step deepened transitions). Ignored by
  non-longitudinal linked children, exactly as before.
* **Track B: cross-table conditioning and longitudinal transitions now combine on
  the same linked child.** Previously a longitudinally-modelled child table
  (`dp_control(longitudinal = ...)`) could not also be cross-conditioned on its
  parent — `longitudinal` silently won. Now setting `cross_table = TRUE` alongside
  a longitudinal model on the same table cross-conditions that table's
  **initial-state** model on the synthetic parent (the first row of each unit draws
  from a parent-conditioned Chow-Liu tree), and the within-unit transition chain
  then carries that parent dependence across the whole trajectory. The parent
  shapes where a trajectory starts; the transitions stay parent-free, so the only
  extra cost is the `nC * nP` parent-by-child initial-state joints, measured at the
  same first-row sensitivity as the initial marginals and folded into the same
  exact (\eqn{\epsilon}, \eqn{\delta}) composition. A non-longitudinal child still
  cross-conditions all its rows exactly as before, and setting neither flag is
  unchanged.
* **Track B: higher-order and cross-variable DP transitions.** The flat DP Markov
  longitudinal engine can now condition a variable's next value on more than its
  own single previous value. `dp_control(transition_order = k)` uses the last `k`
  of a variable's **own** values (\eqn{P(v_t \mid v_{t-1}, \dots, v_{t-k})}), and
  `dp_control(transition_cross = m)` additionally conditions on the lag-1 values of
  the `m` most strongly associated **other** variables, selected automatically and
  budget-neutrally from the pairwise marginals the Chow-Liu tree already measures
  (so `transition_cross > 0` requires `dependence = "tree"`). Each variable's whole
  transition model is measured as one conditional tensor; early rows (before
  position `k + 1`) are generated by marginalising that same tensor, which is
  post-processing and costs no budget. Extra conditioning columns do **not** change
  the (\eqn{\epsilon}, \eqn{\delta}) budget — a transition tuple still lands in one
  cell — while a higher order *lowers* the transition sensitivity to `cap - order`
  (a person contributes fewer, deeper tuples), so composition stays exact and the
  order must be at most `max_rows_per_person - 1`. `transition_order = 1`,
  `transition_cross = 0` (defaults) keep the first-order own-lag model unchanged.
  (These also apply to a longitudinal linked child — see the top entry.)
* **Track B: baseline columns held exactly constant.** `dp_control(baseline =
  c(...))` names subject-invariant columns in a longitudinal `synth()` release (a
  `structure` with a nesting index). Those columns are modelled once in the
  initial-state model — so their joint distribution and their correlation with the
  first visit are preserved — and then broadcast unchanged to every row of a
  synthetic unit, instead of being stepped through a transition matrix that would
  let a baseline covariate drift across visits. Because a baseline column carries
  no transition histogram, declaring it also **removes** that histogram from the
  release, so the remaining measurements are sharper at the same exact
  (\eqn{\epsilon}, \eqn{\delta}). Declaring a column baseline is public schema
  knowledge and costs no budget; names that match no modelled column are ignored.
  `baseline = NULL` (default) treats every column as time-varying. (Also applies to
  a longitudinal linked child — see the top entry.)
* **Track B: budget-efficient structure learning.** `dp_control(structure_frac =
  f)` learns the flat `dependence = "tree"` Chow-Liu structure from a cheap
  all-pairs scan that spends only the fraction `f` of the marginal budget, then
  concentrates the remaining `1 - f` on re-measuring **only** the chosen tree's
  \eqn{d - 1} edges (plus the `d` one-way marginals). The previous single-pass
  fitter measured all \eqn{\binom{d}{2}} pairwise marginals at full fidelity and
  discarded all but the tree's edges; because structure *selection* tolerates
  noise far better than the *parameters* do, moving the bulk of the budget onto
  the retained edges sharpens the conditionals — increasingly so as the number of
  variables grows. Both passes compose into the same exact (\eqn{\epsilon},
  \eqn{\delta}) budget (the cheap scan and the concentrated re-measurement are two
  sequential releases whose pure-\eqn{\epsilon} adds / zCDP \eqn{\rho} adds).
  `structure_frac = NULL` (default) keeps the single-pass behaviour; the knob is
  inert for `dependence = "independent"`, for fewer than three variables, and for
  longitudinal / linked DP releases.
* **Track B: linked + longitudinal DP.** `dp_control(longitudinal = TRUE)` — or a
  character vector of child-table names — makes a linked DP release model a child
  table's repeated rows as a within-unit time series instead of exchangeable
  records. For such a table the children-per-parent count model doubles as a
  trajectory-length model, an initial-state model is measured over each parent
  unit's first (temporally earliest) child row, and a first-order Markov
  transition matrix \eqn{P(v_t \mid v_{t-1})} is measured per variable over
  consecutive within-unit rows (ordered by the child's own key index). Their
  person-sensitivities differ — initial marginals at the parent's path cap,
  transitions at `path_cap[parent] * (branching_cap - 1)` — and fold into the same
  exactly-composed (\eqn{\epsilon}, \eqn{\delta}) budget as the rest of the
  release; over-cap units are prefix-truncated in temporal order so consecutive
  pairs stay intact. This recovers within-unit autocorrelation (a patient's values
  trending across visits) that the exchangeable model drops. A longitudinally
  modelled table is not simultaneously `cross_table`-conditioned; `longitudinal =
  FALSE` (default) is unchanged, and the flag is inert for flat / longitudinal
  `synth()` releases (which pick the DP Markov engine from the `structure`
  formula).
* **Track B: cross-table conditioning under DP.** `dp_control(cross_table =
  TRUE)` makes a linked DP release model each child table's variables on the
  **synthetic parent's** attributes, not just link to it. For every child table
  with a modellable immediate parent, parent-by-child joint marginals are
  measured at the child grain — one observation per child row, at the parent's
  value carried down the foreign key — so their person-sensitivity is the child's
  per-entity *path cap*, identical to a child one-way marginal, and they fold
  into the same exactly-composed (\eqn{\epsilon}, \eqn{\delta}) budget. The parent
  variables enter the child's Chow-Liu structure as fixed context nodes (a seeded
  maximum-weight spanning tree), so each child variable's single strongest
  predictor may be a parent variable or another child variable; under
  `dependence = "independent"` each child variable conditions on its single best
  parent variable only. At generation the synthetic parent's already-drawn value
  conditions the child draw, and a three-level hierarchy chains it (a grandchild
  conditions on its already-conditioned parent). This preserves genuine
  cross-table statistical dependence — the association the first linked-DP release
  dropped — at the cost of the extra joints' share of budget. `cross_table =
  FALSE` (default) is unchanged; the flag is inert for flat / longitudinal
  `synth()` releases, which have no parent table.
* **Track B: linked multi-table DP.** `synth_linked()` now accepts a
  `dp_control()`, producing a differentially private release across a whole key
  hierarchy at once. The privacy unit is the **root entity**: adding or removing
  one individual — its root row and all the descendant rows that cascade from it —
  changes the release within a formal (\eqn{\epsilon}, \eqn{\delta}) budget.
  Contribution is bounded hierarchically; `max_rows_per_person` is reused as the
  maximum children kept per parent for each child table (a single integer for all
  child tables, or a named list keyed by table name), the root cap being 1.
  Capping top-down bounds each entity's per-table row count to a *path cap* (the
  product of branching caps from the root). Each table's variable marginals and a
  children-per-parent **count histogram** are measured under one exactly-composed
  budget: the summed L1 (Laplace) and summed squared L2 (Gaussian zCDP)
  sensitivities add over every histogram, each scaled by its per-entity path cap
  (variable marginals) or the parent's path cap (count histograms). Numeric bin
  edges are DP-estimated per table (`domain = "dp"`) at the same per-table
  sensitivity, sharing the `domain_frac` slice. Generation is parent-first and
  copies each synthetic parent's surrogate key down as the foreign key, so
  referential integrity holds by construction (`check_linkage()` confirms it) at
  no privacy cost. First-version limitations: child variables are modelled by
  their own within-table marginals (cross-table statistical conditioning is not
  preserved under DP — only referential integrity), within-table longitudinal
  structure is not modelled for a linked table, and constraints / `unit = "row"`
  are refused. `dp_control(max_rows_per_person = ...)` now also accepts a named
  vector/list of per-table caps.
* **Track B: longitudinal DP (a DP Markov model).** When the `structure`
  declares a nesting index (`synth(data, ~ id / visit, privacy = dp_control(...))`),
  the differentially private engine now preserves within-unit temporal structure
  instead of flattening to a surrogate id. It is the private analogue of Track
  A's initial-state + lag-1 transition model: under one composed budget it
  measures a **length histogram** (rows per person), the **initial-state**
  one-/two-way marginals (the `t = 1` row), and a **transition matrix**
  \eqn{P(v_t \mid v_{t-1})} per variable, then draws a length, an initial row and
  steps the transitions — so autocorrelation across visits survives the noise.
  Because a person contributes a whole trajectory, `max_rows_per_person` must be
  set to the public maximum number of visits (`>= 2`); it caps each person's
  effect on the transition histograms (a length-\eqn{\le c} trajectory has at most
  \eqn{c-1} consecutive pairs). The composition is exact: total L1 sensitivity
  \eqn{= 1 + n_{\mathrm{init}} + |V|\,(c-1)} (Laplace) and summed squared L2
  \eqn{= 1 + n_{\mathrm{init}} + |V|\,(c-1)^2} (Gaussian zCDP), so the reported
  (\eqn{\epsilon}, \eqn{\delta}) still holds end to end and the accounting print
  breaks the budget down (length + initial + transition). DP-estimated bin edges
  (`domain = "dp"`) compose in as before. The synthetic nesting index is
  regenerated as the within-person position. Two first-cut limitations: every
  non-index column is modelled as time-varying (a constant baseline may drift a
  little between a synthetic person's visits) and transitions are per-variable.
  Flat `~ id` DP releases are unchanged.
* **Track B: rigorous, accounted discretisation.** DP bin edges are no longer
  read silently from the data range. `dp_control()` gains a `domain` argument
  (default `"dp"`): numeric variables without a public range in `bounds` now have
  their edges **estimated under differential privacy** (a clamp-free
  exponential-mechanism quantile at each end), with the cost composed into the
  reported (\eqn{\epsilon}, \eqn{\delta}) via a `domain_frac` budget slice
  (default 10%). `domain = "public"` requires `bounds` for every numeric variable
  and spends no budget on the domain; `domain = "data"` keeps the old, warned,
  unaccounted behaviour for benchmarking only. In the rigorous modes a bare
  `character` column is refused (pass categoricals as `factor`/`logical`, whose
  levels are public). The accounting record and its print method now report how
  the edges were chosen and what they cost. This closes the one data-dependent
  step previously excluded from the DP accounting.
* Phase 0 scaffold: package skeleton, MIT license, and CI (`R-CMD-check`).
* Public API with settled signatures and input validation: `synth()`,
  `synth_linked()`, `synth_control()`, `dp_control()`.
* Two-track privacy design documented: high-utility default (Track A) and an
  opt-in person-level differentially private track (Track B).
* Synthetic cardiac example-data generator (`data-raw/make_toy_cardiac.R`).
* **Phase 1 — single-table Track A engine.** `synth()` now returns real
  synthetic data (a `synth_result`) via sequential conditional synthesis:
  `method = "cart"` (leaf conditional bootstrap) and `method = "sample"`. The
  unit identifier is regenerated by subject bootstrap; subject-invariant
  baseline columns are synthesised once per unit and broadcast; `visit_sequence`,
  per-variable methods, and `proper` synthesis are honoured. Added
  `as.data.frame()` / `print()` methods for `synth_result`.
* **Phase 2 — nested / longitudinal support.** The synthetic skeleton is now
  drawn from an explicit learned rows-per-unit count distribution (rather than
  copying whole real units). Time-varying columns are synthesised with an
  initial-state model plus a Markov transition model that conditions on the
  previous row within the unit (lag-1 predictors of every time-varying
  variable), so within-unit autocorrelation across visits is preserved. Units
  of unequal length are supported; the flat and subject-baseline behaviour from
  Phase 1 is unchanged.
* **Phase 3 — linked multi-table engine.** `synth_linked()` now returns real
  synthetic tables (a `synth_linked_result`). The table hierarchy is read from
  the `keys`; root tables are synthesised with the single-table engine and each
  child table is generated from its synthetic parent. Foreign keys are copied
  from the parent so they always resolve (referential integrity); the number of
  child rows per parent is drawn from a learned count distribution that includes
  parents with no children; the child's own index is regenerated; and child
  variables are conditioned on the parent's synthesised attributes (cross-table
  predictors), the own index and earlier child variables. Added `check_linkage()`
  to verify key uniqueness and the absence of orphan child rows, plus
  `as.list()` / `print()` methods for `synth_linked_result`.
* **Phase 4 — methods, constraints and tuning.**
    * **Method registry.** Synthesis methods are now an extensible registry.
      Built-ins: `sample`, `cart`, `forest` (a bagged CART ensemble, no extra
      dependency), `ctree` (conditional-inference trees, via `partykit`), and
      the parametric numeric methods `norm` and `normrank` (rank-preserving).
      `register_method()` adds custom methods by name; `list_methods()` lists
      them. Numeric-only / categorical-only methods are validated against the
      target type.
    * **Constraints / temporal logic.** `rule()` declares a logical constraint,
      either row-wise (`rule(dbp <= sbp)`) or per unit for temporal logic
      (`rule(all(diff(los) >= 0), scope = "unit")`). `synth()` enforces
      constraints by rejection sampling at the unit grain, so nested structure
      is preserved; `constraint_max_tries` bounds the effort.
    * **`synth_control()` wired end to end.** `smoothing` kernel-smooths numeric
      draws (all numeric, or named variables); `predictor_matrix` restricts which
      variables may predict each target; `forest` hyperparameters (`ntree`,
      `mtry`) are honoured.
* **Phase 5 — diagnostics, utility and disclosure risk.**
    * **`diagnose()`.** Utility diagnostics comparing synthetic to real data:
      per-variable marginal fit (Kolmogorov-Smirnov for numeric, total-variation
      distance for categorical), the difference between the real and synthetic
      correlation matrices, and a general propensity-score utility (pMSE) with a
      null-expectation ratio (~1 when the two are indistinguishable). A
      base-graphics `plot()` method overlays each marginal, and both `diagnose()`
      and `disclosure_risk()` dispatch over a named list of tables (or a
      `synth_linked_result`).
    * **`disclosure_risk()`.** Empirical privacy diagnostics: replicated uniques
      (real records unique on the quasi-identifiers and reproduced exactly),
      distance to closest record with a real-to-real baseline, and — given a
      `holdout` of non-training records — a membership-inference check reporting
      an AUC and attacker advantage. Synthetic data is not anonymisation; these
      make residual risk visible.
    * **Vignettes.** Three vignettes (getting started, nested / longitudinal,
      and multi-table linked cardiac), built pandoc-free with `litedown`.
* **Phase 6 — performance and extensibility.**
    * **data.table fast-path.** When the Suggested `data.table` package is
      installed, the synthetic skeleton and the constraint-filtered unit set are
      row-bound with `data.table::rbindlist` instead of repeated
      `do.call(rbind, ...)`, avoiding the quadratic cost of assembling thousands
      of one-unit blocks. Behaviour is identical without `data.table`; only speed
      changes. Subject-invariant column detection is now a single vectorised pass
      rather than a per-column `tapply()` loop.
    * **Parallel replicates.** `synth_control(parallel = ...)` is wired end to
      end: `TRUE` uses all available workers and a positive integer sets an
      explicit count. The `m` independent synthetic replicates are spread over a
      PSOCK cluster in both `synth()` and `synth_linked()`, with independent
      L'Ecuyer RNG streams so a parallel run is reproducible for a fixed
      (`seed`, worker count). Serial synthesis (the default) is unchanged, and
      the engine falls back to serial if a cluster cannot be started.
* **Phase 7 — Track B differential privacy.** `synth(privacy = dp_control(...))`
  now performs differentially private synthesis with a formal
  (\eqn{\epsilon}, \eqn{\delta}) guarantee at **person level**. The engine is a
  marginal-based synthesiser in the PrivBayes / MST lineage: each person's
  contribution is bounded (`max_rows_per_person`), numeric variables are
  discretised to a public grid, low-order marginals are measured under the
  chosen noise mechanism with correct budget composition, and synthetic records
  are drawn from the fitted model.
    * **Mechanisms and accounting.** `mechanism = "laplace"` gives pure
      \eqn{\epsilon}-DP; `"gaussian"` gives approximate DP composed via
      zero-concentrated DP and calibrated to spend exactly (\eqn{\epsilon},
      \eqn{\delta}). Every result carries a `dp_accounting` record (mechanism,
      composed budget, per-cell noise, contribution bound) for the release file.
    * **Dependence structure.** `dependence = "tree"` (default) learns a
      Chow-Liu tree from the same noisy marginals (no extra budget), retaining
      pairwise correlations; `"independent"` keeps only one-way marginals.
    * **Honesty about scope.** DP mode is a flat / single-table release: it does
      not preserve within-unit longitudinal structure, `synth_linked()` DP is
      not available yet, and `rule()` constraints are refused under DP (their
      data-dependent rejection would leak). Bin edges / category sets default to
      the data range with a warning unless supplied publicly via `bounds`; that
      step is excluded from the (\eqn{\epsilon}, \eqn{\delta}) accounting. See
      `vignette("differential-privacy")`.
* **Phase 8 — polish.** Documentation and release preparation, no behaviour
  change. The README now reflects the finished engine (Phases 1-7) with runnable
  examples and CI / lifecycle / license badges; a grouped `pkgdown` reference
  index (`_pkgdown.yml`) and `cran-comments.md` were added; the package-level
  help now points to `vignette("differential-privacy")` instead of internal
  developer notes. Standard `R CMD check` is clean (0 errors / 0 warnings /
  0 notes) with all four vignettes built.
