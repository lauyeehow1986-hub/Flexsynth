# flexsynth 0.0.0.9000

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
