# CLAUDE.md — flexsynth (R package)

## What this is
`flexsynth` is a production-quality **R package for generating high-quality synthetic data**
from real datasets of *any* structure, with first-class native support for **nested,
longitudinal, and multi-table *linked* data**. It is built for realistic clinical / health-
informatics use cases (e.g. admissions → procedures → labs → medications, linked by patient
and admission id) but works for any tabular, panel, or survival data.

**Design north star:** work **natively in long format**. Users must **never** be forced to
pivot nested or longitudinal data to wide format to synthesise it.

**Privacy model — two tracks, default off.** Privacy is a *mode*, not the whole package:

- **Track A — high-utility (default).** The sequential CART/forest engine below. Makes **no
  formal guarantee**; instead ships honest **empirical disclosure-risk diagnostics** (replicated
  uniques, distance-to-closest-record, a simple membership-inference check) so users can judge
  residual risk. This is the workhorse for method development and teaching.
- **Track B — differentially private (opt-in).** `synth(..., privacy = dp_control(...))` runs a
  **DP-capable synthesiser** and returns the (ε, δ) actually spent, giving a *defensible* formal
  guarantee when governance demands one.

**Unit of privacy for Track B is person-level** (ε-DP protects a *whole patient* — all their
admissions/procedures/labs/meds — not a single row). This is the correct unit for clinical
release and is enforced by bounding each patient's contribution (clamping #admissions, #labs,
etc. before budgeting).

**Be honest about the DP cost.** Person-level DP over rich nested multi-table linked data burns
budget fast, so utility drops and full cross-table linkage is *limited* under Track B — it lands
first on flatter / lower-dimensional / marginal releases (marginal-based methods: PrivBayes /
MST / AIM style), then extends. Never let Track A output be *described* as DP: only Track B
carries the guarantee. The CART engine is **not** DP and must never be labelled as such.

## Non-negotiables
- **NEVER commit real patient data or identifiable clinical studies.** Only synthetic or
  public sample datasets under `data/` / `data-raw/`. Real inputs stay off the repo and out of
  examples, tests, vignettes, and fixtures.
- No network calls in the package itself. No telemetry.
- License is **MIT** (see `LICENSE`); every bundled example dataset must be license-clean.

## Tech stack
- **R (>= 4.1)**. Start on **base R / `data.frame`** for a simpler, correct first cut; the
  public API accepts and returns plain `data.frame`. A **`data.table`** fast-path is a later
  optimisation (Phase 6), swapped in behind the same API — do **not** reach for it early.
- Per-variable synthesis methods (Track A): **CART** (`rpart`), **random forest** (`ranger`),
  **conditional inference trees** (`partykit::ctree`), parametric (norm/logreg/polr/poisson),
  and rank/`sample` — plus user-supplied custom methods.
- Differential privacy (Track B): a separate **DP-capable synthesiser** family — marginal-based
  (PrivBayes / MST / AIM style) with a person-level contribution bound and (ε, δ) accounting.
  Reuse the secure-RNG / Laplace-Gaussian noise ideas from the shinyEncrypt DP work. Keep DP
  engines in `Suggests`; Track A must not depend on them.
- Docs: **roxygen2**; tests: **testthat (3e)**; site: **pkgdown**; checks: `R CMD check`.
- Keep the hard dependency surface small; put heavier engines (`ranger`, `partykit`) in
  `Suggests` and degrade gracefully when absent.

## Conventions
- User-facing functions and arguments: **snake_case** (`synth`, `synth_linked`,
  `synth_control`). Internal helpers may be prefixed `.` or kept un-exported.
- One primary responsibility per file in `R/`; group by concern
  (`synth-engine.R`, `linkage.R`, `constraints.R`, `methods-cart.R`, `diagnostics.R`, …).
- Every exported function has roxygen docs with a **runnable** `@examples` block using bundled
  toy data. Reproducibility: everything stochastic takes a `seed`.
- Validate inputs early with **informative, actionable error messages** (name the offending
  table/variable/key). Prefer `cli`/`rlang::abort` style conditions.
- Plotting must avoid cairo/GTK devices (dev workstation is locked down) — use base/ggplot2
  with default devices; never hard-require a system graphics stack.

## Package layout (target)
```
R/                     Package code (engine, linkage, constraints, methods, diagnostics)
man/                   roxygen-generated docs (do not hand-edit)
tests/testthat/        Unit tests (single-table AND multi-table cases)
vignettes/             Long-form guides (see "Vignettes" below)
data/                  Bundled example datasets (.rda) — synthetic/public only
data-raw/              Scripts that build data/ (kept out of the built package)
inst/                  Extra files shipped with the package
DESCRIPTION            Metadata + dependency tiers (Imports/Suggests)
NAMESPACE              roxygen-generated exports (do not hand-edit)
LICENSE                MIT
docs/                  Design notes, architecture, roadmap
```

## Public API (stable surface)
```r
# Single table (formula defines the nesting hierarchy)
# privacy = NULL  -> Track A (high-utility). privacy = dp_control(...) -> Track B (DP).
synth(data, structure = ~ id / visit / test_number, method = "cart",
      constraints = NULL, tuning = synth_control(), privacy = NULL,
      m = 1, seed = NULL, ...)

# Multiple linked nested-longitudinal tables, synthesised jointly
synth_linked(
  tables      = list(admissions = adm_df, procedures = proc_df, labs = lab_df),
  structures  = list(admissions = ~ id / admission_id,
                     procedures = ~ id / admission_id / procedure_number,
                     labs       = ~ id / admission_id / lab_number),
  keys        = list(admissions = c("id","admission_id"),
                     procedures = c("id","admission_id"),
                     labs       = c("id","admission_id")),
  method      = "cart",
  constraints = NULL,
  tuning      = synth_control(),
  m           = 1,
  seed        = NULL, ...
)

# Fine-grained control (strong defaults for beginners, deep knobs for power users)
synth_control(visit_sequence = NULL, predictor_matrix = NULL, method = NULL,
              smoothing = NULL, proper = FALSE, k = NULL, cart = list(...),
              forest = list(...), parallel = FALSE, ...)
```
Plus helpers: `constraint()` / rule builders, `dp_control(epsilon, delta, unit = "person", ...)`
for Track B, `diagnose()` (+ `plot()` methods), `check_linkage()`, `disclosure_risk()`,
model-inspection and result-extraction accessors, and `register_method()` for custom
synthesisers. A DP run reports the (ε, δ) spent on the returned object.

## Architecture (how joint/linked synthesis works)
1. **Order the tables** parent → child from the declared hierarchy/keys (e.g. patients →
   admissions → procedures/labs/meds). Detect cycles and error clearly.
2. **Synthesise top-level units first** (sequential conditional synthesis per variable within
   the table), then synthesise each child table **conditional on the already-synthesised
   parent rows**, drawing per-parent child counts from the learned structural distribution
   (visits per id, procedures per admission, …).
3. **Maintain referential integrity**: every synthetic child row carries a foreign key that
   matches a synthetic parent; no orphans. `check_linkage()` verifies this.
4. **Cross-table dependence**: allow parent-derived variables into a child's predictor matrix
   so learned relationships survive across tables.
5. **Constraints** are applied during/after generation: temporal order
   (`admission_date <= procedure_date <= discharge_date`), logical consistency, range, and
   hierarchical/referential rules — enforced across tables, not just within one.

Keep the engine modular so a new method or diagnostic is a new file + registration, never a
rewrite of the core loop.

## Vignettes (must exist and knit clean)
1. **Getting started** — flat / cross-sectional + simple longitudinal (long format).
2. **Nested hierarchical data** in long format (patients → visits → tests).
3. **Multi-table linked nested longitudinal** — realistic **cardiac** example modelled on the
   actual domain: patients → cardiac admissions → procedures (e.g. angiography / PCI) → labs
   (e.g. troponin, lipids) → medications, showing linkage preservation and diagnostics. Build
   this from **synthetic/simulated** cardiac data that mirrors the structure — never from real
   NHCS / SCDB records.

## Diagnostics & evaluation (what "good" means)
- Univariate + multivariate distributions, within and **across** tables.
- Correlations/associations; multilevel & longitudinal model coefficients (real vs synthetic).
- Structural fidelity: visits-per-id, procedures-per-admission, etc.
- Linkage quality + referential integrity checks.
- Disclosure-risk metrics via `disclosure_risk()` — replicated uniques, distance-to-closest-
  record, a simple membership-inference check — reported honestly (Track A's safety story).
- For Track B: the (ε, δ) actually spent and the per-patient contribution bounds applied.
- Diagnostic plots for each of the above.

## Build & test
```bash
Rscript -e "devtools::document()"    # regenerate man/ + NAMESPACE from roxygen
Rscript -e "devtools::test()"        # run testthat suite
Rscript -e "devtools::check()"       # full R CMD check (aim: 0 errors/warnings/notes)
Rscript -e "pkgdown::build_site()"   # optional docs site
```
Tests must cover **both** single-table and multi-table paths, referential integrity after
synthesis, constraint enforcement, and reproducibility from a fixed `seed`.

## Implementation guidance
Borrow the best ideas — sequential synthesis / rules / CART from **synthpop**, hierarchical
generation from **simstudy**, and relational/multi-table synthesis patterns — but ship a clean,
modern, **long-format-native, multi-table-aware** API. Prioritise, in order: ease of use for
common cases · deep controllability via `tuning` · correct preservation of nested structure
**and** cross-table linkage · strong defaults with full transparency for power users.

## Roadmap (phased)
- **Phase 0** — scaffold: DESCRIPTION/NAMESPACE/LICENSE, package skeleton, toy datasets, CI.
- **Phase 1** — single-table sequential engine (CART) + `synth()` + basic diagnostics.
- **Phase 2** — nested long-format support + structural distributions.
- **Phase 3** — `synth_linked()` joint multi-table engine + referential integrity.
- **Phase 4** — full method set, constraints/temporal logic, `synth_control()` tuning.
- **Phase 5** — diagnostics/utility suite + `disclosure_risk()` (Track A safety) + plots + three
  vignettes.
- **Phase 6** — performance (`data.table` fast-path, parallel), extensibility
  (`register_method()`), polish.
- **Phase 7** — **Track B (differential privacy)**: `dp_control()`, person-level contribution
  bounding, a marginal-based DP synthesiser with (ε, δ) accounting (flat/marginal releases
  first), and a DP vignette. Opt-in; must not regress Track A.

See [docs/roadmap.md](docs/roadmap.md) for detail as it lands.

## Key references
- synthpop (sequential synthesis / CART / rules): https://cran.r-project.org/package=synthpop
- simstudy (hierarchical data generation): https://cran.r-project.org/package=simstudy
- data.table: https://cran.r-project.org/package=data.table
- Repo: https://github.com/lauyeehow1986-hub/Flexsynth
