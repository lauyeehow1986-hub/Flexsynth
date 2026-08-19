# flexsynth

<!-- badges: start -->
[![R-CMD-check](https://github.com/lauyeehow1986-hub/Flexsynth/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lauyeehow1986-hub/Flexsynth/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Flexible synthetic data for nested, longitudinal and linked multi-table data.**

`flexsynth` generates utility-oriented synthetic data for supported flat,
nested, longitudinal, and tree-linked structures, working **natively in long
format** — no pivoting nested or longitudinal data to wide. It has first-class
support for **tree-linked multi-table** data (e.g. patients → admissions →
procedures / labs / meds), with referential integrity preserved by construction.

It ships **two engines behind one interface**:

- **Track A — high-utility (default).** A sequential conditional-synthesis
  engine (CART, forest, conditional-inference trees, parametric methods) that
  models within-unit longitudinal dependence and cross-table relationships. No
  formal privacy guarantee; ships empirical utility and disclosure-risk
  diagnostics to inform review of a proposed release.
- **Track B — differentially private (opt-in).** `synth(..., privacy =
  dp_control(...))` implements a **person-level** (ε, δ) mechanism and emits its
  budget-accounting record. A governed release still requires independent review
  of the privacy unit, public-domain assumptions, contribution caps, and selected
  mechanism.

> Synthetic data is **not** anonymisation, and Track A output must never be
> described as differentially private.

## Features

- **Supported structures, native long format.** Flat, nested / longitudinal
  (repeated visits), and tree-linked multi-table data — declared with a compact
  `structure` formula (`~ id / visit / test`).
- **Longitudinal dependence.** Learned rows-per-unit count model, initial-state
  model, and a lag-1 Markov transition model, so within-unit autocorrelation
  across visits is preserved. Subject-invariant baseline columns are synthesised
  once per unit and broadcast.
- **Linked multi-table.** Hierarchy inferred from keys; parent-first generation
  with foreign keys copied from the synthetic parent (referential integrity by
  construction), zero-inflated children-per-parent counts, and cross-table
  predictors. `check_linkage()` verifies key uniqueness and the absence of
  orphans.
- **Extensible methods.** `sample`, `cart`, `forest`, `ctree`, `norm`,
  `normrank` built in; add your own with `register_method()`.
- **Constraints / temporal logic.** `rule()` enforces row-wise or per-unit
  constraints (e.g. `dbp <= sbp`, monotone length-of-stay) by unit-grain
  rejection sampling.
- **Pooled inference.** `pool_synth()` / `synth_glm()` implement the published
  fully-synthetic variance rules (synthpop / Reiter) across `m` synthetic
  datasets. They use large-sample normal intervals; calibration still depends on
  the estimand, synthesis model, sample size, and analysis assumptions.
  `compare_estimates()` scores real-vs-synthetic analyses by confidence-interval
  overlap.
- **Diagnostics.** `diagnose()` (marginal fit, correlation-matrix difference,
  categorical association via Cramer's V, propensity pMSE) and
  `disclosure_risk()` (replicated uniques, distance-to-closest-record,
  membership-inference AUC, and TCAP attribute-disclosure via `target =`).
- **Performance.** Optional `data.table` fast-path and parallel replicates
  (`synth_control(parallel = TRUE)`) with reproducible L'Ecuyer streams.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lauyeehow1986-hub/Flexsynth")
```

The package installs with base-R dependencies. The default `method = "cart"`
requires the suggested `rpart` package; `ranger` / `partykit` unlock additional
tree methods, and `data.table` unlocks the row-binding fast-path. Use
`method = "sample"` when only base R is available.

## Quick start

### Single nested / longitudinal table

```r
library(flexsynth)

df <- data.frame(
  id    = rep(1:20, each = 2),
  visit = rep(1:2, times = 20),
  age   = rep(round(rnorm(20, 60, 8)), each = 2),
  sbp   = round(rnorm(40, 130, 15))
)

res <- synth(df, structure = ~ id / visit, method = "cart", seed = 1)
head(as.data.frame(res))
```

### Multiple linked tables, synthesised jointly

```r
patients <- data.frame(id = 1:50, sex = sample(c("F", "M"), 50, TRUE))
adm <- do.call(rbind, lapply(patients$id, function(pid) {
  n <- 1 + rpois(1, 0.6)
  data.frame(id = pid, admission_id = seq_len(n), los = 1L + rpois(n, 3))
}))

res <- synth_linked(
  tables     = list(patients = patients, admissions = adm),
  structures = list(patients   = ~ id,
                    admissions = ~ id / admission_id),
  keys       = list(patients   = "id",
                    admissions = c("id", "admission_id")),
  seed = 1
)
syn <- as.list(res)
check_linkage(syn, keys = list(patients = "id",
                               admissions = c("id", "admission_id")))
```

### Diagnostics and disclosure risk

```r
syn <- as.data.frame(synth(df, structure = ~ id / visit, seed = 1))
analysis_vars <- c("age", "sbp")             # exclude generated structure keys
d <- diagnose(real = df, syn = syn, vars = analysis_vars)
plot(d)                                  # overlaid marginals
disclosure_risk(real = df, syn = syn, quasi = analysis_vars)
```

### Pooled inference from synthetic data

```r
# Analyse all m synthetic sets with a published fully-synthetic pooling rule.
# A single set analysed naively generally under-states synthesis uncertainty.
res <- synth(df, structure = ~ id / visit, m = 10, seed = 1)
synth_glm(res, sbp ~ age)                 # pooled linear model
# any estimator works via pool_synth(res, function(d) <fit returning coef/vcov>)
```

### Opt into differential privacy (Track B)

```r
dp <- dp_control(epsilon = 1, delta = 1e-6, mechanism = "gaussian",
                 bounds = list(visit = c(1, 2), age = c(18, 100), sbp = c(60, 240)))
dp_res <- synth(df, structure = ~ id, privacy = dp, seed = 1)
dp_res$privacy                           # the (ε, δ) accounting record
```

Track B supports three release shapes:

- a flat table (`~ id`);
- a longitudinal table (`~ id / visit`) using a bounded DP Markov model; and
- a tree-linked hierarchy via `synth_linked()`, using the root entity as the
  privacy unit.

`dp_control()` also exposes privately learned domains, cross-table
conditioning, higher-order transitions, adaptive marginal selection,
Private-PGM reconciliation, and AIM-style models. These controls trade
statistical fidelity, cell sparsity, runtime, and privacy budget; they are not
universally beneficial.

A governed DP release requires more than setting `epsilon`: define the privacy
unit and public domain assumptions, justify contribution caps, retain the
accounting record, review utility at the intended analysis grain, and obtain
independent privacy review. See `vignette("differential-privacy")` for the
supported combinations, accounting model, and limitations.

## Learn more

- `vignette("getting-started")` — getting started
- `vignette("valid-inference")` — pooled inference and attribute-disclosure (TCAP)
- `vignette("nested-longitudinal")` — repeated-measures data
- `vignette("linked-cardiac")` — multi-table linked synthesis
- `vignette("differential-privacy")` — Track B
- [`docs/roadmap.md`](https://github.com/lauyeehow1986-hub/Flexsynth/blob/main/docs/roadmap.md) — phased delivery and what's next

## Data & privacy

Bundled example datasets are **fully synthetic** cardiac data (see
`data-raw/make_toy_cardiac.R`). No real patient data ships with this package,
and none should ever be committed.

## License

MIT © 2026 flexsynth authors
