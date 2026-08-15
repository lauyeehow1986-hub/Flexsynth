# flexsynth

<!-- badges: start -->
<!-- badges: end -->

**Flexible synthetic data for nested, longitudinal and linked multi-table data.**

`flexsynth` generates high-quality synthetic data from real datasets of *any*
structure, working **natively in long format** — no pivoting nested or
longitudinal data to wide. It has first-class support for **multi-table linked**
data (e.g. patients → admissions → procedures / labs / meds) with referential
integrity and cross-table relationships preserved.

> ⚠️ **Status: early scaffold (Phase 0).** The public API and validation are in
> place; the synthesis engine is not implemented yet. See
> [`docs/roadmap.md`](docs/roadmap.md).

## Two privacy tracks

- **Track A — high-utility (default).** Sequential CART/forest synthesis. No
  formal guarantee; ships honest empirical disclosure-risk diagnostics.
- **Track B — differentially private (opt-in).** `synth(..., privacy =
  dp_control(...))` gives a formal, **person-level** (ε, δ) guarantee for
  governed release. Utility is lower and linkage support is initially limited —
  by design.

Synthetic data is **not** anonymisation, and Track A output must never be
described as differentially private.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lauyeehow1986-hub/Flexsynth")
```

## Usage (target API)

```r
library(flexsynth)

# Single nested/longitudinal table
synth(data, structure = ~ id / visit / test_number, method = "cart", seed = 1)

# Multiple linked tables, synthesised jointly
synth_linked(
  tables     = list(admissions = adm, procedures = proc, labs = lab),
  structures = list(admissions = ~ id / admission_id,
                    procedures = ~ id / admission_id / procedure_number,
                    labs       = ~ id / admission_id / lab_number),
  keys       = list(admissions = c("id", "admission_id"),
                    procedures = c("id", "admission_id"),
                    labs       = c("id", "admission_id")),
  method = "cart", seed = 1
)

# Opt into differential privacy (Track B)
synth(data, structure = ~ id / visit,
      privacy = dp_control(epsilon = 1, delta = 1e-6, unit = "person"))
```

## Data & privacy

Bundled example datasets are **fully synthetic** cardiac data (see
`data-raw/make_toy_cardiac.R`). No real patient data ships with this package,
and none should ever be committed.

## License

MIT © 2026 flexsynth authors
