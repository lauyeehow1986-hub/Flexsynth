# flexsynth roadmap

Phased delivery. Each phase should leave the package installable and green
(`R CMD check` clean, tests passing). Track A is the default high-utility
engine; Track B is the opt-in differentially private engine.

## Phase 0 — scaffold  *(current)*
- [x] Package skeleton: `DESCRIPTION`, `NAMESPACE`, MIT `LICENSE`.
- [x] Public API stubs with roxygen docs: `synth()`, `synth_linked()`,
      `synth_control()`, `dp_control()` (input validation + settled signatures).
- [x] `testthat` suite for the control constructors + input validation.
- [x] Toy **synthetic cardiac** datasets generator (`data-raw/make_toy_cardiac.R`).
- [x] GitHub Actions `R-CMD-check`.

## Phase 1 — single-table engine (Track A)
- Sequential conditional synthesis with the CART method.
- `synth()` returns a real `synth_result`; basic univariate diagnostics.

## Phase 2 — nested / longitudinal support
- Parse the `structure` formula into a hierarchy; learn structural
  distributions (rows per parent) and synthesise in long format.

## Phase 3 — linked multi-table engine (Track A)
- `synth_linked()`: parent-first joint synthesis, consistent foreign keys,
  referential integrity, cross-table predictors. `check_linkage()`.

## Phase 4 — methods, constraints, tuning
- Full method set (forest, ctree, parametric, rank, custom via
  `register_method()`); constraint/temporal-logic system; `synth_control()`
  wired end to end.

## Phase 5 — diagnostics / utility / safety
- `diagnose()` + plots; `disclosure_risk()` (replicated uniques,
  distance-to-closest-record, membership-inference check).
- Three vignettes: getting started · nested hierarchical · multi-table linked
  cardiac.

## Phase 6 — performance & extensibility
- `data.table` fast-path behind the same API; parallelism; polish.

## Phase 7 — Track B (differential privacy)
- `dp_control()` end to end; **person-level** contribution bounding; a
  marginal-based DP synthesiser (PrivBayes / MST / AIM style) with (ε, δ)
  accounting. Flat / marginal releases first, then extend toward linked data.
- DP vignette. Opt-in; must not regress Track A.
