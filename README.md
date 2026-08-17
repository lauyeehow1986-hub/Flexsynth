# flexsynth

<!-- badges: start -->
[![R-CMD-check](https://github.com/lauyeehow1986-hub/Flexsynth/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lauyeehow1986-hub/Flexsynth/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Flexible synthetic data for nested, longitudinal and linked multi-table data.**

`flexsynth` generates high-quality synthetic data from real datasets of *any*
structure, working **natively in long format** — no pivoting nested or
longitudinal data to wide. It has first-class support for **multi-table linked**
data (e.g. patients → admissions → procedures / labs / meds) with referential
integrity and cross-table relationships preserved.

It ships **two engines behind one interface**:

- **Track A — high-utility (default).** A sequential conditional-synthesis
  engine (CART, forest, conditional-inference trees, parametric methods) that
  models within-unit longitudinal dependence and cross-table relationships. No
  formal privacy guarantee; ships honest empirical disclosure-risk diagnostics
  so residual risk can be judged.
- **Track B — differentially private (opt-in).** `synth(..., privacy =
  dp_control(...))` gives a formal, **person-level** (ε, δ) guarantee for
  governed release, with exact budget accounting.

> Synthetic data is **not** anonymisation, and Track A output must never be
> described as differentially private.

## Features

- **Any structure, long format.** Flat, nested / longitudinal (repeated visits),
  and linked multi-table data — declared with a compact `structure` formula
  (`~ id / visit / test`).
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
- **Valid inference.** `pool_synth()` / `synth_glm()` combine an analysis across
  the `m` synthetic datasets with fully-synthetic variance rules (synthpop /
  Reiter), so confidence intervals and tests are correct rather than
  over-optimistic.
- **Diagnostics.** `diagnose()` (marginal fit, correlation-matrix difference,
  propensity pMSE) and `disclosure_risk()` (replicated uniques,
  distance-to-closest-record, membership-inference AUC, and TCAP
  attribute-disclosure via `target =`).
- **Performance.** Optional `data.table` fast-path and parallel replicates
  (`synth_control(parallel = TRUE)`) with reproducible L'Ecuyer streams.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lauyeehow1986-hub/Flexsynth")
```

The core engine depends only on base R. `rpart` / `ranger` / `partykit` unlock
the tree methods; `data.table` unlocks the row-binding fast-path; all are
Suggested.

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
d <- diagnose(real = df, syn = syn)
plot(d)                                  # overlaid marginals
disclosure_risk(real = df, syn = syn)
```

### Valid inference from synthetic data

```r
# Analyse all m synthetic sets and pool, so the standard errors account for
# synthesis (a single set analysed naively under-states them).
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

DP works for a flat table (`~ id`), a longitudinal one (`~ id / visit`, a DP
Markov model over visits), and a whole linked hierarchy (`synth_linked(...,
privacy = dp_control(...))`, at the root-entity grain) — all with exact
Laplace / zCDP-Gaussian accounting. Under the default `domain = "dp"`, numeric bin
edges are estimated privately and a bare `character` column's category set is
discovered by **DP set-union** (a stability histogram; rare categories fold into an
`"(other)"` catch-all), both funded from `domain_frac` and composed into the exact
budget — so a flat `character` column needs no pre-conversion to a `factor` (this
needs `delta > 0`; a pure-ε release still requires public levels). For a linked release,
`dp_control(cross_table = TRUE)` additionally conditions each child table's
variables on the synthetic parent's attributes (measured parent-by-child joints,
composed into the same budget), so cross-table dependence survives — not just the
key link; and `dp_control(longitudinal = TRUE)` models a child table's repeated
rows as a within-unit DP Markov trajectory (initial-state + per-variable
transitions), so within-unit autocorrelation across visits survives too. Setting
both on the same child **combines** them: its initial-state model is
cross-conditioned on the parent (the only extra cost is the parent-by-child
initial joints) and the transition chain carries that parent dependence across the
trajectory. By default that parent dependence rides the own-lag chain and decays;
`dp_control(transition_parent = p)` instead re-injects the synthetic parent's
subject-invariant attributes into the child's **transition** at every step (each
time-varying column conditions on its `p` most associated parent attributes),
keeping parent→child dependence anchored across the whole trajectory. The parents
are chosen budget-neutrally from the parent-by-child joints the cross-conditioned
initial state already measures, so it costs nothing in the budget and **requires**
`cross_table = TRUE` for that child. In a
longitudinal release, `dp_control(baseline = c(...))` names
subject-invariant columns (e.g. birth sex, a baseline measurement): they are
modelled once in the initial state and held **exactly constant** across a
unit's visits instead of drifting through a transition matrix — which also
drops their transition histograms, sharpening the rest at the same budget. That
same release can deepen its Markov model with
`dp_control(transition_order = k, transition_cross = m)`: each variable's next
value then conditions on its own last `k` values and on the lag-1 values of its
`m` most strongly associated companions (chosen budget-neutrally from the tree's
pairwise marginals) — extra conditioning is free in the (ε, δ) budget, and a
higher order actually lowers the transition sensitivity to `cap - order`. These
two controls apply both to a longitudinal `synth()` release and, per table, to
every longitudinally-modelled linked child. For a
flat tree release, `dp_control(structure_frac = f)` learns the Chow-Liu structure
from a cheap all-pairs scan (a fraction `f` of the budget) and concentrates the
rest on re-measuring only the chosen edges — sharper conditionals at the same
exact budget, more so the more variables there are. Going further,
`dp_control(select = "adaptive", treewidth = w)` replaces the fixed marginal set
with an **AIM-style** selector: it grows the model one marginal at a time, each
round using the exponential mechanism to privately pick the marginal the
model-so-far fits worst, then measuring it (both spends compose into the same
exact budget). At `treewidth = 2` the selected cliques are triangles, so it
captures three-way interactions a tree structurally cannot — the payoff a
tree-only method, `structure_frac` included, cannot reach at any budget
(`treewidth = 3` reaches four-way cliques, e.g. a 3-bit parity no three-way
marginal can see). Adding `anneal = TRUE` turns the fixed `d - treewidth` schedule
into an **AIM-style data-adaptive** one: rounds start noisy and the per-round
budget **doubles** whenever a measurement fails its noise floor, then any surplus
budget re-measures the worst-fit clique — all still composing to the exact same
(ε, δ) over a variable number of rounds. Separately,
`dp_control(dependence = "tree", degree = k)` generalises the Chow-Liu tree (a
degree-1 Bayesian network) to **PrivBayes' GreedyBayes**: each variable may
condition on up to `k` already-generated variables, its parent set chosen greedily
by the exponential mechanism. Because a tree gives a variable one parent, it
cannot represent a variable that depends on two otherwise-unrelated predecessors (a
v-structure, e.g. `C` a noisy XOR of independent `A` and `B`); a degree-2 network
takes both parents and recovers it — all composing to the same exact (ε, δ), and
still forward-sampled with no PGM inference. Finally,
`dp_control(estimator = "pgm")` adds the **Private-PGM inference step** the
PGM-free samplers omit: instead of using each measured marginal locally (each tree
edge as its own `P(child | parent)`), it reconciles the *whole* measured set — one-
ways and edges together — into the single graphical-model distribution that best
fits all of it at once, by belief propagation on the junction tree plus entropic
mirror descent (McKenna et al.'s MST). Reconciliation is pure post-processing of the
already-noised marginals, so it costs **no extra budget** — the (ε, δ) is identical
to the `"local"` release — yet it denoises overlapping marginals into mutual
consistency and lets the otherwise-discarded one-ways constrain the model. It works
for the flat tree and for `select = "adaptive"`. Going all the way,
`dp_control(select = "aim")` is **Full AIM**: it lifts the adaptive selector's
running-intersection constraint, so the exponential mechanism may pick *loopy*
marginals — a pair between two variables already in the model, the cycle a
junction-tree selector structurally cannot close and no tree captures at any budget.
A loopy set has no forward sampler, so AIM always reconciles the whole measured set
into one graphical model over a **triangulated** junction tree with Private-PGM and
samples from that; the `treewidth` cap bounds the triangulated clique size (so
inference stays exact), and at `treewidth = 1` it reduces to an adaptively-selected
tree. Selection and measurement compose to the same exact (ε, δ), and the
reconciliation is budget-neutral. Adding `anneal = TRUE` to `select = "aim"` runs it
on the same data-adaptive σ-halving schedule as the adaptive selector — a baseline of
treewidth-capped new loopy pairs, then surplus budget re-measuring the worst-fit
measured pair, all still triangulated and reconciled with Private-PGM and composing to
the exact same (ε, δ) over a variable round count. By default, `select = "aim"` now
scores each candidate marginal against the current reconciled model's own marginal over
that pair (AIM's actual quality function) rather than the one-way independence product,
so a loopy pair the model already explains (through the marginals measured so far) stops
looking surprising and the budget is steered to the genuinely worst-fit interaction. The
model reference is read from the already-privatised marginals — reconciled each round and
projected onto the candidate, which for a not-yet-measured pair crosses cliques of the
junction tree — so it is pure post-processing: the exponential mechanism's sensitivity
and the exact (ε, δ) are identical either way; only which marginals get selected changes.
The cheaper one-way-product reference remains available as
`dp_control(select = "aim", scoring = "independence")`. See
`vignette("differential-privacy")` for scope and the honest utility trade-off.

## Learn more

- `vignette("getting-started")` — getting started
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
