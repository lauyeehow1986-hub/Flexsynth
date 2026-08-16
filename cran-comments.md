## R CMD check results

Standard `R CMD check` (with vignettes built) passes with **Status: OK** — 0
errors, 0 warnings, 0 notes — on the development platform (Windows 11, R 4.5.2,
full TeX/pandoc/qpdf toolchain) and in CI
(`.github/workflows/R-CMD-check.yaml`, ubuntu/macOS/windows).

Under `R CMD check --as-cran` there are 2 NOTEs, both benign:

1. **New submission.** This is the first submission of `flexsynth` to CRAN.

2. **`checking for non-standard things in the check directory ... NOTE`** —
   lists an empty directory `'NULL'`. This appears only under `--as-cran`, only
   during the examples phase, and only on this local Windows install. It is an
   empty directory that the check driver removes before the run completes. No
   package code creates it: running every example directly (including via
   `devtools::run_examples()` and by sourcing the generated `*-Ex.R` under the
   `--as-cran` example environment) never produces it, and a minimal one-example
   package run under the identical `--as-cran` driver does not produce it either.
   Standard `R CMD check` and CI never show this NOTE. It is flagged here so it
   can be confirmed harmless on CRAN's own builders.

## Test environments

* Windows 11, R 4.5.2 (local; full toolchain — TeX Live, pandoc 3.10, qpdf 12)
* win-builder: R-devel and R-release
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest,
  windows-latest (via r-lib/actions)

## Downstream dependencies

None (new package).
