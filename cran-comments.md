## Submission

This is the first CRAN submission of `flexsynth` (version 0.1.1). Version 0.1.0
was tagged on GitHub only; 0.1.1 fixes a bug in the optional `ctree` method that
win-builder R-devel surfaced (see NEWS.md), where `partykit` is available and a
previously-skipped test runs.

## R CMD check results

Standard `R CMD check` (with vignettes built) passes with **Status: OK** — 0
errors, 0 warnings, 0 notes — on the development platform (Windows 11, R 4.5.2,
full TeX/pandoc/qpdf toolchain, with `partykit` installed) and in CI
(`.github/workflows/R-CMD-check.yaml`, ubuntu/macOS/windows).

Under `R CMD check --as-cran` any NOTEs are benign:

1. **New submission.** This is the first submission of `flexsynth` to CRAN.

2. **Possibly misspelled words in DESCRIPTION** — "anonymisation". This is the
   intended British-English spelling (recorded in `inst/WORDLIST`).

3. **`checking for non-standard things in the check directory ... NOTE`** (local
   Windows only) — an empty directory `'NULL'` created by the `--as-cran` example
   driver, not by any package code; it is removed before the run completes and
   never appears in standard `R CMD check`, CI, or on win-builder.

## Test environments

* Windows 11, R 4.5.2 (local; full toolchain — TeX Live, pandoc 3.10, qpdf 12;
  partykit installed so the full suite runs)
* win-builder: R-devel and R-release
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest,
  windows-latest (via r-lib/actions)

## Downstream dependencies

None (new package).
