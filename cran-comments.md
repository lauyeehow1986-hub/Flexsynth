## Submission

This is the first CRAN submission of `flexsynth` (version 0.2.1). Previous
versions were distributed through GitHub only.

The package provides utility-oriented synthesis for supported flat, nested,
longitudinal and tree-linked multi-table data, plus an opt-in person-level
differentially private engine. Version 0.2.1 includes multiple-synthesis pooling,
missing-data and conditional-count models, categorical-association and TCAP
diagnostics, and corrections to linked holdout handling, key encoding, TCAP
comparison, pMSE calibration, the correlation-matrix norm and custom-estimator
validation. See NEWS.md for details.

## R CMD check results

On Windows 11 with R 4.5.1, the source tarball passes
`R CMD check --as-cran --no-manual` with 0 errors, 0 warnings and 1 NOTE:

1. **New submission.** This is the first submission of `flexsynth` to CRAN.

The same run with PDF-manual checking enabled reaches the manual step with all
earlier checks clean, but the local environment does not have `pdflatex`.
Vignettes rebuild and the HTML manual check passes. PDF-manual confirmation is
therefore an external release gate rather than represented as a completed local
check.

The package test suite reports 965 passes, 0 failures and 2 skips. Both skips
exercise optional reproducible PSOCK parallelism, which is unavailable to the
test process on this local Windows configuration; serial behaviour is covered.

## Test environments

* Windows 11 x64 (build 26200), R 4.5.1 (local source-tarball check)
* GitHub Actions and win-builder results will be added only after the final
  release commit has completed those services

## Downstream dependencies

None (new package).
