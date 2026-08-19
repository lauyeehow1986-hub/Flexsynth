# flexsynth 0.2.1 release checklist

This checklist is the release record for version 0.2.1. Do not tag, publish, or
submit a source tree that differs from the commit that passed these gates.

## Scope and metadata

- [x] `DESCRIPTION` version is 0.2.1 and its claims match the supported scope.
- [x] The first `NEWS.md` heading is `flexsynth 0.2.1`; post-0.2.0 corrections
      are not attributed to the already-published 0.2.0 release.
- [x] README and vignettes distinguish descriptive diagnostics from privacy or
      release-safety guarantees.
- [x] README examples exclude regenerated identifiers from utility and
      disclosure-risk comparisons.
- [x] CRAN and win-builder helper scripts resolve the package root dynamically.
- [x] `cran-comments.md` records only checks actually run on the prepared source
      tree; re-confirm it against the final commit before submission.

## Verification gates

- [x] `R CMD build` succeeds from a clean staged source tree.
- [x] `R CMD check --as-cran --no-manual flexsynth_0.2.1.tar.gz` completes with
      0 errors, 0 warnings and the expected new-submission NOTE.
- [x] The complete `testthat` suite passes: 965 passes, 0 failures or warnings,
      and 2 recorded skips for unavailable reproducible PSOCK parallelism.
- [x] Vignettes build and rebuild during `R CMD check`; README code was reviewed
      against the public API.
- [ ] PDF-manual generation passes in an environment with `pdflatex`; the local
      R 4.5.1 environment does not provide it, while its HTML manual check passes.
- [ ] GitHub Actions passes all five jobs on the final commit: Ubuntu
      R-release, R-devel and R-oldrel-1; Windows R-release; macOS R-release.
- [ ] win-builder R-devel and R-release results are reviewed and acceptable.

## Statistical and privacy review

- [ ] A reviewer checks the pooling formulas and the revised pMSE, correlation,
      TCAP, membership-inference and key-encoding tests.
- [ ] A reviewer confirms Track A is never represented as anonymised or
      differentially private.
- [ ] Any Track B release has an independent review of the privacy unit, public
      domain assumptions, contribution bounds, mechanism, composed budget, and
      retained `dp_accounting` record.
- [ ] Release examples and bundled data contain no confidential or real patient
      data.

## Publish

- [ ] Merge the release branch after review, then re-run or confirm all checks on
      the resulting `main` commit.
- [ ] Create annotated tag `v0.2.1` on that exact commit and publish a GitHub
      release using the 0.2.1 section of `NEWS.md`.
- [ ] Attach the source tarball built from the tagged commit and verify its
      checksum against the locally checked tarball.
- [ ] If submitting to CRAN, source `tools/submit-cran.R` interactively from the
      tagged package root and confirm the emailed submission link.
- [ ] Record the GitHub release URL, tag commit, tarball checksum, CRAN result,
      and release date below.

## Commands

Run from the package root with the intended R installation on `PATH`:

```powershell
R CMD build .
R CMD check --as-cran --no-manual flexsynth_0.2.1.tar.gz
Rscript -e "testthat::test_local('.')"
```

External checks:

```r
source("tools/winbuilder.R")
source("tools/submit-cran.R")  # only after every pre-flight gate is complete
```

## Release record

- Final commit: pending
- GitHub release: pending
- Source SHA-256: `0E0622C2B0E1765E87F9D8EE9C308459F2ADB099649E4AD988244C42E0D27D47`
- CRAN submission/result: pending
- Release date: pending

## Recovery

Before CRAN acceptance, withdraw a faulty submission and delete an incorrect
GitHub release/tag only after confirming nobody relies on it; then rebuild from
the corrected commit. After public distribution or CRAN acceptance, do not
rewrite the tag. Deprecate or fix forward in a new patch release and document
the affected versions in `NEWS.md`.
