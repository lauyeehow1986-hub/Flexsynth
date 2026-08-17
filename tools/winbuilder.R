#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# flexsynth win-builder submission helper
#
# Uploads the current source tree to CRAN's win-builder service, which builds
# and `R CMD check`s it under R-devel and R-release on Windows and emails the
# results to the package maintainer (the Maintainer: address in DESCRIPTION).
# This is the pre-CRAN gate: run it and wait for BOTH result emails (0 errors /
# 0 warnings, only the benign NOTEs listed in cran-comments.md) before tagging a
# release or calling devtools::submit_cran().
#
# This is a network action; run it yourself from an interactive R session:
#
#     source("tools/winbuilder.R")           # from the package root, or
#     Rscript tools/winbuilder.R              # from a shell
#
# Nothing is submitted to CRAN by this script — only to the win-builder check
# service. Results arrive by email, usually within ~30-60 minutes.
# ---------------------------------------------------------------------------

pkg <- "C:/Users/lauye/Downloads/flexsynth"   # adjust if the checkout moved

if (!requireNamespace("devtools", quietly = TRUE))
  stop("Install 'devtools' first: install.packages('devtools')")

## Confirm the version we are about to send, and that the maintainer email is
## set (that is where win-builder mails the results). DESCRIPTION uses
## Authors@R, so there is no literal Maintainer field to read -- derive the
## creator ("cre") from Authors@R instead.
desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
maintainer <- local({
  cols <- colnames(desc)
  if ("Maintainer" %in% cols) return(desc[, "Maintainer"])
  if ("Authors@R" %in% cols) {
    people <- eval(parse(text = desc[, "Authors@R"]))
    cre <- Filter(function(p) "cre" %in% p$role, people)
    if (length(cre))
      return(format(cre[[1L]], include = c("given", "family", "email")))
  }
  "(unknown -- check DESCRIPTION)"
})
message("Package : ", desc[, "Package"])
message("Version : ", desc[, "Version"])
message("Maint.  : ", maintainer)
message("Uploading to win-builder (R-devel and R-release). ",
        "Results will be emailed to the maintainer address above.\n")

## R-devel is the one CRAN cares about most; R-release catches release-only
## regressions. Both build the vignettes (litedown, no pandoc required).
devtools::check_win_devel(pkg, quiet = FALSE)
devtools::check_win_release(pkg, quiet = FALSE)

message("\nSubmitted. Watch the maintainer inbox for two result emails.\n",
        "Expected: Status OK or only the benign NOTEs in cran-comments.md.\n",
        "Only after both come back clean: tag v", desc[, "Version"],
        " and run devtools::submit_cran(\"", pkg, "\").")
