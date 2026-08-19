#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# flexsynth CRAN submission helper
#
# Wraps devtools::submit_cran(): it rebuilds the source tarball and uploads it
# to CRAN's submission portal, then CRAN emails the maintainer a link to CONFIRM
# the submission (a second, manual step -- nothing reaches the CRAN queue until
# you click that link).
#
# RUN THIS IN AN INTERACTIVE R SESSION (R console or RStudio), NOT via Rscript:
# submit_cran() asks "Ready to submit? [y/N]" and you must answer it.
#
#     source("tools/submit-cran.R")     # from the package root
#
# Pre-flight (do NOT submit until all are true):
#   * win-builder R-devel AND R-release both returned clean
#     (Status OK or only the benign NOTEs in cran-comments.md);
#   * `R CMD check --as-cran` is clean locally;
#   * NEWS.md, cran-comments.md and the DESCRIPTION Version agree;
#   * the matching v<Version> tag / GitHub release are published.
#
# After you run it: watch the maintainer inbox for the CRAN confirmation email
# and click the confirmation link. Then wait for the CRAN incoming checks.
# ---------------------------------------------------------------------------

pkg <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(pkg, "DESCRIPTION")))
  stop("Run this script from the flexsynth package root.")

if (!requireNamespace("devtools", quietly = TRUE))
  stop("Install 'devtools' first: install.packages('devtools')")

if (!interactive())
  stop("Run this in an interactive R session -- submit_cran() needs to prompt ",
       "you for confirmation. Open R/RStudio and source(\"tools/submit-cran.R\").")

## Show what is about to go up, and surface cran-comments.md so you can eyeball
## it one last time (CRAN reads it).
desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
version <- unname(desc[, "Version"])
news_heading <- grep("^# flexsynth ", readLines(file.path(pkg, "NEWS.md")),
                     value = TRUE)[1L]
if (!identical(news_heading, paste("# flexsynth", version)))
  stop("The first NEWS.md heading does not match DESCRIPTION Version ", version, ".")
cran_comments <- readLines(file.path(pkg, "cran-comments.md"))
if (!any(grepl(paste("version", version), cran_comments, fixed = TRUE)))
  stop("cran-comments.md does not name DESCRIPTION Version ", version, ".")
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
message("Package    : ", desc[, "Package"])
message("Version    : ", desc[, "Version"])
message("Maintainer : ", maintainer)
message("Date       : ", format(Sys.Date()))
message("Expected tag: v", version)
message("\n--- cran-comments.md ---")
cat(cran_comments, sep = "\n")
message("\n------------------------\n")

## This will ask you to confirm, then upload and trigger the CRAN confirmation
## email. Nothing is queued at CRAN until you click the link in that email.
devtools::submit_cran(pkg)

message("\nUploaded to the CRAN submission portal.\n",
        "NEXT: open the CRAN confirmation email sent to the maintainer address ",
        "and click the link to actually submit. Then wait for the incoming ",
        "checks / the CRAN team's reply.")
