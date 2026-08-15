# Build the bundled toy cardiac example datasets (SYNTHETIC — not real data).
#
# Structure mirrors the real domain but is fully simulated:
#   patients (id)
#     -> admissions (id / admission_id)         cardiac admissions
#          -> procedures (id / admission_id / procedure_number)  e.g. angiography, PCI
#          -> labs       (id / admission_id / lab_number)        e.g. troponin, LDL
#          -> meds       (id / admission_id / med_number)        e.g. aspirin, statin
#
# Run from the package root:  Rscript data-raw/make_toy_cardiac.R
# It writes data/*.rda via usethis::use_data() (or save() as a fallback).

set.seed(20260815)

n_patients <- 60

cardiac_patients <- data.frame(
  id     = seq_len(n_patients),
  sex    = sample(c("F", "M"), n_patients, replace = TRUE, prob = c(0.45, 0.55)),
  age    = round(rnorm(n_patients, mean = 64, sd = 11)),
  smoker = sample(c(FALSE, TRUE), n_patients, replace = TRUE, prob = c(0.7, 0.3)),
  stringsAsFactors = FALSE
)

# One row per admission (nested under patient).
adm_list <- lapply(cardiac_patients$id, function(pid) {
  n_adm <- 1 + rpois(1, 0.6)
  base_date <- as.Date("2024-01-01") + sample(0:600, 1)
  admit <- sort(base_date + cumsum(sample(30:400, n_adm, replace = TRUE)))
  data.frame(
    id           = pid,
    admission_id = seq_len(n_adm),
    admit_date   = admit,
    los_days     = 1 + rpois(n_adm, 3),
    stringsAsFactors = FALSE
  )
})
cardiac_admissions <- do.call(rbind, adm_list)
cardiac_admissions$discharge_date <-
  cardiac_admissions$admit_date + cardiac_admissions$los_days

# Helper: expand child rows per admission.
expand_children <- function(adm, n_fun, row_fun) {
  out <- vector("list", nrow(adm))
  for (i in seq_len(nrow(adm))) {
    n <- n_fun()
    if (n > 0) out[[i]] <- row_fun(adm[i, ], n)
  }
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

cardiac_procedures <- expand_children(
  cardiac_admissions,
  n_fun = function() rpois(1, 1.1),
  row_fun = function(a, n) data.frame(
    id               = a$id,
    admission_id     = a$admission_id,
    procedure_number = seq_len(n),
    procedure        = sample(c("angiography", "PCI", "echocardiogram", "stress_test"),
                              n, replace = TRUE),
    procedure_date   = a$admit_date + sample(0:max(0, a$los_days), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
)

cardiac_labs <- expand_children(
  cardiac_admissions,
  n_fun = function() 1 + rpois(1, 2),
  row_fun = function(a, n) data.frame(
    id           = a$id,
    admission_id = a$admission_id,
    lab_number   = seq_len(n),
    analyte      = sample(c("troponin", "LDL", "HDL", "creatinine", "hba1c"),
                          n, replace = TRUE),
    value        = round(abs(rnorm(n, 2, 1.5)), 2),
    lab_date     = a$admit_date + sample(0:max(0, a$los_days), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
)

cardiac_meds <- expand_children(
  cardiac_admissions,
  n_fun = function() 1 + rpois(1, 1.5),
  row_fun = function(a, n) data.frame(
    id           = a$id,
    admission_id = a$admission_id,
    med_number   = seq_len(n),
    drug         = sample(c("aspirin", "statin", "beta_blocker", "ACE_inhibitor",
                            "clopidogrel"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
)

datasets <- c("cardiac_patients", "cardiac_admissions", "cardiac_procedures",
              "cardiac_labs", "cardiac_meds")

if (requireNamespace("usethis", quietly = TRUE)) {
  for (nm in datasets) {
    do.call(usethis::use_data, list(as.name(nm), overwrite = TRUE))
  }
} else {
  dir.create("data", showWarnings = FALSE)
  for (nm in datasets) {
    save(list = nm, file = file.path("data", paste0(nm, ".rda")),
         compress = "bzip2")
  }
}

message("Wrote ", length(datasets), " toy cardiac datasets to data/.")
