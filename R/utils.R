# Internal input-validation helpers. Not exported.

validate_control <- function(tuning) {
  if (!inherits(tuning, "synth_control")) {
    stop("`tuning` must be a synth_control() object.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_privacy <- function(privacy) {
  if (!is.null(privacy) && !inherits(privacy, "dp_control")) {
    stop("`privacy` must be NULL or a dp_control() object.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_m <- function(m) {
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 1 || m != as.integer(m)) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_structure <- function(structure) {
  if (!inherits(structure, "formula") || length(structure) != 2L) {
    stop("`structure` must be a one-sided formula, e.g. ~ id / visit.",
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_synth_inputs <- function(data, structure, tuning, privacy, m) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame in long format.", call. = FALSE)
  }
  validate_structure(structure)
  validate_control(tuning)
  validate_privacy(privacy)
  validate_m(m)
  invisible(TRUE)
}
