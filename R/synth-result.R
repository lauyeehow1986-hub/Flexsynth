# Constructor + methods for the object returned by synth().

new_synth_result <- function(syn, m, n, method, structure, visit_sequence,
                             fixed, subject, privacy, seed, call,
                             proper = FALSE) {
  structure(
    list(
      syn = syn,                     # data.frame (m == 1) or list of them
      m = m,
      n = n,
      method = method,
      structure = structure,
      visit_sequence = visit_sequence,
      fixed = fixed,
      subject = subject,
      privacy = privacy,
      proper = isTRUE(proper),       # whether proper (posterior-predictive) synthesis
      seed = seed,
      call = call
    ),
    class = "synth_result"
  )
}

#' @export
print.synth_result <- function(x, ...) {
  cat("<synth_result>\n")
  if (is.null(x$privacy)) {
    cat("  track        : A (high-utility; NOT differentially private)\n")
  } else {
    p <- x$privacy
    cat(sprintf("  track        : B (differentially private; epsilon = %s%s, %s level)\n",
                p$epsilon,
                if (p$delta > 0) paste0(", delta = ", format(p$delta, scientific = TRUE))
                else "", p$unit))
  }
  cat("  datasets (m) :", x$m, "\n")
  nr <- if (x$m == 1L) nrow(x$syn) else nrow(x$syn[[1L]])
  cat("  rows each    :", nr, " (input:", x$n, ")\n")
  cat("  synthesised  :", if (length(x$visit_sequence))
    paste(x$visit_sequence, collapse = ", ") else "(none)", "\n")
  if (length(x$subject)) {
    cat("  unit-level   :", paste(x$subject, collapse = ", "),
        "(once per unit)\n")
  }
  if (length(x$fixed)) {
    cat("  carried      :", paste(x$fixed, collapse = ", "), "\n")
  }
  cat("  method       :", paste(unique(x$method), collapse = ", "), "\n")
  cat("\nGet the data with as.data.frame(x)",
      if (x$m > 1L) "or x$syn[[i]]" else "", "\n")
  invisible(x)
}

#' Extract the synthetic data frame from a synth_result
#'
#' @param x A `synth_result`.
#' @param ... Unused.
#' @return The synthetic `data.frame` (for `m == 1`).
#' @export
as.data.frame.synth_result <- function(x, ...) {
  if (x$m == 1L) {
    x$syn
  } else {
    stop("m > 1: pick a dataset with x$syn[[i]].", call. = FALSE)
  }
}
