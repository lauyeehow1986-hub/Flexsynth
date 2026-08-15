# Constructor + methods for the object returned by synth_linked().

new_synth_linked_result <- function(syn, m, n, hierarchy, method, privacy,
                                    seed, call) {
  structure(
    list(
      syn = syn,                     # named list (m == 1) or list of them
      m = m,
      n = n,                         # named vector of input row counts
      hierarchy = hierarchy,
      method = method,
      privacy = privacy,
      seed = seed,
      call = call
    ),
    class = "synth_linked_result"
  )
}

#' @export
print.synth_linked_result <- function(x, ...) {
  h <- x$hierarchy
  first <- if (x$m == 1L) x$syn else x$syn[[1L]]
  dp <- inherits(x$privacy, "dp_accounting")
  cat("<synth_linked_result>\n")
  if (dp) {
    cat("  track        : B (differentially private; root-entity level)\n")
    cat(sprintf("  guarantee    : (epsilon = %s%s)-DP\n", x$privacy$epsilon,
                if (x$privacy$delta > 0)
                  paste0(", delta = ", format(x$privacy$delta, scientific = TRUE))
                else ", pure eps"))
  } else {
    cat("  track        : A (high-utility; NOT differentially private)\n")
  }
  cat("  datasets (m) :", x$m, "\n")
  cat("  tables       :", length(h$names), "\n")
  for (t in h$order) {
    p <- h$parent[[t]]
    rel <- if (is.na(p)) "(root)" else paste0("child of ", p)
    cat(sprintf("    - %-14s %6d rows  (input %d)  %s\n",
                t, nrow(first[[t]]), x$n[[t]], rel))
  }
  cat("  method       :", x$method, "\n")
  cat("\nGet the tables with as.list(x)",
      if (x$m > 1L) "or x$syn[[i]]" else "", "\n")
  cat("Verify linkage with check_linkage(x)\n")
  invisible(x)
}

#' Extract the synthetic tables from a synth_linked_result
#'
#' @param x A `synth_linked_result`.
#' @param ... Unused.
#' @return A named list of synthetic `data.frame`s (for `m == 1`).
#' @export
as.list.synth_linked_result <- function(x, ...) {
  if (x$m == 1L) {
    x$syn
  } else {
    stop("m > 1: pick a collection with x$syn[[i]].", call. = FALSE)
  }
}
