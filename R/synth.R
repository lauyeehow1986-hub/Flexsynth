#' Synthesise a single (optionally nested / longitudinal) table
#'
#' Generate synthetic data for one table, working natively in long format. The
#' `structure` formula declares the nesting hierarchy (e.g.
#' `~ id / visit / test_number`) so repeated-measures and nested designs are
#' preserved without pivoting to wide format.
#'
#' Passing `privacy = dp_control(...)` switches from the default high-utility
#' engine (Track A) to differentially private synthesis (Track B).
#'
#' @param data A `data.frame` in long format.
#' @param structure A one-sided formula giving the nesting hierarchy, e.g.
#'   `~ id / visit`.
#' @param method Synthesis method (e.g. `"cart"`), applied per variable unless
#'   overridden in `tuning`.
#' @param constraints Optional constraints / rules (temporal, logical, range).
#' @param tuning A [synth_control()] object.
#' @param privacy `NULL` for Track A, or a [dp_control()] object for Track B.
#' @param m Number of synthetic datasets to produce.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Reserved for future use.
#'
#' @return A `synth_result` object (once the engine lands). See the roadmap.
#' @export
#' @examples
#' \dontrun{
#' synth(my_long_df, structure = ~ id / visit, method = "cart", seed = 1)
#' }
synth <- function(data,
                  structure,
                  method = "cart",
                  constraints = NULL,
                  tuning = synth_control(),
                  privacy = NULL,
                  m = 1,
                  seed = NULL,
                  ...) {
  validate_synth_inputs(data, structure, tuning, privacy, m)

  track <- if (is.null(privacy)) "A (high-utility)" else "B (differentially private)"
  stop(sprintf(
    paste0("synth(): the synthesis engine is not implemented yet ",
           "(Phase 1). Inputs validated for Track %s. ",
           "See CLAUDE.md roadmap."),
    track
  ), call. = FALSE)
}
