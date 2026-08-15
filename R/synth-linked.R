#' Jointly synthesise multiple linked nested / longitudinal tables
#'
#' Synthesise several related tables together so that referential integrity and
#' cross-table statistical relationships are preserved. Parent tables (e.g.
#' patients, admissions) are synthesised first; child tables (procedures, labs,
#' meds) are then generated conditionally on the synthetic parent records, with
#' consistent foreign keys.
#'
#' @param tables Named list of input `data.frame`s (one per table).
#' @param structures Named list of one-sided formulas, one per table, giving
#'   each table's nesting hierarchy.
#' @param keys Named list of character vectors giving each table's key columns
#'   (including foreign keys used for linkage).
#' @param method Synthesis method applied per variable unless overridden.
#' @param constraints Optional cross-table constraints / rules.
#' @param tuning A [synth_control()] object.
#' @param privacy `NULL` for Track A, or a [dp_control()] object for Track B.
#' @param m Number of synthetic dataset collections to produce.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Reserved for future use.
#'
#' @return A `synth_linked_result` object (once the engine lands). See the roadmap.
#' @export
#' @examples
#' \dontrun{
#' synth_linked(
#'   tables     = list(admissions = adm, procedures = proc, labs = lab),
#'   structures = list(admissions = ~ id / admission_id,
#'                     procedures = ~ id / admission_id / procedure_number,
#'                     labs       = ~ id / admission_id / lab_number),
#'   keys       = list(admissions = c("id", "admission_id"),
#'                     procedures = c("id", "admission_id"),
#'                     labs       = c("id", "admission_id")),
#'   method = "cart", seed = 1
#' )
#' }
synth_linked <- function(tables,
                         structures,
                         keys,
                         method = "cart",
                         constraints = NULL,
                         tuning = synth_control(),
                         privacy = NULL,
                         m = 1,
                         seed = NULL,
                         ...) {
  if (!is.list(tables) || is.null(names(tables)) || any(names(tables) == "")) {
    stop("`tables` must be a *named* list of data frames.", call. = FALSE)
  }
  if (!all(vapply(tables, is.data.frame, logical(1)))) {
    stop("every element of `tables` must be a data.frame.", call. = FALSE)
  }
  for (nm in c("structures", "keys")) {
    obj <- get(nm)
    if (!is.list(obj) || !setequal(names(obj), names(tables))) {
      stop(sprintf("`%s` must be a named list matching names(tables).", nm),
           call. = FALSE)
    }
  }
  validate_control(tuning)
  validate_privacy(privacy)
  validate_m(m)

  track <- if (is.null(privacy)) "A (high-utility)" else "B (differentially private)"
  stop(sprintf(
    paste0("synth_linked(): the joint synthesis engine is not implemented yet ",
           "(Phase 3). Inputs validated for %d tables, Track %s. ",
           "See CLAUDE.md roadmap."),
    length(tables), track
  ), call. = FALSE)
}
