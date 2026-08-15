#' Jointly synthesise multiple linked nested / longitudinal tables
#'
#' Synthesise several related tables together so that referential integrity and
#' cross-table statistical relationships are preserved. The table hierarchy is
#' read from the `keys`: a table is a *child* of the table whose full key equals
#' its own key with the last column dropped (so `c("id", "admission_id")` is a
#' child of `c("id")`). Root tables are synthesised first with the single-table
#' engine ([synth()]); each child table is then generated from its synthetic
#' parent.
#'
#' For every synthetic parent record the number of child rows is drawn from a
#' learned count distribution (including parents with **no** children), the
#' child's foreign key is copied from the parent so it always resolves
#' (referential integrity), the child's own structural index is regenerated, and
#' the child's variables are synthesised conditionally on the parent's
#' synthesised attributes (**cross-table predictors**), the own index and earlier
#' child variables. Conditioning is on the immediate parent; a grandparent
#' reaches a child only through the parent's synthesised values. Use
#' [check_linkage()] to verify the result. Track A carries **no** formal privacy
#' guarantee.
#'
#' @param tables Named list of input `data.frame`s (one per table).
#' @param structures Named list of one-sided formulas, one per table, giving
#'   each table's nesting hierarchy (e.g. `~ id / admission_id`).
#' @param keys Named list of character vectors giving each table's key columns.
#'   The last column is the table's own index; the leading columns are the
#'   foreign key into its parent.
#' @param method Synthesis method applied per variable unless overridden.
#' @param constraints Optional cross-table constraints / rules. Constraint
#'   enforcement currently applies to single-table [synth()] only; a message is
#'   emitted if supplied here.
#' @param tuning A [synth_control()] object.
#' @param privacy `NULL` for Track A, or a [dp_control()] object for Track B.
#' @param m Number of synthetic dataset collections to produce.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Reserved for future use.
#'
#' @return A `synth_linked_result`. Use [as.list()] (for `m == 1`) or `$syn` to
#'   get the named list of synthetic tables.
#' @export
#' @examples
#' patients <- data.frame(id = 1:20,
#'                        sex = sample(c("F", "M"), 20, TRUE),
#'                        stringsAsFactors = FALSE)
#' adm <- do.call(rbind, lapply(patients$id, function(pid) {
#'   n <- 1 + rpois(1, 0.6)
#'   data.frame(id = pid, admission_id = seq_len(n),
#'              los = 1L + rpois(n, 3))
#' }))
#' res <- synth_linked(
#'   tables     = list(patients = patients, admissions = adm),
#'   structures = list(patients = ~ id, admissions = ~ id / admission_id),
#'   keys       = list(patients = "id", admissions = c("id", "admission_id")),
#'   seed = 1
#' )
#' check_linkage(res)
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

  if (!is.null(privacy)) {
    stop(paste0("synth_linked(): differentially private synthesis (Track B) is ",
                "not implemented yet (Phase 7). See CLAUDE.md roadmap."),
         call. = FALSE)
  }
  if (!is.null(constraints)) {
    message("`constraints` are not enforced for linked synthesis yet; ignoring ",
            "them. (Single-table synth() enforces rule() constraints.)")
  }
  for (t in names(tables)) validate_structure(structures[[t]])

  hierarchy <- link_hierarchy(tables, structures, keys)

  if (!is.null(seed)) set.seed(seed)
  gen_fun <- function() synth_linked_once(tables, hierarchy, method, tuning)
  syns <- run_replicates(m, gen_fun, tuning, seed)
  syn <- if (m == 1L) syns[[1L]] else syns

  new_synth_linked_result(
    syn = syn, m = as.integer(m), n = vapply(tables, nrow, integer(1)),
    hierarchy = hierarchy, method = method, privacy = NULL,
    seed = seed, call = match.call()
  )
}
