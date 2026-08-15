# Referential-integrity check for linked tables.

# Verify one collection of tables against a parent / foreign-key map. Returns a
# per-table report data.frame.
verify_collection <- function(tabs, keys, parent, fk) {
  parts <- lapply(names(tabs), function(t) {
    df <- tabs[[t]]
    k  <- keys[[t]]
    dup <- anyDuplicated(key_string(df, k)) > 0L
    p <- parent[[t]]
    orphan <- NA_integer_
    if (!is.na(p)) {
      ck <- key_string(df, fk[[t]])
      pk <- key_string(tabs[[p]], fk[[t]])
      orphan <- as.integer(sum(!(ck %in% pk)))
    }
    data.frame(
      table = t,
      parent = if (is.na(p)) NA_character_ else p,
      rows = nrow(df),
      duplicate_keys = dup,
      orphan_rows = orphan,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parts)
}

#' Check referential integrity of linked tables
#'
#' Verify that a set of linked tables is internally consistent: every table's
#' key identifies its rows uniquely, and every child row's foreign key resolves
#' to an existing parent row (no orphans). Works on a [synth_linked()] result or
#' on a raw named list of tables (supply `keys` in that case).
#'
#' @param object A `synth_linked_result`, or a named list of `data.frame`s.
#' @param keys When `object` is a raw list of tables, a named list of key
#'   columns (as passed to [synth_linked()]). Ignored for a result object.
#' @param verbose Logical; print a short summary. Defaults to `TRUE`.
#'
#' @return Invisibly, a report `data.frame` (one row per table, or per table and
#'   collection when `m > 1`) with an `ok` attribute that is `TRUE` when there
#'   are no duplicate keys and no orphan rows.
#' @export
#' @examples
#' tabs <- list(
#'   patients   = data.frame(id = 1:3),
#'   admissions = data.frame(id = c(1, 2, 9), admission_id = c(1, 1, 1))
#' )
#' check_linkage(tabs, keys = list(patients = "id",
#'                                 admissions = c("id", "admission_id")))
check_linkage <- function(object, keys = NULL, verbose = TRUE) {
  if (inherits(object, "synth_linked_result")) {
    h <- object$hierarchy
    keys <- h$keys; parent <- h$parent; fk <- h$fk
    collections <- if (object$m == 1L) list(object$syn) else object$syn
  } else if (is.list(object) &&
             length(object) &&
             all(vapply(object, is.data.frame, logical(1)))) {
    if (is.null(keys) || !is.list(keys) || !setequal(names(keys), names(object))) {
      stop("supply `keys` as a named list matching the tables.", call. = FALSE)
    }
    pm <- derive_parent_map(keys)
    parent <- pm$parent; fk <- pm$fk
    collections <- list(object)
  } else {
    stop("`object` must be a synth_linked_result or a named list of data frames.",
         call. = FALSE)
  }

  reps <- lapply(seq_along(collections), function(i) {
    r <- verify_collection(collections[[i]], keys, parent, fk)
    if (length(collections) > 1L) r <- cbind(collection = i, r)
    r
  })
  report <- do.call(rbind, reps)

  ok <- all(!report$duplicate_keys) &&
        all(report$orphan_rows == 0L | is.na(report$orphan_rows))
  attr(report, "ok") <- ok

  if (isTRUE(verbose)) {
    if (ok) {
      message("Linkage OK: keys unique, no orphan child rows.")
    } else {
      bad <- report[report$duplicate_keys |
                    (!is.na(report$orphan_rows) & report$orphan_rows > 0L), ,
                    drop = FALSE]
      message("Linkage problems found:")
      print(bad)
    }
  }
  invisible(report)
}
