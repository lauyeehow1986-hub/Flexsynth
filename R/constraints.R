# Constraint / temporal-logic system for Track A synthesis.
#
# A constraint is an R expression that must hold in the synthetic data. Rules are
# enforced by rejection sampling at the *unit* grain: a synthetic unit is kept
# only if all of its rows satisfy every rule, and units are regenerated until
# enough valid ones accumulate. Working at the unit grain means row-wise rules
# (e.g. dbp <= sbp) and within-unit temporal rules (e.g. los non-decreasing
# across visits) are both expressible without breaking nested structure.

#' Declare a synthesis constraint
#'
#' Capture a logical rule the synthetic data must satisfy. Rules are enforced by
#' [synth()] via rejection sampling at the unit grain (see *Details*).
#'
#' The expression is written in terms of the data's column names and is **not**
#' evaluated when `rule()` is called — it is captured and checked later against
#' the synthetic data.
#'
#' @param expr A logical expression over the column names, e.g. `dbp <= sbp`.
#'   With `scope = "row"` it is evaluated across a unit's rows and must be `TRUE`
#'   for every row. With `scope = "unit"` it is evaluated once per unit (the
#'   columns are that unit's rows in positional order) and must return `TRUE` —
#'   use this for temporal logic such as `all(diff(los) >= 0)`. Missing (`NA`)
#'   results are treated as satisfied.
#' @param scope `"row"` (default) or `"unit"`.
#' @param label Optional human-readable label; defaults to the deparsed
#'   expression.
#'
#' @return A `flexsynth_rule` object.
#' @export
#' @examples
#' # A row rule and a within-unit temporal rule.
#' rule(dbp <= sbp)
#' rule(all(diff(visit_time) > 0), scope = "unit")
rule <- function(expr, scope = c("row", "unit"), label = NULL) {
  ex <- substitute(expr)
  scope <- match.arg(scope)
  if (is.null(label)) label <- paste(deparse(ex), collapse = " ")
  structure(list(expr = ex, scope = scope, label = label),
            class = "flexsynth_rule")
}

#' @export
print.flexsynth_rule <- function(x, ...) {
  cat(sprintf("<flexsynth_rule [%s]> %s\n", x$scope, x$label))
  invisible(x)
}

# Coerce the user-supplied `constraints` argument into a validated list of rules
# (or NULL). Accepts a single rule or a list of rules; every variable a rule
# mentions must exist in `data`.
normalise_constraints <- function(constraints, data) {
  if (is.null(constraints)) return(NULL)
  if (inherits(constraints, "flexsynth_rule")) constraints <- list(constraints)
  if (is.list(constraints) && !length(constraints)) return(NULL)  # empty == none
  if (!is.list(constraints) ||
      !all(vapply(constraints, inherits, logical(1), "flexsynth_rule"))) {
    stop("`constraints` must be a rule() or a list of rule() objects.",
         call. = FALSE)
  }
  nm <- names(data)
  for (r in constraints) {
    miss <- setdiff(all.vars(r$expr), nm)
    if (length(miss)) {
      stop(sprintf("constraint [%s] references unknown column(s): %s",
                   r$label, paste(miss, collapse = ", ")), call. = FALSE)
    }
  }
  constraints
}

# Is a single synthetic unit (a block of rows in positional order) valid under
# all rules? Row rules must hold for every row; unit rules must return all-TRUE.
# NA results count as satisfied. A rule that errors marks the unit invalid.
unit_ok <- function(block, rules) {
  for (r in rules) {
    res <- tryCatch(eval(r$expr, envir = block, enclos = baseenv()),
                    error = function(e) FALSE)
    res <- as.logical(res)
    if (r$scope == "row") {
      if (any(!is.na(res) & !res)) return(FALSE)
    } else {
      if (!all(res | is.na(res))) return(FALSE)
    }
  }
  TRUE
}

# Enforce `rules` on a generator `gen_one()` (which returns one full synthetic
# table). Valid units are accumulated across regenerations until as many as the
# first unconstrained draw produced are collected, then relabelled 1..N (the
# unit identifier stays a clean sequence) and rebound into one table.
apply_constraints <- function(gen_one, id_col, rules, control) {
  first  <- gen_one()
  target <- length(unique(first[[id_col]]))
  valid  <- Filter(function(b) unit_ok(b, rules), split(first, first[[id_col]]))

  max_tries <- control$constraint_max_tries %||% 50L
  tries <- 1L
  while (length(valid) < target && tries < max_tries) {
    nxt <- gen_one()
    ok  <- Filter(function(b) unit_ok(b, rules), split(nxt, nxt[[id_col]]))
    valid <- c(valid, ok)
    tries <- tries + 1L
  }

  if (length(valid) < target) {
    warning(sprintf(paste0("constraints: kept %d of %d units after %d attempt(s); ",
                           "returning the valid units only. Relax the rules or ",
                           "raise constraint_max_tries."),
                    length(valid), target, tries), call. = FALSE)
    target <- length(valid)
  }
  if (!target) {                                   # nothing satisfied the rules
    out <- first[0, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  valid <- valid[seq_len(target)]
  for (i in seq_along(valid)) valid[[i]][[id_col]] <- i
  out <- do.call(rbind, valid)
  rownames(out) <- NULL
  out
}
