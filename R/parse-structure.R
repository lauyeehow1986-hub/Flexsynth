# Internal: interpret the `structure` formula against a data frame.
# Not exported.

# A one-sided formula such as ~ id / visit / test_number declares the nesting
# hierarchy. The first variable is the unit identifier (the "person"); any
# further variables are structural / time indices carried alongside it. Every
# named variable must exist in `data`.
parse_structure <- function(structure, data) {
  vars <- all.vars(structure)
  if (length(vars) == 0L) {
    stop("`structure` must name at least one identifier, e.g. ~ id.",
         call. = FALSE)
  }
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop(sprintf(
      "structure variable(s) not found in `data`: %s",
      paste(missing_vars, collapse = ", ")
    ), call. = FALSE)
  }
  list(
    vars   = vars,
    id     = vars[1L],       # unit identifier (regenerated, never modelled)
    nested = vars[-1L]       # structural indices carried from the skeleton
  )
}
