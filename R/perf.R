# Internal: performance helpers for the hot paths.
#
# `rbind_rows()` is an optional data.table fast-path: it uses
# `data.table::rbindlist` when the (Suggested) data.table package is installed
# and falls back to `do.call(rbind, ...)` otherwise, so behaviour never depends
# on data.table being present — only speed does. `constant_within()` is a
# vectorised base-R replacement for a per-column `tapply()` loop. Both avoid the
# quadratic / per-group-closure costs that dominate on large longitudinal data.

# Is the data.table fast-path available in this session?
have_data_table <- function() {
  requireNamespace("data.table", quietly = TRUE)
}

# Row-bind a list of data.frames that share the same columns, returning a plain
# data.frame. The synthetic skeleton and the constraint-filtered unit set are
# each assembled from thousands of one-unit blocks, where `do.call(rbind, ...)`
# copies the growing result on every call. `data.table::rbindlist` binds in a
# single pass. NULL / empty parts are dropped; column names and types come from
# the parts (which are constructed to match), so the result equals the base
# `do.call(rbind, ...)` value.
rbind_rows <- function(parts) {
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(NULL)
  if (length(parts) == 1L) return(parts[[1L]])
  if (have_data_table()) {
    out <- data.table::rbindlist(parts, use.names = TRUE, fill = FALSE)
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  do.call(rbind, parts)
}

# Which of `cols` are constant within every level of `gid`? Returns a named
# logical vector over `cols`. Rather than a per-column `tapply()` closure (one
# grouped split per column), this sorts the rows by group once and then, for each
# column, asks in a single vectorised pass whether any two adjacent same-group
# rows differ. NA is treated as an ordinary value (NA vs NA is "same", NA vs a
# value is "different"), matching `length(unique(v)) == 1L`, so the result equals
# the earlier `tapply`-based computation.
constant_within <- function(data, gid, cols) {
  if (!length(cols)) return(stats::setNames(logical(0), character(0)))
  n <- length(gid)
  if (n <= 1L) return(stats::setNames(rep(TRUE, length(cols)), cols))
  ord <- order(gid)
  g <- gid[ord]
  same_group <- g[-1L] == g[-n]          # adjacent rows in the same unit
  vapply(cols, function(col) {
    v <- data[[col]][ord]
    a <- v[-1L]; b <- v[-n]
    na_a <- is.na(a); na_b <- is.na(b)
    differ <- (na_a != na_b) | (!na_a & !na_b & a != b)
    !any(same_group & differ)            # constant iff no within-unit change
  }, logical(1))
}
