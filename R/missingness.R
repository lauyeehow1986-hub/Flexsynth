# Missing-data model for Track A synthesis.
#
# A variable with NAs is split into (1) a companion *missingness indicator* -- a
# factor with levels "FALSE"/"TRUE" -- synthesised in sequence just before the
# variable, so whether a value is missing can depend on the variables already
# synthesised (and, for a time-varying column, carry its own within-unit
# dependence), and (2) the value itself, whose model is fitted on the observed
# rows only. The variable's NAs are meanwhile replaced by a placeholder so its
# role as a *predictor* of later variables never propagates missingness. After
# generation the indicator drives which synthetic rows are set back to NA and is
# then dropped. This mirrors synthpop's approach and preserves both the
# missingness rate and its association with other variables.

# Placeholder for a variable's missing entries: the median for numeric columns,
# the modal observed value otherwise (kept in the column's own type).
impute_placeholder <- function(x) {
  obs <- x[!is.na(x)]
  if (!length(obs)) return(x[[1L]])                 # all-missing: any typed value
  if (is.numeric(obs)) {
    ph <- stats::median(obs)
    if (is.integer(x)) ph <- as.integer(round(ph))
    return(ph)
  }
  ux <- unique(obs)
  ux[which.max(tabulate(match(obs, ux)))]           # modal value, original type
}

# Augment `data` with a missingness indicator for every `synth_cols` variable
# that has any NA, impute those NAs with a placeholder, and return the widened
# data, the new synthesis order (indicator inserted before each variable), and a
# map value-column -> indicator-column.
prepare_missingness <- function(data, synth_cols) {
  na_map    <- list()
  new_order <- character(0)
  data2     <- data
  for (v in synth_cols) {
    if (anyNA(data[[v]])) {
      ind <- paste0(".na_", v)
      while (ind %in% names(data2)) ind <- paste0(ind, "_")
      data2[[ind]] <- factor(ifelse(is.na(data[[v]]), "TRUE", "FALSE"),
                             levels = c("FALSE", "TRUE"))
      col <- data2[[v]]; col[is.na(col)] <- impute_placeholder(data[[v]])
      data2[[v]] <- col
      na_map[[v]] <- list(ind = ind)
      new_order <- c(new_order, ind, v)
    } else {
      new_order <- c(new_order, v)
    }
  }
  list(data = data2, synth_cols = new_order, na_map = na_map)
}

# Names of the indicator columns in a na_map.
indicator_cols <- function(na_map) {
  if (!length(na_map)) return(character(0))
  vapply(na_map, `[[`, character(1), "ind")
}

# The training rows on which a variable's *value* model should be fitted: all
# rows when the variable has no indicator, otherwise only its observed rows.
fit_rows <- function(train, v, na_map) {
  if (is.null(na_map) || is.null(na_map[[v]])) return(seq_len(nrow(train)))
  obs <- which(train[[na_map[[v]]$ind]] == "FALSE")
  if (length(obs)) obs else seq_len(nrow(train))    # guard: no observed rows
}

# Impose the synthesised missingness on `syn`: set each value column to NA where
# its indicator says "TRUE", then drop the indicator columns.
impose_missingness <- function(syn, na_map) {
  for (v in names(na_map)) {
    ind  <- na_map[[v]]$ind
    miss <- !is.na(syn[[ind]]) & syn[[ind]] == "TRUE"
    syn[[v]][miss] <- NA
    syn[[ind]] <- NULL
  }
  syn
}
