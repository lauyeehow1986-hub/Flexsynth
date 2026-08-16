# Internal: build a discrete domain for DP marginal synthesis, and enforce the
# per-person contribution bound. Not exported.
#
# DP marginals live over a discrete grid. Numeric variables are cut into
# equal-width bins; categorical variables use their level set. Each column's
# domain object knows how to encode real values to integer cell indices and how
# to decode synthetic indices back to typed values.

# Resolve the single-table (flat / longitudinal) per-person row cap from a
# dp_control. A named per-table cap is only meaningful for linked DP, so reject
# it here with a clear message rather than letting a vector leak into `cap`.
dp_scalar_cap <- function(dp) {
  mrp <- dp$max_rows_per_person
  if (is.null(mrp)) return(1L)
  if (!is.null(names(mrp)) || length(mrp) != 1L) {
    stop(paste0("A named per-table `max_rows_per_person` is only for linked DP ",
                "(synth_linked()); pass a single positive integer to synth()."),
         call. = FALSE)
  }
  as.integer(mrp)
}

# Subsample each person's rows down to at most `cap`, so one person changes any
# marginal by at most `cap`. Returns the capped data and the number of rows
# dropped. Rows are chosen at random within each over-cap person.
dp_contribution_bound <- function(data, id, cap) {
  gid <- data[[id]]
  idx <- seq_len(nrow(data))
  groups <- split(idx, gid)
  over <- lengths(groups) > cap
  if (!any(over)) return(list(data = data, dropped = 0L))
  keep <- unlist(lapply(groups, function(rows) {
    if (length(rows) <= cap) rows
    else rows[sort(sample.int(length(rows), cap))]
  }), use.names = FALSE)
  keep <- sort(keep)
  list(data = data[keep, , drop = FALSE], dropped = nrow(data) - length(keep))
}

# Build the domain object for one column. `levels` supplies a DP-discovered
# category set for a character column (from dp_discover_categories); an "(other)"
# catch-all cell is appended so undiscovered categories have a home.
dp_domain_column <- function(v, nbin, bound = NULL, levels = NULL) {
  if (is.numeric(v)) {
    is_int <- is.integer(v) || all(v == round(v), na.rm = TRUE)
    if (!is.null(bound)) {
      lo <- bound[1L]; hi <- bound[2L]; derived <- FALSE
    } else {
      rng <- range(v, na.rm = TRUE)
      lo <- rng[1L]; hi <- rng[2L]; derived <- TRUE
    }
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
      # Degenerate / constant column: a single cell holding the value.
      lo <- if (is.finite(lo)) lo else 0
      hi <- lo + 1
      edges <- c(lo, hi)
    } else {
      edges <- seq(lo, hi, length.out = nbin + 1L)
    }
    list(kind = "numeric", is_integer = is_int, lo = lo, hi = hi,
         edges = edges, nbin = length(edges) - 1L, derived = derived)
  } else if (is.factor(v)) {
    lv <- levels(v)
    list(kind = "factor", levels = lv, nbin = length(lv), derived = FALSE)
  } else if (is.logical(v)) {
    list(kind = "logical", levels = c("FALSE", "TRUE"), nbin = 2L,
         derived = FALSE)
  } else if (!is.null(levels)) {
    # DP-discovered character domain: the private level set plus an "(other)"
    # catch-all for undiscovered categories. Public (accounted), so not derived.
    lv <- c(as.character(levels), "(other)")
    list(kind = "character", levels = lv, nbin = length(lv), derived = FALSE,
         other_cell = length(lv))
  } else {
    lv <- sort(unique(as.character(v)))
    list(kind = "character", levels = lv, nbin = length(lv), derived = TRUE)
  }
}

# Estimate one quantile of `x` under the exponential mechanism, clamp-free: the
# output is one of the observed values, chosen with probability proportional to
# exp(eps_q * u / (2 * cap)) where the utility u = -|rank - target| favours the
# value whose rank is closest to `prob * n`. This is eps_q-DP at person
# sensitivity `cap` (adding / removing one person's <= cap rows shifts any rank
# by at most cap). Using an inner quantile (small `prob`, or 1 - small) rather
# than the exact extreme keeps a lone outlier from dominating the choice; values
# beyond the returned edge are clamped into the end bin by dp_encode().
dp_exp_quantile <- function(x, prob, eps_q, cap) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(NA_real_)
  if (n == 1L) return(x[1L])
  xs <- sort(x)
  target <- prob * n
  u <- -abs(seq_len(n) - target)              # rank utility per sorted value
  logw <- (eps_q * u) / (2 * cap)
  logw <- logw - max(logw)                    # stabilise before exp()
  xs[sample.int(n, 1L, prob = exp(logw))]
}

# Inner tail probability used to estimate a numeric variable's working range
# under DP. The 5th / 95th percentiles are robust to a single extreme record;
# values outside are clamped into the end bins.
DP_DOMAIN_ALPHA <- 0.05

# DP-estimate public-style bin edges c(lo, hi) for one numeric variable, at a
# per-query budget of `eps_q` (two queries: the lower and upper inner quantile).
# Guarantees lo < hi so seq() produces valid edges.
dp_estimate_bounds <- function(x, eps_q, cap) {
  lo <- dp_exp_quantile(x, DP_DOMAIN_ALPHA, eps_q, cap)
  hi <- dp_exp_quantile(x, 1 - DP_DOMAIN_ALPHA, eps_q, cap)
  b <- sort(c(lo, hi))
  if (!all(is.finite(b)) || b[1L] == b[2L]) {
    lo <- if (is.finite(b[1L])) b[1L] else 0
    b <- c(lo, lo + 1)
  }
  b
}

# DP set-union: privately discover the category set of a character column via a
# stability histogram. Each present category's count gets Laplace(cap / eps_op)
# noise; a category is kept only if its noisy count clears a threshold that hides
# any category a single person (contributing <= cap rows) could have created:
#   tau = cap + (cap / eps_op) * log(C * cap / (2 * delta_cat))
# where `C` is the number of character columns sharing the `delta_cat` budget (the
# union bound runs over those columns and over the <= cap categories one person can
# seed in a column). This is (eps_op, delta_cat)-DP by the standard unknown-domain
# histogram argument. Returns the surviving levels, sorted; NA is dropped. Needs
# delta_cat > 0 (a threshold cannot hide a lone category's presence at delta = 0).
dp_discover_categories <- function(x, eps_op, delta_cat, C, cap) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(character(0))
  tab <- table(x)
  counts <- as.numeric(tab)
  noisy <- counts + rlaplace(length(counts), scale = cap / eps_op)
  tau <- cap + (cap / eps_op) * log((C * cap) / (2 * delta_cat))
  sort(names(tab)[noisy >= tau])
}

# Build the full domain (a named list of column-domain objects) for `vars`.
# `est_bounds` supplies DP-estimated numeric ranges (from dp_estimate_bounds);
# together with the public `dp$bounds` these cover every numeric variable in the
# rigorous modes. A warning fires only in the legacy `domain = "data"` mode,
# where some edges are still read from the data and left out of the accounting.
dp_build_domain <- function(data, vars, dp, est_bounds = NULL, est_levels = NULL) {
  dom <- stats::setNames(vector("list", length(vars)), vars)
  for (v in vars) {
    bound <- dp$bounds[[v]]
    if (is.null(bound)) bound <- est_bounds[[v]]
    dom[[v]] <- dp_domain_column(data[[v]], dp$bins, bound, est_levels[[v]])
  }
  derived <- vapply(dom, function(d) isTRUE(d$derived), logical(1))
  if (any(derived)) {
    warning(sprintf(
      paste0("DP: bin edges / category sets for %s were derived from the data. ",
             "These are not included in the (epsilon, delta) accounting; supply ",
             "`bounds` (and factor levels) from public knowledge, or use the ",
             "default `domain = \"dp\"`, for a fully accounted release."),
      paste(names(dom)[derived], collapse = ", ")), call. = FALSE)
  }
  dom
}

# Encode a column's values to integer cell indices in 1..nbin (clamped).
dp_encode <- function(col_dom, v) {
  if (col_dom$kind == "numeric") {
    x <- as.numeric(v)
    x[is.na(x)] <- col_dom$lo
    x <- pmin(pmax(x, col_dom$lo), col_dom$hi)
    code <- findInterval(x, col_dom$edges, rightmost.closed = TRUE,
                         all.inside = TRUE)
    as.integer(pmin(pmax(code, 1L), col_dom$nbin))
  } else {
    key <- if (col_dom$kind == "logical") as.character(as.logical(v))
           else as.character(v)
    code <- match(key, col_dom$levels)
    # Unseen / NA -> the "(other)" catch-all for a DP-discovered character domain,
    # else the first cell (public factor / logical / legacy character).
    unseen <- if (!is.null(col_dom$other_cell)) col_dom$other_cell else 1L
    code[is.na(code)] <- unseen
    as.integer(code)
  }
}

# Decode integer cell indices back to typed values. Numeric cells draw a value
# uniformly within the bin (rounded for integer columns).
dp_decode <- function(col_dom, codes) {
  n <- length(codes)
  if (col_dom$kind == "numeric") {
    lo <- col_dom$edges[codes]
    hi <- col_dom$edges[codes + 1L]
    val <- stats::runif(n, lo, hi)
    if (col_dom$is_integer) {
      val <- round(val)
      val <- pmin(pmax(val, round(col_dom$lo)), round(col_dom$hi))
      storage.mode(val) <- "integer"
    }
    val
  } else if (col_dom$kind == "logical") {
    col_dom$levels[codes] == "TRUE"
  } else if (col_dom$kind == "factor") {
    factor(col_dom$levels[codes], levels = col_dom$levels)
  } else {
    col_dom$levels[codes]
  }
}
