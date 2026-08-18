# Disclosure-risk diagnostics. Synthetic data is not anonymisation; these
# empirical measures let residual identity / attribute / membership risk be
# judged honestly. None of them is a formal guarantee (that is Track B).

# Coerce a synth object to plain data / list; leave data.frames and lists alone.
as_risk_frame <- function(x) {
  if (inherits(x, "synth_result")) return(as.data.frame(x))
  if (inherits(x, "synth_linked_result")) return(as.list(x))
  x
}

# Row keys: encode every field with an explicit byte length.  A plain separator
# is not sufficient because user values may contain that separator (and a
# literal "NA" must remain distinct from a missing value).
row_keys <- function(df, quasi) {
  parts <- lapply(quasi, function(v) {
    raw <- df[[v]]
    val <- enc2utf8(as.character(raw))
    missing <- is.na(raw)
    encoded <- rep.int("N;", length(val))
    present <- !missing
    encoded[present] <- paste0(
      "V", nchar(val[present], type = "bytes"), ":", val[present], ";"
    )
    encoded
  })
  do.call(paste0, parts)
}

# Sample down to at most `n` rows (reproducible via the caller's seed).
subsample <- function(df, n) {
  if (nrow(df) <= n) return(df)
  df[sort(sample.int(nrow(df), n)), , drop = FALSE]
}

# TCAP: Target Correct Attribution Probability (Taub, Elliot et al.). An attacker
# who knows a real record's quasi-identifier `keys` looks up the synthetic
# records sharing those key values and reads off the conditional distribution of
# the sensitive `target`; the CAP for that record is the synthetic probability of
# its true target value. TCAP averages this over real records that have any
# synthetic key match (`coverage`). It is compared with a marginal-only attacker
# (`baseline`, the synthetic marginal probability of the true target), so the
# `lift` is the excess attribution the key-conditioning buys. Meant for
# categorical keys and target; continuous columns should be coarsened first.
tcap_risk <- function(real, syn, keys, target, max_records) {
  rl <- subsample(real, max_records)
  rk <- row_keys(rl, keys); rt <- as.character(rl[[target]])
  sk <- row_keys(syn, keys); st <- as.character(syn[[target]])

  syn_by_key <- split(st, sk)                  # synthetic targets per key value
  cap <- vapply(seq_along(rk), function(i) {
    g <- syn_by_key[[rk[i]]]
    if (is.null(g)) NA_real_ else mean(g == rt[i])
  }, numeric(1))
  covered <- !is.na(cap)

  marg <- table(st) / length(st)               # marginal-only attacker
  base <- as.numeric(marg[rt]); base[is.na(base)] <- 0

  tcap <- if (any(covered)) mean(cap[covered]) else NA_real_
  baseline_covered <- if (any(covered)) mean(base[covered]) else NA_real_
  list(target = target, keys = keys,
       tcap = tcap, baseline = baseline_covered,
       baseline_unconditional = mean(base),
       lift = tcap - baseline_covered,
       coverage = mean(covered), n_scored = length(rk),
       n_covered = sum(covered))
}

#' Empirical disclosure-risk diagnostics for synthetic data
#'
#' Reports how much a synthetic dataset could leak about the real records it was
#' trained on. Four complementary measures are computed:
#'
#' \itemize{
#'   \item **Replicated uniques** — records that are unique in the real data on
#'     the quasi-identifiers and are nonetheless reproduced exactly in the
#'     synthetic data. These are the classic identity-disclosure risks.
#'   \item **Distance to closest record (DCR)** — for each synthetic record, the
#'     Gower distance to the nearest real record. A DCR of 0 is an exact copy.
#'     The synthetic-to-real distances are compared against the real-to-real
#'     nearest-neighbour distances: if synthetic records are not systematically
#'     closer to real records than real records are to each other, identity risk
#'     is low.
#'   \item **Membership inference** — if a `holdout` of records *not* used for
#'     synthesis is supplied, an attacker who guesses "member" for records close
#'     to the synthetic data is simulated. The reported AUC (0.5 = no advantage)
#'     and advantage (`2 * AUC - 1`) measure how well training membership can be
#'     inferred.
#'   \item **Attribute disclosure (TCAP)** — if a sensitive `target` column is
#'     named, the Target Correct Attribution Probability (Taub, Elliot et al.):
#'     an attacker who knows a real record's quasi-identifier keys reads the
#'     conditional distribution of the target off the synthetic records that share
#'     those keys. `tcap` is the mean synthetic probability of the true target
#'     over matched records, `baseline` the marginal-only attacker evaluated on
#'     those same covered records, and `lift` the excess the key-conditioning
#'     buys. `baseline_unconditional` is also returned for context. Meant for
#'     categorical keys / target.
#' }
#'
#' Quasi-identifiers should be the genuinely identifying columns; exclude
#' surrogate keys such as a regenerated unit id (a synthetic id never matches a
#' real one, which would understate risk). If `real` and `syn` are named lists
#' (or `syn` is a `synth_linked_result`), each table is assessed separately.
#'
#' @param real The real `data.frame` (or named list) the data was synthesised
#'   from.
#' @param syn The synthetic data: a `data.frame`, [synth()] `synth_result`,
#'   named list, or [synth_linked()] `synth_linked_result`.
#' @param quasi Character vector of quasi-identifier columns; defaults to all
#'   columns present in both frames.
#' @param target Optional single sensitive column for the attribute-disclosure
#'   (TCAP) measure; the quasi-identifiers minus `target` are used as the keys.
#'   `NULL` (default) skips it. For a linked (list) input it is assessed only on
#'   the table that contains it.
#' @param holdout Optional `data.frame` of real records *excluded* from
#'   synthesis, enabling the membership-inference check. For linked input, a
#'   named list containing a holdout `data.frame` for every evaluated table.
#' @param max_records Cap on the number of rows used for the distance
#'   computations (each of real / synthetic is sampled down to this); keeps the
#'   pairwise distances tractable. Default 2000.
#' @param seed Optional integer seed for the subsampling.
#' @param ... Unused.
#'
#' @return A `flexsynth_disclosure` object (or `flexsynth_disclosure_list`) with
#'   a `print()` method.
#' @seealso [diagnose()] for utility.
#' @export
#' @examples
#' df <- data.frame(
#'   id  = 1:200,
#'   age = round(rnorm(200, 60, 10)),
#'   sex = sample(c("F", "M"), 200, replace = TRUE)
#' )
#' res <- synth(df, ~ id, seed = 1)
#' disclosure_risk(df, res, quasi = c("age", "sex"))
disclosure_risk <- function(real, syn, quasi = NULL, target = NULL,
                            holdout = NULL, max_records = 2000L, seed = NULL,
                            ...) {
  syn <- as_risk_frame(syn)

  if (is.list(real) && !is.data.frame(real) &&
      is.list(syn) && !is.data.frame(syn)) {
    tbls <- intersect(names(real), names(syn))
    if (!length(tbls)) stop("`real` and `syn` share no table names.", call. = FALSE)
    if (!is.null(holdout)) {
      if (!is.list(holdout) || is.data.frame(holdout))
        stop("for linked data, `holdout` must be a named list of data.frames.",
             call. = FALSE)
      missing_holdout <- setdiff(tbls, names(holdout))
      if (length(missing_holdout))
        stop(sprintf("`holdout` is missing table(s): %s",
                     paste(missing_holdout, collapse = ", ")), call. = FALSE)
    }
    out <- stats::setNames(lapply(tbls, function(t) {
      # A target applies only to the table that contains it.
      tgt <- if (!is.null(target) && target %in% names(real[[t]])) target else NULL
      disclosure_risk(real[[t]], syn[[t]], quasi = quasi, target = tgt,
                      holdout = if (is.null(holdout)) NULL else holdout[[t]],
                      max_records = max_records, seed = seed)
    }), tbls)
    return(structure(out, class = "flexsynth_disclosure_list"))
  }

  if (!is.data.frame(real)) stop("`real` must be a data.frame.", call. = FALSE)
  if (!is.data.frame(syn))  stop("`syn` must be a data.frame or synth result.",
                                 call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  cols <- intersect(names(real), names(syn))
  quasi <- quasi %||% cols
  miss <- setdiff(quasi, cols)
  if (length(miss))
    stop(sprintf("`quasi` not present in both frames: %s",
                 paste(miss, collapse = ", ")), call. = FALSE)
  if (!length(quasi)) stop("no quasi-identifier columns to compare.", call. = FALSE)

  ## --- attribute disclosure (TCAP), if a sensitive target is named --------
  attr_res <- NULL
  if (!is.null(target)) {
    if (length(target) != 1L || !target %in% cols)
      stop(sprintf("`target` must be a single column present in both frames: %s",
                   target), call. = FALSE)
    keys <- setdiff(quasi, target)
    if (!length(keys))
      stop("attribute disclosure needs at least one quasi-identifier key besides the target.",
           call. = FALSE)
    attr_res <- tcap_risk(real, syn, keys, target, max_records)
  }

  ## --- replicated uniques (over the full data, exact match) ---------------
  rk <- row_keys(real, quasi)
  sk <- row_keys(syn, quasi)
  rt <- table(rk)
  unique_keys <- names(rt)[rt == 1L]
  syn_set     <- unique(sk)
  repl_unique <- intersect(unique_keys, syn_set)
  exact_any   <- sum(sk %in% names(rt))          # syn rows matching any real row

  ## --- distances (on a subsample) -----------------------------------------
  rs <- subsample(real, max_records)
  ss <- subsample(syn,  max_records)
  scales <- distance_scales(rs, quasi)
  dcr_syn  <- nn_distance(ss, rs, quasi, scales)              # syn -> real
  dcr_real <- nn_distance(rs, rs, quasi, scales, exclude_self = TRUE)  # real -> real

  ## --- membership inference (needs a holdout) -----------------------------
  mi <- NULL
  if (!is.null(holdout)) {
    if (!is.data.frame(holdout)) stop("`holdout` must be a data.frame.", call. = FALSE)
    if (length(setdiff(quasi, names(holdout))))
      stop("`holdout` is missing quasi-identifier columns.", call. = FALSE)
    hs <- subsample(holdout, max_records)
    members    <- rs
    d_member   <- nn_distance(members, ss, quasi, scales)   # member -> syn
    d_nonmemb  <- nn_distance(hs,      ss, quasi, scales)   # non-member -> syn
    score  <- c(-d_member, -d_nonmemb)                      # closer = more "member"
    is_mem <- c(rep(TRUE, length(d_member)), rep(FALSE, length(d_nonmemb)))
    auc <- auc_mw(score, is_mem)
    mi <- list(auc = auc, advantage = 2 * auc - 1,
               n_member = length(d_member), n_nonmember = length(d_nonmemb))
  }

  structure(
    list(
      quasi = quasi,
      n_real = nrow(real), n_syn = nrow(syn),
      replicated_uniques = list(
        n_real_unique = length(unique_keys),
        n_replicated  = length(repl_unique),
        prop_of_uniques = if (length(unique_keys))
          length(repl_unique) / length(unique_keys) else NA_real_,
        prop_of_real  = length(repl_unique) / nrow(real),
        n_exact_syn   = exact_any,
        prop_exact_syn = exact_any / nrow(syn)
      ),
      dcr = list(
        syn_to_real  = dcr_syn,
        real_to_real = dcr_real,
        prop_syn_zero = mean(dcr_syn == 0),
        median_syn   = stats::median(dcr_syn),
        median_real  = stats::median(dcr_real),
        q05_syn      = stats::quantile(dcr_syn, 0.05, names = FALSE)
      ),
      membership = mi,
      attribute = attr_res
    ),
    class = "flexsynth_disclosure"
  )
}

#' @export
print.flexsynth_disclosure <- function(x, ...) {
  cat("<flexsynth_disclosure>\n")
  cat("  rows            : real", x$n_real, " synthetic", x$n_syn, "\n")
  cat("  quasi-identifiers:", paste(x$quasi, collapse = ", "), "\n\n")

  ru <- x$replicated_uniques
  cat("Replicated uniques (identity risk):\n")
  cat(sprintf("  real sample-uniques : %d\n", ru$n_real_unique))
  cat(sprintf("  reproduced in syn   : %d  (%.2f%% of uniques, %.2f%% of real rows)\n",
              ru$n_replicated, 100 * ru$prop_of_uniques, 100 * ru$prop_of_real))
  cat(sprintf("  syn rows copying a real row: %d  (%.2f%% of syn)\n",
              ru$n_exact_syn, 100 * ru$prop_exact_syn))

  d <- x$dcr
  cat("\nDistance to closest record (Gower, 0 = exact copy):\n")
  cat(sprintf("  syn->real  : median %.4f   5th pct %.4f   exact copies %.2f%%\n",
              d$median_syn, d$q05_syn, 100 * d$prop_syn_zero))
  cat(sprintf("  real->real : median %.4f   (baseline)\n", d$median_real))
  verdict <- if (d$median_syn + 1e-9 >= d$median_real)
    "median syn distance is not smaller than the real-neighbour baseline"
  else
    "median syn distance is smaller than the real-neighbour baseline"
  cat("  ", verdict, "\n", sep = "")
  cat("  descriptive only: inspect lower-tail distances and exact copies; ",
      "this is not a safety guarantee\n", sep = "")

  if (!is.null(x$membership)) {
    m <- x$membership
    cat("\nMembership inference (0.5 AUC = no advantage):\n")
    cat(sprintf("  AUC %.3f   attacker advantage %.3f   (%d members / %d non-members)\n",
                m$auc, m$advantage, m$n_member, m$n_nonmember))
  } else {
    cat("\nMembership inference: not run (supply `holdout` of non-training records).\n")
  }

  if (!is.null(x$attribute)) {
    a <- x$attribute
    cat(sprintf("\nAttribute disclosure (TCAP, target = %s):\n", a$target))
    cat(sprintf("  TCAP %.3f   baseline %.3f   lift %+.3f   (coverage %.1f%%, keys: %s)\n",
                a$tcap, a$baseline, a$lift, 100 * a$coverage,
                paste(a$keys, collapse = ", ")))
    cat("  ", if (isTRUE(a$lift > 0.1))
      "key-conditioning attributes the target well above the covered-record margin -> inspect"
      else "little lift over the covered-record marginal baseline in this sample",
      "\n", sep = "")
  } else {
    cat("\nAttribute disclosure: not run (supply `target` = a sensitive column).\n")
  }
  invisible(x)
}

#' @export
print.flexsynth_disclosure_list <- function(x, ...) {
  cat("<flexsynth_disclosure_list>", length(x), "tables\n\n")
  for (nm in names(x)) {
    cat("== ", nm, " ==\n", sep = "")
    print(x[[nm]])
    cat("\n")
  }
  invisible(x)
}
