# Internal: within-table longitudinal (DP Markov) modelling for a linked child
# table. Not exported. Enabled by dp_control(longitudinal = TRUE | c(names)).
#
# Without it a child table's repeated rows are exchangeable: each of a parent
# unit's children is drawn i.i.d. from the child's marginals, so a patient's lab
# values across visits carry no autocorrelation. This module models those rows as
# a first-order Markov trajectory under differential privacy - the linked analogue
# of the flat DP longitudinal engine (dp-longitudinal.R), with the parent unit
# playing the role of the "person" and the child's own key index the time index.
#
# For a longitudinally-modelled child table T (immediate parent P, branching cap
# lc = local_cap[T], parent path cap pcp = path_cap[P]) we measure, in place of
# T's plain variable marginals:
#   * an initial-state model over each unit's FIRST (temporally earliest) child
#     row - one-way (+ pairwise, for a Chow-Liu tree) marginals, person
#     sensitivity pcp (an entity owns <= pcp parent units, so <= pcp first rows);
#   * one transition matrix P(v_t | v_{t-1}) per variable over consecutive
#     within-unit rows, person sensitivity pcp * (lc - 1) (each of an entity's
#     <= pcp units of length <= lc contributes <= lc - 1 consecutive pairs).
# The children-per-parent COUNT histogram (already released, sensitivity pcp)
# doubles as the trajectory-length model, so no extra length histogram is added.
# Because every histogram concatenates, these fold into the same exact composition
# as the rest of the release (dp_make_noise): the heterogeneous sensitivities are
# captured entirely in the release totals total_l1 / sum_sq, and one uniform
# per-cell noise draw is applied to every histogram (initial marginals,
# transitions and the count model alike).
#
# Rows must be prefix-truncated in temporal order (not randomly) so consecutive
# pairs are intact; dp_cap_hierarchy() does this for longitudinal tables.

# Resolve dp$longitudinal (FALSE / TRUE / character table names) into a named
# logical vector over the hierarchy: which child tables get a within-unit DP
# Markov model. A table is eligible when it is a child, its branching cap is
# >= 2 (so a transition can be measured) and it has >= 1 modellable variable.
# `TRUE` selects every eligible child; a character vector selects exactly the
# named tables and errors if any is unknown, a root, cap < 2, or variable-free.
dp_resolve_longitudinal <- function(dp, hierarchy, caps, tables) {
  order <- hierarchy$order
  use   <- stats::setNames(logical(length(order)), order)
  req   <- dp$longitudinal
  if (is.logical(req) && !isTRUE(req)) return(use)          # FALSE / default
  eligible <- function(t) {
    !is.na(hierarchy$parent[[t]]) && isTRUE(caps$local[[t]] >= 2L) &&
      length(dp_table_vars(hierarchy, t, tables)) > 0L
  }
  if (isTRUE(req)) {
    for (t in order) use[[t]] <- eligible(t)
    return(use)
  }
  unknown <- setdiff(req, order)
  if (length(unknown)) {
    stop(sprintf(paste0("linked DP longitudinal: unknown table(s) %s. Tables in ",
                        "this release: %s."), paste(unknown, collapse = ", "),
                 paste(order, collapse = ", ")), call. = FALSE)
  }
  for (t in req) {
    if (is.na(hierarchy$parent[[t]])) {
      stop(sprintf(paste0("linked DP longitudinal: table '%s' is a root; only a ",
                          "child table has a within-unit trajectory to model."),
                   t), call. = FALSE)
    }
    if (!isTRUE(caps$local[[t]] >= 2L)) {
      stop(sprintf(paste0("linked DP longitudinal: table '%s' needs a branching ",
                          "cap >= 2 (max children per parent) so within-unit ",
                          "transitions can be measured. Set its ",
                          "`max_rows_per_person` in dp_control() to >= 2."), t),
           call. = FALSE)
    }
    if (length(dp_table_vars(hierarchy, t, tables)) == 0L) {
      stop(sprintf(paste0("linked DP longitudinal: table '%s' has no modellable ",
                          "variable (only keys) to model over time."), t),
           call. = FALSE)
    }
    use[[t]] <- TRUE
  }
  use
}

# Number of initial-state marginals a longitudinally-modelled child contributes:
# one-way for every variable, plus all pairwise when a Chow-Liu tree is requested
# and there are >= 2 variables. (Transitions are counted separately - one per
# variable - because their sensitivity differs.)
dp_longi_n_init <- function(nC, dp) {
  if (dp$dependence == "tree" && nC > 1L) nC + nC * (nC - 1L) / 2L else nC
}

# Fit the initial-state + transition model for one longitudinal child table from
# its capped rows. `cdata_t` is the (already prefix-capped) child data frame;
# `fk` the foreign-key columns identifying the parent unit; `own` the child's own
# key index (the temporal order within a unit); `dom` / `vars` / `nbins` the
# child's discretisation. Returns a model consumed by dp_markov_codes().
#
# When `parent_ctx` / `parent_nbins` are supplied (combined cross-table +
# longitudinal, dp_control(cross_table = TRUE) on a longitudinal child) the
# INITIAL-STATE model is cross-conditioned on the synthetic parent's attributes:
# the first row of each unit draws from a seeded Chow-Liu tree over the parent's
# variables plus the child's, instead of the child's own marginals. `parent_ctx`
# is the named per-parent-var code list aligned to `cdata_t`'s original rows (from
# dp_parent_ctx_codes()); it is reordered to the temporal layout and restricted to
# first rows here. The transition model is unchanged - the parent shapes where a
# trajectory starts, and within-unit autocorrelation carries that dependence
# forward. Cross-conditioning the initial state adds parent-by-child joints at the
# same first-row (parent path-cap) sensitivity as the initial marginals, so it
# folds into the same exact composition.
dp_fit_child_longitudinal <- function(cdata_t, fk, own, dom, vars, nbins,
                                      dp, calib,
                                      parent_ctx = NULL, parent_nbins = NULL) {
  # Order rows so each parent unit's children are contiguous and temporal.
  ord <- do.call(order, c(lapply(fk, function(c) cdata_t[[c]]),
                          list(cdata_t[[own]])))
  ot   <- cdata_t[ord, , drop = FALSE]
  unit <- key_string(ot, fk)
  pos  <- stats::ave(seq_along(unit), match(unit, unique(unit)), FUN = seq_along)

  codes <- stats::setNames(
    lapply(vars, function(v) dp_encode(dom[[v]], ot[[v]])), vars)

  first_rows <- which(pos == 1L)
  cur_rows   <- which(pos > 1L)
  prev_rows  <- cur_rows - 1L                 # contiguous & ordered -> same unit

  # Initial-state model over first rows. Plain: the flat marginal fitter, which
  # honours dp$dependence exactly like the root/child variable models. Cross: the
  # parent-conditioned child fitter over the same first rows, with the parent
  # context aligned to the temporal layout and restricted to first rows.
  init_codes <- stats::setNames(
    lapply(vars, function(v) codes[[v]][first_rows]), vars)
  init_cross <- !is.null(parent_ctx) && length(parent_ctx) > 0L
  if (init_cross) {
    pctx_first <- lapply(parent_ctx, function(x) x[ord][first_rows])
    init_model <- dp_fit_child_cross(init_codes, nbins, pctx_first, parent_nbins,
                                     dp, calib)
  } else {
    init_model <- dp_fit_model(init_codes, nbins, dp, calib)
  }

  # Transition matrices over consecutive within-unit pairs.
  if (length(cur_rows)) {
    prev_codes <- stats::setNames(
      lapply(vars, function(v) codes[[v]][prev_rows]), vars)
    cur_codes  <- stats::setNames(
      lapply(vars, function(v) codes[[v]][cur_rows]), vars)
    tran <- dp_fit_transitions(prev_codes, cur_codes, nbins, calib$add_noise)
  } else {
    # No unit kept >= 2 rows: transitions unidentified but unused (every drawn
    # length collapses to 1). Fall back to a uniform kernel.
    tran <- stats::setNames(
      lapply(vars, function(v) matrix(1 / nbins[[v]], nbins[[v]], nbins[[v]])),
      vars)
  }

  list(kind = "child-longi", vars = vars, nbins = nbins,
       init_model = init_model, tran = tran, init_cross = init_cross)
}
