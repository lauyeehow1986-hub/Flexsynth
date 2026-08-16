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
#
# `held` (a logical vector over `vars`, default all FALSE) marks subject-invariant
# BASELINE columns (dp_control(baseline = ...)): they are modelled once in the
# initial state - so their joint distribution and their correlation with the first
# row are kept - and then held exactly constant within a unit, contributing NO
# transition histogram. `ord` / `cross` (dp_control(transition_order /
# transition_cross)) deepen the transition model of the time-varying columns: each
# steps on its own last `ord` values plus the lag-1 values of its `cross` most
# strongly associated companions (chosen budget-neutrally from the initial model's
# pairwise mutual information). Extra conditioning is budget-free (a tuple lands in
# one cell); a higher order lowers the transition sensitivity to path_cap * (cap -
# ord). These mirror the flat DP longitudinal engine (dp-longitudinal.R), applied
# per child table.
#
# `tran_parent` (dp_control(transition_parent), requires the cross-conditioned
# initial state so `parent_ctx` is present) additionally conditions each
# time-varying column's transition on its `tran_parent` most strongly associated
# immediate-PARENT attributes - reusing the parent-by-child joints the init model
# already measured (init_model$parent_mi) to pick them, so it stays budget-neutral
# and re-anchors the parent dependence at every step instead of only at the start.
dp_fit_child_longitudinal <- function(cdata_t, fk, own, dom, vars, nbins,
                                      dp, calib,
                                      parent_ctx = NULL, parent_nbins = NULL,
                                      held = NULL, ord = 1L, cross = 0L,
                                      tran_parent = 0L) {
  # Order rows so each parent unit's children are contiguous and temporal.
  ord_ix <- do.call(order, c(lapply(fk, function(c) cdata_t[[c]]),
                             list(cdata_t[[own]])))
  ot   <- cdata_t[ord_ix, , drop = FALSE]
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
    pctx_first <- lapply(parent_ctx, function(x) x[ord_ix][first_rows])
    init_model <- dp_fit_child_cross(init_codes, nbins, pctx_first, parent_nbins,
                                     dp, calib)
  } else {
    init_model <- dp_fit_model(init_codes, nbins, dp, calib)
  }

  # Baseline (held) columns get no transition; only the time-varying columns do.
  if (is.null(held)) held <- logical(length(vars))
  tv_vars <- vars[!held]
  nT      <- length(tv_vars)

  # Parent-attribute conditioning of the transitions (transition_parent) needs the
  # cross-conditioned initial state: it reuses that model's parent-by-child joints
  # to pick each variable's parent predictors, and its parent context (reordered to
  # the temporal layout) to look them up. Silently inert when the init is not
  # cross-conditioned (parent_ctx absent) - the caller only sets tran_parent > 0
  # when the child is cross-conditioned.
  tp <- if (init_cross) as.integer(tran_parent) else 0L
  parent_parents <- NULL; parent_ctx_ord <- NULL
  if (tp > 0L && nT > 0L) {
    parent_parents <- dp_select_parent_lags(init_model$parent_mi, tv_vars, tp)
    parent_ctx_ord <- lapply(parent_ctx, function(x) x[ord_ix])
  }

  # Transition model over consecutive within-unit tuples for the time-varying
  # columns. A first-order own-lag-only model uses the simple per-variable
  # matrices; a higher order, any cross-parent, or any parent-attribute uses
  # conditional tensors. Both draw from `calib$add_noise` at the totals already
  # composed by the caller.
  use_tensor <- (ord > 1L || cross > 0L || tp > 0L) && nT > 0L
  tran <- NULL; tensors <- NULL; cross_parents <- NULL
  if (use_tensor) {
    cross_parents <- dp_select_cross_parents(init_model$pairwise_mi, tv_vars,
                                             vars, cross)
    tensors <- dp_fit_transition_tensors(codes, nbins, pos, ord, cross_parents,
                                         tv_vars, calib$add_noise,
                                         parent_ctx = parent_ctx_ord,
                                         parent_parents = parent_parents,
                                         parent_nbins = parent_nbins)
  } else if (nT && length(cur_rows)) {
    prev_codes <- stats::setNames(
      lapply(tv_vars, function(v) codes[[v]][prev_rows]), tv_vars)
    cur_codes  <- stats::setNames(
      lapply(tv_vars, function(v) codes[[v]][cur_rows]), tv_vars)
    tran <- dp_fit_transitions(prev_codes, cur_codes, nbins[tv_vars],
                               calib$add_noise)
  } else {
    # No unit kept >= 2 rows (or nothing is time-varying): transitions
    # unidentified but unused. Fall back to a uniform kernel over tv columns.
    tran <- stats::setNames(
      lapply(tv_vars, function(v) matrix(1 / nbins[[v]], nbins[[v]], nbins[[v]])),
      tv_vars)
  }

  list(kind = "child-longi", vars = vars, nbins = nbins,
       init_model = init_model, tran = tran, tensors = tensors,
       init_cross = init_cross, held = held, tv_vars = tv_vars,
       use_tensor = use_tensor, order = as.integer(ord),
       cross = as.integer(cross), cross_parents = cross_parents,
       tran_parent = tp, parent_parents = parent_parents)
}
