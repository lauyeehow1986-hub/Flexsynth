# Internal: the Track A sequential-synthesis engine (synthpop lineage).
# Not exported. Everything here works on plain data.frames.
#
# Phase 1 gave flat + subject-level-baseline synthesis. Phase 2 adds two things:
#   * an explicit structural (rows-per-unit) count model, so the synthetic
#     skeleton is drawn from the learned size distribution rather than by
#     copying whole real units, and
#   * within-unit temporal dependence for time-varying variables, via an
#     initial-state model (first row of each unit) plus a Markov transition
#     model that conditions on the previous row (lag-1 predictors), so
#     autocorrelation across visits is preserved.

# ---------------------------------------------------------------------------
# Method resolution
# ---------------------------------------------------------------------------

# Resolve the per-variable method vector from the call-level `method` and any
# per-variable override carried on the control object.
resolve_methods <- function(method, per_var, cols) {
  m <- stats::setNames(rep(method, length(cols)), cols)
  if (length(cols) && !is.null(per_var)) {
    if (length(per_var) == 1L && is.null(names(per_var))) {
      m[] <- per_var
    } else {
      common <- intersect(names(per_var), cols)
      m[common] <- per_var[common]
    }
  }
  supported <- list_methods()
  bad <- setdiff(unique(m), supported)
  if (length(bad)) {
    stop(sprintf(
      "unsupported method(s): %s (registered: %s).",
      paste(bad, collapse = ", "), paste(supported, collapse = ", ")
    ), call. = FALSE)
  }
  m
}

# ---------------------------------------------------------------------------
# Structural / count model + skeleton
# ---------------------------------------------------------------------------

# Row-ordering key: sort by unit identifier, then by each structural index in
# turn, so rows within a unit are in temporal / positional order.
order_rows <- function(data, st) {
  keys <- c(list(data[[st$id]]), lapply(st$nested, function(k) data[[k]]))
  do.call(order, keys)
}

# Learn the rows-per-unit ("count") distribution empirically and draw counts for
# the synthetic units. With `k` the loop accumulates whole units until the row
# budget is met. `sizes` is sampled by position (never with sample()'s scalar
# shortcut) so a length-one size set is handled correctly.
draw_counts <- function(sizes, n_units, k) {
  draw1 <- function(m) sizes[sample.int(length(sizes), m, replace = TRUE)]
  if (is.null(k)) return(draw1(n_units))
  counts <- integer(0)
  total <- 0L
  while (total < k) {
    cc <- draw1(1L)
    counts <- c(counts, cc)
    total <- total + cc
  }
  counts
}

# Build a synthetic skeleton (identifier + structural indices). The number of
# rows per synthetic unit is an explicit draw from the learned count model; the
# structural-index sequence for a unit is then drawn from a real unit of that
# size (so the index values stay in-distribution and correctly ordered).
synth_skeleton <- function(data, st, control) {
  id <- st$id
  ord <- order_rows(data, st)
  rdat <- data[ord, , drop = FALSE]
  blocks <- split(seq_len(nrow(rdat)), rdat[[id]])   # row indices per unit, in order
  sizes <- lengths(blocks)
  by_size <- split(seq_along(blocks), sizes)         # unit positions grouped by size

  counts <- draw_counts(sizes, length(blocks), control$k)

  keep <- c(id, st$nested)
  parts <- vector("list", length(counts))
  for (j in seq_along(counts)) {
    pool <- by_size[[as.character(counts[[j]])]]     # real units with this many rows
    u <- pool[sample.int(length(pool), 1L)]
    block <- rdat[blocks[[u]], keep, drop = FALSE]
    block[[id]] <- j                                 # fresh sequential identifier
    parts[[j]] <- block
  }
  skel <- rbind_rows(parts)
  if (!is.null(control$k) && nrow(skel) > control$k) {
    skel <- skel[seq_len(control$k), , drop = FALSE]
  }
  rownames(skel) <- NULL
  skel
}

# ---------------------------------------------------------------------------
# CART primitives (fit / apply are split so a fitted model can be reused across
# time positions without refitting — essential for the transition model).
# ---------------------------------------------------------------------------

# Coerce character/logical predictors to factors; optionally align factor levels
# to a reference frame so rpart does not choke on unseen levels in new data.
prep_predictors <- function(x, levels_from = NULL) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    if (is.character(x[[nm]]) || is.logical(x[[nm]])) {
      x[[nm]] <- factor(x[[nm]])
    }
  }
  if (!is.null(levels_from)) {
    for (nm in names(x)) {
      if (is.factor(levels_from[[nm]])) {
        x[[nm]] <- factor(x[[nm]], levels = levels(levels_from[[nm]]))
      }
    }
  }
  x
}

# Leaf (terminal-node) id for each training row.
train_leaves <- function(fit) {
  as.integer(rownames(fit$frame))[fit$where]
}

# Leaf id for each row of new data. Overwriting `yval` with the node numbers and
# asking for the "vector" prediction makes rpart return the terminal-node id for
# both regression and classification trees.
predict_leaves <- function(fit, newdata) {
  fit$frame$yval <- as.numeric(rownames(fit$frame))
  as.integer(predict(fit, newdata = newdata, type = "vector"))
}

# Draw a synthetic value for each new row by sampling a real observation from the
# training rows that share its terminal node (conditional bootstrap). Preserves
# the target's type (numeric / factor / character) automatically.
draw_from_leaves <- function(y, tr_leaf, syn_leaf) {
  by_leaf <- split(seq_along(y), tr_leaf)
  idx <- integer(length(syn_leaf))
  groups <- split(seq_along(syn_leaf), syn_leaf)
  for (lf in names(groups)) {
    pool <- by_leaf[[lf]]
    if (is.null(pool)) pool <- seq_along(y)   # defensive: unseen node
    pos <- groups[[lf]]
    idx[pos] <- pool[sample.int(length(pool), length(pos), replace = TRUE)]
  }
  y[idx]
}

# Fit a CART model of y on the given predictor frame. Returns everything
# apply-time needs (the tree, the training target, per-row leaf ids, and the
# training frame for factor-level alignment). `bootstrap` resamples the training
# rows once, here, so the model is stable across every later apply; it defaults
# to `proper` (posterior-predictive approximation) but the forest turns it off
# because bagging already supplies the resampling.
cart_fit <- function(y, xtrain, control, bootstrap = isTRUE(control$proper)) {
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("method = \"cart\" needs the 'rpart' package. install.packages(\"rpart\").",
         call. = FALSE)
  }
  xtrain <- prep_predictors(xtrain)
  yy <- y
  if (isTRUE(bootstrap)) {                      # posterior-predictive approximation
    bi <- sample.int(length(y), length(y), replace = TRUE)
    yy <- y[bi]
    xtrain <- xtrain[bi, , drop = FALSE]
  }
  method <- if (is.numeric(yy)) "anova" else "class"
  ctrl_args <- utils::modifyList(list(minbucket = 5L, cp = 1e-4), control$cart)
  df <- data.frame(.y = yy, xtrain, check.names = FALSE)
  fit <- rpart::rpart(.y ~ ., data = df, method = method,
                      control = do.call(rpart::rpart.control, ctrl_args))
  list(fit = fit, y = yy, tr_leaf = train_leaves(fit), xref = xtrain)
}

# Draw synthetic values for new rows from a model returned by cart_fit().
cart_apply <- function(model, xsyn) {
  xsyn <- prep_predictors(xsyn, levels_from = model$xref)
  draw_from_leaves(model$y, model$tr_leaf, predict_leaves(model$fit, xsyn))
}

# Convenience: fit + apply in one shot (used by the subject-grain sequence).
cart_draw <- function(y, xtrain, xsyn, control) {
  cart_apply(cart_fit(y, xtrain, control), xsyn)
}

# Unconditional (empirical) draw; `proper` adds a bootstrap of the donor pool.
sample_draw <- function(y, n, proper = FALSE) {
  pool <- if (isTRUE(proper)) y[sample.int(length(y), length(y), replace = TRUE)] else y
  pool[sample.int(length(pool), n, replace = TRUE)]
}

# A per-variable synthesiser is resolved from the method registry (see
# methods.R): `fit_var` builds a fitted model (falling back to the marginal
# sampler when a predictor-hungry method has no predictors), and `apply_var`
# draws from it. Both live in methods.R so the registry is the single source of
# truth for method dispatch.

# ---------------------------------------------------------------------------
# Column-role detection
# ---------------------------------------------------------------------------

# Which of `cols` are constant within every unit (i.e. subject-level / baseline
# rather than time-varying)? Those are synthesised once per unit.
subject_level_cols <- function(data, id, cols) {
  if (!length(cols)) return(character(0))
  cols[constant_within(data, data[[id]], cols)]
}

# ---------------------------------------------------------------------------
# Sequence synthesis (subject grain: no lags, one row per unit)
# ---------------------------------------------------------------------------

# Synthesise `cols` sequentially at a given grain (`train` = real rows at that
# grain, `syn` = the synthetic frame to fill). Each variable is drawn from its
# registered method, conditioning on the columns synthesised before it (subject
# to any `predictor_matrix` restriction). Returns `syn` with the new cols.
synth_sequence <- function(train, syn, cols, available, methods, control) {
  n_syn  <- nrow(syn)
  pm     <- control$predictor_matrix
  smooth <- smoothing_targets(control, cols)
  for (v in cols) {
    y <- train[[v]]
    preds <- allowed_predictors(pm, v, available)
    model <- fit_var(y, train[preds], preds, methods[[v]], control)
    val   <- apply_var(model, syn[preds], n_syn, control)
    if (v %in% smooth) val <- smooth_draw(val, y)
    syn[[v]] <- val
    available <- c(available, v)
  }
  syn
}

# ---------------------------------------------------------------------------
# Temporal synthesis (row grain: lag-1 predictors, autoregressive within unit)
# ---------------------------------------------------------------------------

# Lag-1 of `x` within each unit (NA at the first row of a unit). Assumes `ids`
# is sorted so that a unit's rows are contiguous and in positional order.
unit_lag <- function(x, ids) {
  n <- length(x)
  idx <- c(NA_integer_, seq_len(n - 1L))
  idx[!duplicated(ids)] <- NA_integer_
  x[idx]
}

# Assemble a predictor frame from an ordered vector of predictor names. Plain
# names are pulled from `nonlag_df` at `rows`; names of the form ".lag_<v>" are
# supplied by `lag_fun("<v>")`. Building the training and the generation frames
# through this single assembler guarantees identical columns and column order,
# which matters for the parametric (model-matrix) methods.
assemble_frame <- function(pred, nonlag_df, rows, lag_fun = NULL) {
  cols <- lapply(pred, function(p) {
    if (startsWith(p, ".lag_")) lag_fun(sub("^\\.lag_", "", p))
    else nonlag_df[[p]][rows]
  })
  as.data.frame(stats::setNames(cols, pred), stringsAsFactors = FALSE,
                check.names = FALSE)
}

# Fill the time-varying columns of `syn` with an initial-state + Markov
# transition model. For each variable we fit two models on the real data:
#   * initial  — the first row of each unit, conditioned on baseline + structural
#                indices + earlier-in-sequence current-row variables;
#   * transition — every later row, additionally conditioned on the lag-1 value
#                of every time-varying variable (autocorrelation + cross-lags).
# Generation then proceeds position by position (t = 1, 2, ...): the lag columns
# for position t are read from the already-synthesised rows at position t - 1.
synth_temporal <- function(data, syn, st, subj_cols, time_cols, fixed_cols,
                           methods, control) {
  id     <- st$id
  pm     <- control$predictor_matrix
  smooth <- smoothing_targets(control, time_cols)

  ## Real side: sort, positions, and lag-1 predictors.
  rdat <- data[order_rows(data, st), , drop = FALSE]
  rids <- rdat[[id]]
  rpos <- ave(seq_len(nrow(rdat)), rids, FUN = seq_along)
  first <- rpos == 1L
  nonfirst <- !first
  have_nonfirst <- any(nonfirst)

  lag_names <- paste0(".lag_", time_cols)
  rlags <- stats::setNames(
    lapply(time_cols, function(v) unit_lag(rdat[[v]], rids)), lag_names
  )

  ## Per-variable predictor sets (structural indices in `fixed_cols` are always
  ## kept; baseline, earlier current-row and lag predictors honour any
  ## predictor_matrix restriction). Stored so generation rebuilds the exact same
  ## columns the model was fitted on.
  init_pred <- stats::setNames(vector("list", length(time_cols)), time_cols)
  tran_pred <- stats::setNames(vector("list", length(time_cols)), time_cols)
  init_models <- stats::setNames(vector("list", length(time_cols)), time_cols)
  tran_models <- stats::setNames(vector("list", length(time_cols)), time_cols)

  earlier <- character(0)
  for (v in time_cols) {
    meth  <- methods[[v]]
    subj_a    <- allowed_predictors(pm, v, subj_cols)
    earlier_a <- allowed_predictors(pm, v, earlier)
    lag_a     <- paste0(".lag_", allowed_predictors(pm, v, time_cols))
    ipred <- c(subj_a, fixed_cols, earlier_a)
    tpred <- c(ipred, lag_a)
    init_pred[[v]] <- ipred
    tran_pred[[v]] <- tpred

    iframe <- assemble_frame(ipred, rdat, which(first))
    init_models[[v]] <- fit_var(rdat[[v]][first], iframe, ipred, meth, control)
    if (have_nonfirst) {
      rows_nf <- which(nonfirst)
      tframe <- assemble_frame(tpred, rdat, rows_nf,
                               lag_fun = function(w) rlags[[paste0(".lag_", w)]][rows_nf])
      tran_models[[v]] <- fit_var(rdat[[v]][nonfirst], tframe, tpred, meth, control)
    }
    earlier <- c(earlier, v)
  }

  ## Generate: syn is grouped and ordered by unit (from synth_skeleton).
  sids <- syn[[id]]
  spos <- ave(seq_len(nrow(syn)), sids, FUN = seq_along)
  prev_idx <- c(NA_integer_, seq_len(nrow(syn) - 1L))
  prev_idx[!duplicated(sids)] <- NA_integer_

  for (v in time_cols) syn[[v]] <- rdat[[v]][rep(NA_integer_, nrow(syn))]  # typed NA

  for (t in seq_len(max(spos))) {
    rows_t <- which(spos == t)
    if (!length(rows_t)) next
    use_transition <- t > 1L && have_nonfirst
    for (v in time_cols) {
      if (use_transition) {
        model <- tran_models[[v]]
        xsyn  <- assemble_frame(tran_pred[[v]], syn, rows_t,
                                lag_fun = function(w) syn[[w]][prev_idx[rows_t]])
      } else {
        model <- init_models[[v]]
        xsyn  <- assemble_frame(init_pred[[v]], syn, rows_t)
      }
      val <- apply_var(model, xsyn, length(rows_t), control)
      if (v %in% smooth) val <- smooth_draw(val, rdat[[v]])
      syn[[v]][rows_t] <- val
    }
  }
  syn
}

# ---------------------------------------------------------------------------
# Column filling (shared by single-table and linked-child synthesis)
# ---------------------------------------------------------------------------

# Fill the variable columns of an already-built skeleton `syn` (which carries the
# unit identifier, structural indices, and any predictor columns attached ahead
# of time). Subject-invariant columns are synthesised once per unit and
# broadcast; time-varying columns get the autoregressive temporal model.
#   * `subj_fixed` — predictor columns available at the subject grain (constant
#     within a unit, e.g. a linked child's parent attributes). For a plain single
#     table this is empty.
#   * `time_fixed` — predictor columns available at the row grain (structural
#     indices, plus any carried attributes).
# Both must already be present in `train` and in `syn`.
fill_columns <- function(train, syn, st, subj_cols, time_cols,
                         subj_fixed, time_fixed, methods, control) {
  id_col <- st$id

  if (length(subj_cols)) {
    keepc <- c(id_col, subj_fixed)
    sdat <- train[!duplicated(train[[id_col]]), c(keepc, subj_cols), drop = FALSE]
    ssyn <- syn[!duplicated(syn[[id_col]]), keepc, drop = FALSE]
    ssyn <- synth_sequence(sdat, ssyn, subj_cols, subj_fixed, methods, control)
    pos <- match(syn[[id_col]], ssyn[[id_col]])
    for (v in subj_cols) syn[[v]] <- ssyn[[v]][pos]
  }

  if (length(time_cols)) {
    syn <- synth_temporal(train, syn, st, subj_cols, time_cols, time_fixed,
                          methods, control)
  }
  syn
}

# ---------------------------------------------------------------------------
# One synthetic dataset (single table)
# ---------------------------------------------------------------------------

# Produce one synthetic data.frame. The skeleton (unit id + structural indices)
# comes from the learned count model; subject-invariant columns are synthesised
# once per unit and broadcast; remaining time-varying columns are synthesised
# with the autoregressive temporal model.
synthesise_once <- function(data, st, subj_cols, time_cols, fixed_cols,
                            methods, control) {
  syn <- synth_skeleton(data, st, control)
  syn <- fill_columns(data, syn, st, subj_cols, time_cols,
                      subj_fixed = character(0), time_fixed = fixed_cols,
                      methods, control)
  syn <- syn[names(data)]          # restore original column order
  rownames(syn) <- NULL
  syn
}
