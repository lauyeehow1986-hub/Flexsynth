# Internal: the Track A sequential-synthesis engine (synthpop lineage).
# Not exported. Everything here works on plain data.frames.

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
  supported <- c("cart", "sample")
  bad <- setdiff(unique(m), supported)
  if (length(bad)) {
    stop(sprintf(
      "unsupported method(s) for Phase 1: %s (supported: %s).",
      paste(bad, collapse = ", "), paste(supported, collapse = ", ")
    ), call. = FALSE)
  }
  m
}

# Build a synthetic "skeleton": the identifier + structural columns, produced by
# resampling whole units (subject bootstrap) so realistic per-unit block sizes
# are preserved. Measurement columns are synthesised on top of this later.
build_skeleton <- function(data, st, k) {
  id_col <- st$id
  blocks <- split(seq_len(nrow(data)), data[[id_col]])
  n_subj <- length(blocks)

  if (is.null(k)) {
    chosen <- sample.int(n_subj, n_subj, replace = TRUE)
  } else {
    chosen <- integer(0)
    total <- 0L
    while (total < k) {
      s <- sample.int(n_subj, 1L)
      chosen <- c(chosen, s)
      total <- total + length(blocks[[s]])
    }
  }

  keep <- c(id_col, st$nested)
  parts <- vector("list", length(chosen))
  for (j in seq_along(chosen)) {
    block <- data[blocks[[chosen[j]]], keep, drop = FALSE]
    block[[id_col]] <- j                 # fresh sequential identifier
    parts[[j]] <- block
  }
  skel <- do.call(rbind, parts)
  if (!is.null(k) && nrow(skel) > k) {
    skel <- skel[seq_len(k), , drop = FALSE]
  }
  rownames(skel) <- NULL
  skel
}

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

# One CART-based conditional draw: fit y ~ predictors on the real data, then draw
# a synthetic value per new row from the matching leaf.
cart_draw <- function(y, xtrain, xsyn, control) {
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("method = \"cart\" needs the 'rpart' package. install.packages(\"rpart\").",
         call. = FALSE)
  }
  xtrain <- prep_predictors(xtrain)
  xsyn   <- prep_predictors(xsyn, levels_from = xtrain)

  yy <- y
  if (isTRUE(control$proper)) {                # posterior-predictive approximation
    bi <- sample.int(length(y), length(y), replace = TRUE)
    yy <- y[bi]
    xtrain <- xtrain[bi, , drop = FALSE]
  }

  method <- if (is.numeric(yy)) "anova" else "class"
  ctrl_args <- utils::modifyList(list(minbucket = 5L, cp = 1e-4), control$cart)
  df <- data.frame(.y = yy, xtrain, check.names = FALSE)
  fit <- rpart::rpart(.y ~ ., data = df, method = method,
                      control = do.call(rpart::rpart.control, ctrl_args))

  draw_from_leaves(yy, train_leaves(fit), predict_leaves(fit, xsyn))
}

# Unconditional (empirical) draw; `proper` adds a bootstrap of the donor pool.
sample_draw <- function(y, n, proper = FALSE) {
  pool <- if (isTRUE(proper)) y[sample.int(length(y), length(y), replace = TRUE)] else y
  pool[sample.int(length(pool), n, replace = TRUE)]
}

# Which of `cols` are constant within every unit (i.e. subject-level / baseline
# rather than time-varying)? Those are synthesised once per unit.
subject_level_cols <- function(data, id, cols) {
  if (!length(cols)) return(character(0))
  gid <- data[[id]]
  keep <- vapply(cols, function(col) {
    all(tapply(data[[col]], gid, function(v) length(unique(v)) == 1L))
  }, logical(1))
  cols[keep]
}

# Synthesise `cols` sequentially at a given grain (`train` = real rows at that
# grain, `syn` = the synthetic frame to fill). Returns `syn` with the new cols.
synth_sequence <- function(train, syn, cols, available, methods, control) {
  n_syn <- nrow(syn)
  for (v in cols) {
    y <- train[[v]]
    if (length(available) == 0L) {
      syn[[v]] <- sample_draw(y, n_syn, proper = control$proper)
    } else {
      syn[[v]] <- switch(
        methods[[v]],
        cart   = cart_draw(y, train[available], syn[available], control),
        sample = sample_draw(y, n_syn, proper = control$proper)
      )
    }
    available <- c(available, v)
  }
  syn
}

# Produce one synthetic data.frame. Subject-invariant columns are synthesised
# once per unit (subject grain) and broadcast across that unit's rows; remaining
# time-varying columns are synthesised at row grain, conditioning on the
# structural indices, the subject-level columns and any earlier time-varying
# columns.
synthesise_once <- function(data, st, subj_cols, time_cols, fixed_cols,
                            methods, control) {
  id_col <- st$id
  syn <- build_skeleton(data, st, control$k)   # id + structural indices

  if (length(subj_cols)) {
    sdat <- data[!duplicated(data[[id_col]]), c(id_col, subj_cols), drop = FALSE]
    ssyn <- data.frame(sort(unique(syn[[id_col]])))
    names(ssyn) <- id_col
    ssyn <- synth_sequence(sdat, ssyn, subj_cols, character(0), methods, control)
    pos <- match(syn[[id_col]], ssyn[[id_col]])
    for (v in subj_cols) syn[[v]] <- ssyn[[v]][pos]
  }

  syn <- synth_sequence(data, syn, time_cols, c(fixed_cols, subj_cols),
                        methods, control)

  syn <- syn[names(data)]          # restore original column order
  rownames(syn) <- NULL
  syn
}
