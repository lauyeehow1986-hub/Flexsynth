# Method registry: the single source of truth for per-variable synthesis
# methods. A method is a (fit, draw) pair plus a little metadata. Built-in
# methods are registered in .onLoad (see zzz.R); users add their own with
# register_method().
#
# fit(y, x, control)   -> a fitted model object (opaque to the engine)
# draw(model, x, n, control) -> a length-n vector of synthetic values
#
# `x` is a data.frame of predictor columns (possibly zero columns); `n` is the
# number of rows to synthesise.

# ---------------------------------------------------------------------------
# Registry storage + engine dispatch
# ---------------------------------------------------------------------------

.method_registry <- new.env(parent = emptyenv())

# Look up a method spec by name (used throughout the engine).
get_method <- function(name) {
  if (!exists(name, envir = .method_registry, inherits = FALSE)) {
    stop(sprintf("unsupported method '%s' (registered: %s).",
                 name, paste(list_methods(), collapse = ", ")), call. = FALSE)
  }
  get(name, envir = .method_registry, inherits = FALSE)
}

#' List the available synthesis methods
#'
#' Return the names of every method currently registered (built-in plus any
#' added with [register_method()]).
#'
#' @return A character vector of method names.
#' @export
#' @examples
#' list_methods()
list_methods <- function() sort(ls(envir = .method_registry))

# Build a fitted per-variable synthesiser. A predictor-hungry method with no
# predictors falls back to the marginal sampler; a type check keeps numeric-only
# / categorical-only methods from being applied to the wrong kind of variable.
fit_var <- function(y, xtrain, preds, method, control) {
  spec <- get_method(method)
  if (length(preds) == 0L && isTRUE(spec$needs_predictors)) {
    spec <- get_method("sample")
  }
  if (is.numeric(y) && !isTRUE(spec$numeric)) {
    stop(sprintf("method '%s' cannot synthesise the numeric variable.", method),
         call. = FALSE)
  }
  if (!is.numeric(y) && !isTRUE(spec$categorical)) {
    stop(sprintf("method '%s' cannot synthesise a categorical variable.", method),
         call. = FALSE)
  }
  list(spec = spec, model = spec$fit(y, xtrain, control))
}

# Draw synthetic values from a model built by fit_var().
apply_var <- function(model, xsyn, n, control) {
  model$spec$draw(model$model, xsyn, n, control)
}

#' Register a synthesis method
#'
#' Add a custom per-variable synthesiser to the method registry so it can be
#' selected by name via the `method` argument of [synth()] / [synth_linked()] or
#' the per-variable `method` in [synth_control()]. Registering a name that
#' already exists overwrites it, so you can override a built-in.
#'
#' A method is a pair of functions:
#' \describe{
#'   \item{`fit(y, x, control)`}{fits a model of the target `y` on the predictor
#'     `data.frame` `x` (which may have zero columns) and returns any object;}
#'   \item{`draw(model, x, n, control)`}{returns a length-`n` vector of synthetic
#'     values for new predictor rows `x`.}
#' }
#' `control` is the [synth_control()] object, so a method can honour `proper`,
#' read its own hyperparameters, and so on.
#'
#' @param name Method name (a single string) used to select it.
#' @param fit A function `function(y, x, control)` returning a fitted model.
#' @param draw A function `function(model, x, n, control)` returning `n` values.
#' @param numeric,categorical Logical flags declaring which target types the
#'   method supports (used for validation). Both default to `TRUE`.
#' @param needs_predictors Logical; if `TRUE` (default) the engine falls back to
#'   marginal sampling when the variable has no predictors.
#'
#' @return Invisibly, the method name.
#' @export
#' @examples
#' # A trivial "mean" method for numeric variables.
#' register_method(
#'   "constant_mean",
#'   fit  = function(y, x, control) mean(y),
#'   draw = function(model, x, n, control) rep(model, n),
#'   categorical = FALSE
#' )
#' "constant_mean" %in% list_methods()
register_method <- function(name, fit, draw, numeric = TRUE,
                            categorical = TRUE, needs_predictors = TRUE) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || name == "") {
    stop("`name` must be a single non-empty string.", call. = FALSE)
  }
  if (!is.function(fit) || !is.function(draw)) {
    stop("`fit` and `draw` must both be functions.", call. = FALSE)
  }
  assign(name, list(fit = fit, draw = draw,
                    numeric = isTRUE(numeric),
                    categorical = isTRUE(categorical),
                    needs_predictors = isTRUE(needs_predictors)),
         envir = .method_registry)
  invisible(name)
}

# ---------------------------------------------------------------------------
# Predictor selection + smoothing (synth_control wiring)
# ---------------------------------------------------------------------------

# Restrict `candidates` to the predictors allowed for `target` by a 0/1
# predictor matrix (rows = targets, cols = predictors). A NULL matrix or an
# absent target row imposes no restriction.
allowed_predictors <- function(pm, target, candidates) {
  if (is.null(pm) || !length(candidates)) return(candidates)
  if (is.null(rownames(pm)) || !target %in% rownames(pm)) return(candidates)
  keep <- colnames(pm)[pm[target, ] != 0]
  intersect(candidates, keep)
}

# Which of `cols` should have their numeric draws smoothed? `smoothing` may be
# NULL (none), TRUE / "density" (all numeric), or a character vector of names.
smoothing_targets <- function(control, cols) {
  sm <- control$smoothing
  if (is.null(sm)) return(character(0))
  if (isTRUE(sm) || identical(sm, "density")) return(cols)
  if (is.character(sm)) return(intersect(sm, cols))
  character(0)
}

# Kernel-smooth a numeric synthetic draw: add Gaussian noise with a rule-of-thumb
# bandwidth, then clamp to the observed range. Degenerate (few-valued) or
# non-numeric draws are returned unchanged; integer targets stay integer.
smooth_draw <- function(val, yref) {
  if (!is.numeric(val) || length(unique(val[!is.na(val)])) < 3L) return(val)
  bw <- stats::bw.nrd0(val)
  if (!is.finite(bw) || bw <= 0) return(val)
  out <- val + stats::rnorm(length(val), 0, bw)
  rng <- range(yref, na.rm = TRUE)
  out <- pmin(pmax(out, rng[1L]), rng[2L])
  if (is.integer(yref)) out <- as.integer(round(out))
  out
}

# ---------------------------------------------------------------------------
# Built-in methods
# ---------------------------------------------------------------------------

# sample: marginal empirical bootstrap (no predictors needed).
method_sample <- list(
  fit  = function(y, x, control) list(pool = y),
  draw = function(model, x, n, control)
    sample_draw(model$pool, n, proper = isTRUE(control$proper)),
  numeric = TRUE, categorical = TRUE, needs_predictors = FALSE
)

# cart: rpart leaf conditional bootstrap.
method_cart <- list(
  fit  = function(y, x, control) cart_fit(y, x, control),
  draw = function(model, x, n, control) cart_apply(model, x),
  numeric = TRUE, categorical = TRUE, needs_predictors = TRUE
)

# forest: a bagged ensemble of CART trees (a random forest built on rpart, so no
# extra dependency). Each tree is grown on a bootstrap sample of the rows and,
# when `mtry` is set, a random subset of the predictors; each synthetic row is
# routed through one randomly chosen tree and leaf-bootstrapped.
forest_fit <- function(y, x, control) {
  fp    <- control$forest
  ntree <- if (!is.null(fp$ntree)) max(1L, as.integer(fp$ntree)) else 10L
  allp  <- names(x)
  mtry  <- if (!is.null(fp$mtry)) min(as.integer(fp$mtry), length(allp)) else length(allp)
  n     <- length(y)
  trees <- lapply(seq_len(ntree), function(b) {
    bi   <- sample.int(n, n, replace = TRUE)
    cols <- if (length(allp) && mtry < length(allp)) sample(allp, mtry) else allp
    cart_fit(y[bi], x[bi, cols, drop = FALSE], control, bootstrap = FALSE)
  })
  list(trees = trees, template = y[integer(0)])
}
forest_draw <- function(model, x, n, control) {
  trees <- model$trees
  out   <- model$template[rep(NA_integer_, n)]        # typed NA vector, length n
  pick  <- sample.int(length(trees), n, replace = TRUE)
  for (b in sort(unique(pick))) {
    rows <- which(pick == b)
    tr   <- trees[[b]]
    out[rows] <- cart_apply(tr, x[rows, names(tr$xref), drop = FALSE])
  }
  out
}
method_forest <- list(fit = forest_fit, draw = forest_draw,
                      numeric = TRUE, categorical = TRUE, needs_predictors = TRUE)

# ctree: partykit conditional-inference tree, leaf-bootstrapped like cart.
ctree_fit <- function(y, x, control) {
  if (!requireNamespace("partykit", quietly = TRUE)) {
    stop("method = \"ctree\" needs the 'partykit' package. ",
         "install.packages(\"partykit\").", call. = FALSE)
  }
  x  <- prep_predictors(x)
  yy <- y
  if (isTRUE(control$proper)) {
    bi <- sample.int(length(y), length(y), replace = TRUE)
    yy <- y[bi]; x <- x[bi, , drop = FALSE]
  }
  # partykit::ctree needs a factor/numeric response; a bare character (or
  # logical) column trips ".y2infl: unknown response class". Coerce for the fit
  # only -- draws still sample the original values held in `yy`.
  resp <- if (is.character(yy) || is.logical(yy)) factor(yy) else yy
  df  <- data.frame(.y = resp, x, check.names = FALSE)
  fit <- partykit::ctree(stats::as.formula(".y ~ ."), data = df)
  node <- predict(fit, type = "node")
  list(fit = fit, y = yy, tr_leaf = as.integer(node), xref = x)
}
ctree_draw <- function(model, x, n, control) {
  x    <- prep_predictors(x, levels_from = model$xref)
  node <- as.integer(predict(model$fit, newdata = x, type = "node"))
  draw_from_leaves(model$y, model$tr_leaf, node)
}
method_ctree <- list(fit = ctree_fit, draw = ctree_draw,
                     numeric = TRUE, categorical = TRUE, needs_predictors = TRUE)

# --- parametric numeric methods (norm / normrank) --------------------------

# Fit a linear model of a numeric y on predictors, keeping what the draw needs
# to add either Gaussian residual noise (norm) or a rank map (normrank).
norm_fit <- function(y, x, control) {
  if (!is.numeric(y)) {
    stop("method = \"norm\" needs a numeric target.", call. = FALSE)
  }
  x  <- prep_predictors(x)
  df <- data.frame(.y = y, x, check.names = FALSE)
  mt <- if (ncol(x)) stats::terms(stats::as.formula(".y ~ ."), data = df)
        else          stats::terms(stats::as.formula(".y ~ 1"), data = df)
  X    <- stats::model.matrix(mt, df)
  qrX  <- qr(X)
  beta <- qr.coef(qrX, y)
  beta[is.na(beta)] <- 0
  resid  <- as.vector(y - X %*% beta)
  dfres  <- max(nrow(X) - qrX$rank, 1L)
  s2     <- sum(resid^2) / dfres
  list(mt = mt, beta = beta, s2 = s2, df = dfres,
       R = qr.R(qrX), full_rank = qrX$rank == ncol(X),
       xref = x, y = y, rng = range(y), as_int = is.integer(y))
}

# Continuous prediction mu + noise; `proper` also draws (sigma^2, beta) from the
# Bayesian posterior. Shared by norm (clamp / integerise) and normrank (rank map).
norm_continuous <- function(model, x, control) {
  x    <- prep_predictors(x, levels_from = model$xref)
  Xs   <- stats::model.matrix(model$mt, data.frame(.y = 0, x, check.names = FALSE))
  beta <- model$beta; s2 <- model$s2
  if (isTRUE(control$proper) && model$full_rank) {
    s2   <- model$s2 * model$df / stats::rchisq(1L, model$df)
    Rinv <- backsolve(model$R, diag(length(beta)))
    beta <- beta + sqrt(s2) * as.vector(Rinv %*% stats::rnorm(length(beta)))
  }
  as.vector(Xs %*% beta) + stats::rnorm(nrow(Xs), 0, sqrt(s2))
}

norm_draw <- function(model, x, n, control) {
  out <- norm_continuous(model, x, control)
  out <- pmin(pmax(out, model$rng[1L]), model$rng[2L])
  if (model$as_int) out <- as.integer(round(out))
  out
}
method_norm <- list(fit = norm_fit, draw = norm_draw,
                    numeric = TRUE, categorical = FALSE, needs_predictors = TRUE)

# normrank: predict on the continuous scale, then map each synthetic value to the
# observed value at the matching rank, so the marginal is preserved exactly.
normrank_draw <- function(model, x, n, control) {
  z   <- norm_continuous(model, x, control)
  obs <- sort(model$y)
  q   <- (rank(z, ties.method = "first") - 0.5) / length(z)
  idx <- pmax(1L, pmin(length(obs), ceiling(q * length(obs))))
  obs[idx]
}
method_normrank <- list(fit = norm_fit, draw = normrank_draw,
                        numeric = TRUE, categorical = FALSE, needs_predictors = TRUE)

# Register every built-in method. Called from .onLoad so the registry is ready
# before any synthesis runs.
register_builtin_methods <- function() {
  register_method("sample", method_sample$fit, method_sample$draw,
                  numeric = TRUE, categorical = TRUE, needs_predictors = FALSE)
  register_method("cart", method_cart$fit, method_cart$draw)
  register_method("forest", method_forest$fit, method_forest$draw)
  register_method("ctree", method_ctree$fit, method_ctree$draw)
  register_method("norm", method_norm$fit, method_norm$draw, categorical = FALSE)
  register_method("normrank", method_normrank$fit, method_normrank$draw,
                  categorical = FALSE)
  invisible(TRUE)
}
