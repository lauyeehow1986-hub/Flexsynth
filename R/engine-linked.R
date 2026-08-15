# Internal: Phase 3 linked multi-table joint synthesis (Track A). Not exported.
#
# Tables are synthesised parent-first. A root table (one whose key has no parent
# among the supplied tables) is synthesised with the single-table engine. Each
# child table is then generated from its already-synthesised parent so that:
#   * every child foreign key equals a real synthetic parent key (referential
#     integrity is guaranteed by construction);
#   * the number of child rows per parent is drawn from a learned count
#     distribution that includes parents with *no* children (zero-inflation);
#   * the child's own structural index is regenerated 1..count within each unit;
#   * child variables are conditioned on the parent's synthesised attributes
#     (cross-table predictors), the child's own structural index, and earlier
#     child variables (immediate-parent conditioning — deeper ancestors reach a
#     child only through its parent's synthesised values).

# ---------------------------------------------------------------------------
# Hierarchy
# ---------------------------------------------------------------------------

# Row-wise composite key string for a set of columns.
key_string <- function(df, cols) {
  if (length(cols) == 1L) return(as.character(df[[cols]]))
  do.call(paste, c(lapply(cols, function(cc) df[[cc]]), sep = "\r"))
}

# Derive the parent map from the key lists alone. Table P is the parent of table
# T when P's full key equals T's key with its last (own) column dropped.
derive_parent_map <- function(keys) {
  tn  <- names(keys)
  sig <- vapply(tn, function(t) paste(keys[[t]], collapse = "\r"), character(1))
  by_sig <- stats::setNames(tn, sig)
  parent <- stats::setNames(rep(NA_character_, length(tn)), tn)
  fk  <- stats::setNames(vector("list", length(tn)), tn)
  own <- stats::setNames(rep(NA_character_, length(tn)), tn)
  for (t in tn) {
    k <- keys[[t]]
    if (length(k) >= 2L) {
      psig <- paste(k[-length(k)], collapse = "\r")
      if (psig %in% names(by_sig)) {
        parent[[t]] <- by_sig[[psig]]
        fk[[t]]  <- k[-length(k)]
        own[[t]] <- k[length(k)]
      }
    }
  }
  list(parent = parent, fk = fk, own = own)
}

# Topological order: every parent precedes its children.
topo_order <- function(tn, parent) {
  done <- character(0)
  remaining <- tn
  while (length(remaining)) {
    ready <- remaining[vapply(remaining, function(t) {
      p <- parent[[t]]
      is.na(p) || p %in% done
    }, logical(1))]
    if (!length(ready)) {
      stop("table hierarchy is cyclic or references a missing parent.",
           call. = FALSE)
    }
    done <- c(done, ready)
    remaining <- setdiff(remaining, ready)
  }
  done
}

# Parse + validate structures and keys, then build the full hierarchy.
link_hierarchy <- function(tables, structures, keys) {
  tn  <- names(tables)
  st  <- stats::setNames(vector("list", length(tn)), tn)
  key <- stats::setNames(vector("list", length(tn)), tn)
  for (t in tn) {
    st[[t]] <- parse_structure(structures[[t]], tables[[t]])
    k <- keys[[t]]
    if (!is.character(k) || !length(k)) {
      stop(sprintf("keys[['%s']] must be a non-empty character vector.", t),
           call. = FALSE)
    }
    miss <- setdiff(k, names(tables[[t]]))
    if (length(miss)) {
      stop(sprintf("keys[['%s']] names column(s) not in that table: %s",
                   t, paste(miss, collapse = ", ")), call. = FALSE)
    }
    if (anyDuplicated(key_string(tables[[t]], k))) {
      stop(sprintf(paste0("table '%s' has duplicate rows for its key (%s); ",
                          "keys must identify rows uniquely."),
                   t, paste(k, collapse = ", ")), call. = FALSE)
    }
    key[[t]] <- k
  }
  sig <- vapply(tn, function(t) paste(key[[t]], collapse = "\r"), character(1))
  if (anyDuplicated(sig)) {
    stop("two tables share identical key columns; cannot infer the hierarchy.",
         call. = FALSE)
  }
  pm  <- derive_parent_map(key)
  ord <- topo_order(tn, pm$parent)
  list(names = tn, structure = st, keys = key,
       parent = pm$parent, fk = pm$fk, own = pm$own, order = ord)
}

# ---------------------------------------------------------------------------
# Root + child synthesis
# ---------------------------------------------------------------------------

# Synthesise a root table with the single-table engine. `st` is the parsed
# structure for this table.
synth_root_table <- function(data, st, method, control) {
  fixed_cols <- st$nested
  synth_cols <- setdiff(names(data), st$vars)
  vs <- control$visit_sequence
  if (!is.null(vs)) {
    synth_cols <- c(intersect(vs, synth_cols), setdiff(synth_cols, vs))
  }
  methods   <- resolve_methods(method, control$method, synth_cols)
  subj_cols <- subject_level_cols(data, st$id, synth_cols)
  time_cols <- setdiff(synth_cols, subj_cols)
  synthesise_once(data, st, subj_cols, time_cols, fixed_cols, methods, control)
}

# Synthesise a child table conditionally on its already-synthesised parent.
synth_child_table <- function(child_real, parent_real, parent_syn,
                              fk_cols, own_col, keys_child, method, control) {
  parent_attr <- setdiff(names(parent_syn), fk_cols)   # cross-table predictors

  ## Count model: child rows per real parent unit, including zeros.
  pk_real <- key_string(parent_real, fk_cols)          # one per parent row (unique)
  ck_real <- key_string(child_real,  fk_cols)          # one per child row
  ctab <- table(ck_real)
  size_pool <- as.integer(ctab[match(pk_real, names(ctab))])
  size_pool[is.na(size_pool)] <- 0L

  ## Synthetic parent units and their drawn child counts.
  sk <- key_string(parent_syn, fk_cols)
  unit_first <- which(!duplicated(sk))
  n_units <- length(unit_first)
  draw <- size_pool[sample.int(length(size_pool), n_units, replace = TRUE)]

  keep <- draw > 0L
  if (!any(keep)) {                                    # no child rows at all
    out <- child_real[0, names(child_real), drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  ## Child skeleton: copy each parent unit's foreign key `count` times (this is
  ## what guarantees referential integrity), regenerate the own index 1..count,
  ## and attach the parent's synthesised attributes as carried predictors.
  rep_rows <- rep(unit_first[keep], draw[keep])
  skel <- parent_syn[rep_rows, fk_cols, drop = FALSE]
  skel[[own_col]]  <- sequence(draw[keep])             # 1,2,..,count within unit
  skel[[".unit"]]  <- rep(which(keep), draw[keep])     # integer unit id
  for (a in parent_attr) skel[[a]] <- parent_syn[[a]][rep_rows]
  rownames(skel) <- NULL

  ## Training frame: real child rows with the real parent attributes joined on.
  train <- child_real
  pidx <- match(ck_real, pk_real)
  for (a in parent_attr) train[[a]] <- parent_real[[a]][pidx]
  train[[".unit"]] <- match(ck_real, unique(ck_real))

  ## Synthesise the child's non-key variables, conditioning on parent attributes
  ## (subject grain) plus the own index and lags (row grain).
  child_vars <- setdiff(names(child_real), keys_child)
  if (length(child_vars)) {
    st_child <- list(id = ".unit", vars = c(".unit", own_col), nested = own_col)
    methods   <- resolve_methods(method, control$method, child_vars)
    subj_cols <- subject_level_cols(train, ".unit", child_vars)
    time_cols <- setdiff(child_vars, subj_cols)
    skel <- fill_columns(train, skel, st_child, subj_cols, time_cols,
                         subj_fixed = parent_attr,
                         time_fixed = c(own_col, parent_attr),
                         methods, control)
  }

  out <- skel[names(child_real)]                       # drop helpers, restore order
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# One synthetic collection
# ---------------------------------------------------------------------------

synth_linked_once <- function(tables, hierarchy, method, control) {
  syn <- stats::setNames(vector("list", length(hierarchy$names)), hierarchy$names)
  for (t in hierarchy$order) {
    p <- hierarchy$parent[[t]]
    if (is.na(p)) {
      syn[[t]] <- synth_root_table(tables[[t]], hierarchy$structure[[t]],
                                   method, control)
    } else {
      syn[[t]] <- synth_child_table(
        child_real  = tables[[t]],
        parent_real = tables[[p]],
        parent_syn  = syn[[p]],
        fk_cols     = hierarchy$fk[[t]],
        own_col     = hierarchy$own[[t]],
        keys_child  = hierarchy$keys[[t]],
        method      = method,
        control     = control)
    }
  }
  syn
}
