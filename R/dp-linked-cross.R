# Internal: cross-table conditioning for the Track B linked DP synthesiser.
# Not exported. Enabled by dp_control(cross_table = TRUE).
#
# Without cross-table conditioning a child table's variables are modelled by
# their own within-table marginals, so the synthetic child links to a synthetic
# parent (referential integrity) but is statistically independent of it. This
# module conditions each child variable on the *synthetic parent's* attributes.
#
# For a child table T with a modellable immediate parent P we measure, in
# addition to T's own one-way (and, for a tree, pairwise) marginals:
#   * one parent-by-child joint histogram for every (parent var u, child var v),
#     counted at the CHILD grain - one observation per child row, at the parent's
#     value carried down the foreign key. Its person-sensitivity is path_cap[T]
#     (an entity owns <= path_cap[T] child rows, each moving one cell), identical
#     to a child one-way marginal, so it folds into the same exact composition.
# The child's Chow-Liu structure is then learned over a graph whose nodes are the
# parent variables (fixed context, pre-placed in the spanning tree) plus the
# child variables: each child variable's single strongest predictor may be a
# parent variable or another child variable. At generation the synthetic parent's
# already-drawn value conditions the child draw, so cross-table dependence - not
# just the key link - survives the noise.
#
# Under dependence = "independent" the child-child edges are suppressed (weight
# -Inf), so each child variable conditions on its single best parent variable and
# nothing else. Deeper ancestors reach a child only through its immediate parent's
# synthesised values (as in Track A). If the parent has no modellable variable the
# caller falls back to the plain within-table model.

# Number of variable marginals a CHILD table contributes when conditioned on its
# parent: the within-table count (one-way, + pairwise for a tree) plus one
# parent-by-child joint for every (parent var, child var) pair.
dp_child_nvarmarg <- function(nC, nP, dp) {
  dp_n_var_marg(nC, dp) + nC * nP
}

# Integer joint-count matrix of two aligned code vectors (na x nb), NA-dropping.
dp_joint_count <- function(a, b, na, nb) {
  ok <- !is.na(a) & !is.na(b)
  tt <- table(factor(a[ok], levels = seq_len(na)),
              factor(b[ok], levels = seq_len(nb)))
  matrix(as.integer(tt), na, nb)
}

# Maximum-weight spanning tree with a set of nodes `seeded` already in the tree
# (fixed context). Grows only the remaining nodes, each attaching to its
# strongest already-reachable neighbour. Returns directed c(parent, child) edges
# in add order (each edge's parent is in the tree when the edge is used), so it is
# a valid conditional sampling order given the seeded nodes' values.
dp_mst_seeded <- function(W, seeded) {
  d <- nrow(W)
  in_tree <- logical(d)
  in_tree[seeded] <- TRUE
  if (all(in_tree) || !length(seeded)) return(list())
  best <- rep(-Inf, d)
  best_from <- rep(NA_integer_, d)
  for (s in seeded) {
    upd <- !in_tree & (W[s, ] > best)
    best[upd] <- W[s, upd]
    best_from[upd] <- s
  }
  n_add <- sum(!in_tree)
  edges <- vector("list", n_add)
  for (step in seq_len(n_add)) {
    cand <- which(!in_tree)
    j <- cand[which.max(best[cand])]
    in_tree[j] <- TRUE
    edges[[step]] <- c(best_from[j], j)
    upd <- !in_tree & (W[j, ] > best)
    best[upd] <- W[j, upd]
    best_from[upd] <- j
  }
  edges
}

# Parent context codes at the CHILD grain: for every real child row of `t`, the
# code (in the parent's own domain grid) of each parent variable, looked up along
# the foreign key. Returned as a named list (per parent var) of aligned code
# vectors. Orphans (no matching parent, should not occur after capping) map to the
# first cell.
dp_parent_ctx_codes <- function(cdata, hierarchy, t, p, parent_dom, parent_vars) {
  fk <- hierarchy$fk[[t]]
  pk <- key_string(cdata[[p]], fk)
  ck <- key_string(cdata[[t]], fk)
  mi <- match(ck, pk)                                  # parent row per child row
  stats::setNames(lapply(parent_vars, function(u) {
    pc <- dp_encode(parent_dom[[u]], cdata[[p]][[u]])  # code per parent row
    code <- pc[mi]
    code[is.na(code)] <- 1L
    as.integer(code)
  }), parent_vars)
}

# Fit a child variable model conditioned on parent context. `child_codes` /
# `child_nbins` are the child's own discretised variables; `parent_codes` /
# `parent_nbins` the parent's variables at the child grain (from
# dp_parent_ctx_codes). Measures child one-way marginals, parent-by-child joints
# (always) and child-child joints (tree only) under `calib`, then builds a
# seeded max-spanning tree (parent nodes fixed) and per-edge conditionals.
dp_fit_child_cross <- function(child_codes, child_nbins,
                               parent_codes, parent_nbins, dp, calib) {
  cvars <- names(child_codes)
  pvars <- names(parent_codes)
  nC <- length(cvars)
  nP <- length(pvars)
  tree <- dp$dependence == "tree" && nC > 1L

  # Child one-way (clipped nonneg), reused as root/fallback marginals.
  c1 <- stats::setNames(lapply(cvars, function(v)
    pmax(calib$add_noise(tabulate(child_codes[[v]], child_nbins[[v]])), 0)), cvars)
  n_est <- max(1, round(mean(vapply(c1, sum, numeric(1)))))

  # Node layout: 1..nP parent context, nP+1..nP+nC child. Joints stored at the
  # lower-index slot with rows = lower-index variable's codes.
  nNode <- nP + nC
  W <- matrix(-Inf, nNode, nNode)
  joints <- matrix(list(), nNode, nNode)

  # Parent-by-child joints (rows = parent code, cols = child code).
  for (ci in seq_len(nC)) {
    cnode <- nP + ci
    cv <- cvars[ci]
    for (pj in seq_len(nP)) {
      pv <- pvars[pj]
      tab <- dp_joint_count(parent_codes[[pv]], child_codes[[cv]],
                            parent_nbins[[pv]], child_nbins[[cv]])
      noisy <- pmax(matrix(calib$add_noise(as.vector(tab)),
                           parent_nbins[[pv]], child_nbins[[cv]]), 0)
      joints[[pj, cnode]] <- noisy
      W[pj, cnode] <- W[cnode, pj] <- dp_mutual_information(noisy)
    }
  }
  # Child-child joints (tree only; rows = earlier child var's codes).
  if (tree) {
    for (i in seq_len(nC - 1L)) for (j in (i + 1L):nC) {
      ni <- nP + i; nj <- nP + j
      tab <- dp_joint_count(child_codes[[cvars[i]]], child_codes[[cvars[j]]],
                            child_nbins[[cvars[i]]], child_nbins[[cvars[j]]])
      noisy <- pmax(matrix(calib$add_noise(as.vector(tab)),
                           child_nbins[[cvars[i]]], child_nbins[[cvars[j]]]), 0)
      joints[[ni, nj]] <- noisy
      W[ni, nj] <- W[nj, ni] <- dp_mutual_information(noisy)
    }
  }

  edges <- dp_mst_seeded(W, seq_len(nP))
  cond <- lapply(edges, function(e) {
    a <- e[1L]; b <- e[2L]                              # b is always a child node
    jt <- if (a < b) joints[[a, b]] else t(joints[[b, a]])  # orient rows = a
    t(apply(jt, 1L, dp_normalise))                     # P(b | a)
  })

  # Child-by-child mutual information (post-processing of the already-measured
  # child-child joints), dimnamed by child variable. When this model seeds the
  # initial state of a longitudinally-modelled child (combined cross + longi), the
  # transition tensors pick cross-variable parents from it exactly as the flat
  # engine reads dp_fit_model()$pairwise_mi. Only populated for a tree (the child-
  # child joints are unmeasured under dependence = "independent").
  child_mi <- if (tree) {
    cc <- W[nP + seq_len(nC), nP + seq_len(nC), drop = FALSE]
    dimnames(cc) <- list(cvars, cvars)
    cc
  } else NULL

  list(kind = "child-cross", cvars = cvars, vars = cvars, pvars = pvars, nP = nP,
       child_nbins = child_nbins, edges = edges, cond = cond,
       marginals = lapply(c1, dp_normalise), n_est = n_est,
       pairwise_mi = child_mi)
}

# Draw child variable codes conditioned on parent context. `parent_context` is a
# named list (per parent var) of code vectors aligned to the `n` child rows (from
# the synthetic parent). Returns an n x nC integer code matrix (cols = cvars).
dp_sample_child_codes <- function(model, parent_context, n) {
  nNode <- model$nP + length(model$cvars)
  node_codes <- vector("list", nNode)
  for (pj in seq_len(model$nP)) {
    node_codes[[pj]] <- parent_context[[model$pvars[pj]]]
  }
  # Any child never reached by an edge (only if it has no context and no picked
  # neighbour) falls back to its one-way marginal.
  reached <- logical(nNode)
  reached[seq_len(model$nP)] <- TRUE
  for (k in seq_along(model$edges)) {
    e <- model$edges[[k]]; a <- e[1L]; b <- e[2L]
    cond <- model$cond[[k]]
    av <- node_codes[[a]]
    child <- integer(n)
    for (val in unique(av)) {
      rows <- which(av == val)
      child[rows] <- dp_sample_cat(length(rows), cond[val, ])
    }
    node_codes[[b]] <- child
    reached[b] <- TRUE
  }
  out <- matrix(NA_integer_, n, length(model$cvars),
                dimnames = list(NULL, model$cvars))
  for (ci in seq_along(model$cvars)) {
    cnode <- model$nP + ci
    out[, ci] <- if (reached[cnode]) node_codes[[cnode]]
                 else dp_sample_cat(n, model$marginals[[ci]])
  }
  out
}
