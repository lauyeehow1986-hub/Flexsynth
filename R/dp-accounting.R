# Internal: differential-privacy noise calibration and budget accounting.
# Not exported. All composition is exact for the mechanisms we implement.
#
# We release a fixed set of `n_marginals` marginal histograms, each measured on
# the (contribution-bounded) data. Under person-level privacy with a per-person
# row cap of `cap`, adding or removing one person changes any single marginal by
# at most `cap` in L1 (total mass moved) and at most `cap` in L2 (worst case all
# `cap` rows fall in one cell). So each marginal query has sensitivity
# Delta1 = Delta2 = cap.
#
#   * Laplace (pure eps-DP): the marginals are released together, so their
#     L1 sensitivities add. Total Delta1 = n_marginals * cap, and independent
#     Laplace(scale = Delta1 / eps) noise on every cell gives eps-DP.
#   * Gaussian ((eps, delta)-DP via zero-concentrated DP): each Gaussian
#     release with noise sd = sigma satisfies rho_i = Delta2^2 / (2 sigma^2)
#     zCDP; they compose additively to rho = n_marginals * cap^2 / (2 sigma^2).
#     rho-zCDP implies (eps, delta)-DP with
#         eps = rho + 2 * sqrt(rho * ln(1/delta)).
#     We invert that to spend exactly (eps, delta): with L = ln(1/delta),
#         sqrt(rho) = sqrt(L + eps) - sqrt(L),  rho = (sqrt(L + eps) - sqrt(L))^2,
#     then sigma = cap * sqrt(n_marginals / (2 rho)).

# Total zCDP budget rho that is (eps, delta)-DP-equivalent (delta > 0).
zcdp_rho_for <- function(epsilon, delta) {
  L <- log(1 / delta)
  (sqrt(L + epsilon) - sqrt(L))^2
}

# Core noise calibration from explicit release sensitivities. A release is a
# concatenation of histograms; `total_l1` is the total L1 person-sensitivity
# (used by Laplace) and `sum_sq` is the sum of squared L2 person-sensitivities
# over the histograms (used by the Gaussian / zCDP path). `budget_frac` (in
# (0, 1]) is the share of the total budget these measurements may spend; the
# remainder is reserved for DP domain estimation (see `dp_quantile_eps()`). The
# split is exact: pure eps adds and zCDP rho adds, so measurements + domain never
# exceed (eps, delta). Returns the mechanism description plus a closure
# `add_noise(counts)` that draws one histogram's worth of calibrated noise.
#
# Uniform per-cell noise with these totals is correct for heterogeneous
# sensitivities: for Laplace, Laplace(total_l1 / eps) on every cell gives eps-DP
# because the concatenated L1 sensitivity is total_l1; for Gaussian, sd
# sqrt(sum_sq / (2 rho)) gives rho-zCDP because the summed per-release zCDP is
# sum_sq / (2 sigma^2).
dp_make_noise <- function(dp, total_l1, sum_sq, budget_frac = 1) {
  eps <- dp$epsilon
  if (dp$mechanism == "laplace") {
    eps_marg <- budget_frac * eps          # measurements' share of pure-eps budget
    scale <- total_l1 / eps_marg
    add_noise <- function(counts)
      counts + rlaplace(length(counts), scale = scale)
    list(mechanism = "laplace", scale = scale, epsilon = eps, delta = 0,
         rho = NA_real_, budget_frac = budget_frac, add_noise = add_noise)
  } else {
    rho_total <- zcdp_rho_for(eps, dp$delta)
    rho_marg <- budget_frac * rho_total     # measurements' share of zCDP rho
    sigma <- sqrt(sum_sq / (2 * rho_marg))
    add_noise <- function(counts)
      counts + stats::rnorm(length(counts), sd = sigma)
    list(mechanism = "gaussian", sigma = sigma, epsilon = eps, delta = dp$delta,
         rho = rho_total, rho_marginals = rho_marg, budget_frac = budget_frac,
         add_noise = add_noise)
  }
}

# Calibrate the noise for a flat release of `n_marginals` marginals, each at
# person-sensitivity `cap` (adding / removing one person moves at most `cap`
# counts in any single marginal). Total L1 = n_marginals * cap, and the summed
# squared L2 is n_marginals * cap^2, so this is a thin wrapper over
# `dp_make_noise()` that reproduces the original flat calibration exactly.
dp_calibrate <- function(dp, n_marginals, cap, budget_frac = 1) {
  n_marginals <- as.integer(n_marginals)
  cap <- as.numeric(cap)
  core <- dp_make_noise(dp, total_l1 = n_marginals * cap,
                        sum_sq = n_marginals * cap^2, budget_frac = budget_frac)
  c(core, list(n_marginals = n_marginals, cap = cap))
}

# Per-query epsilon for the exponential-mechanism quantiles used to estimate bin
# edges under DP, given `n_queries` such queries share a `budget_frac` slice of
# the budget. For Laplace (pure eps) the slice is split evenly and adds. For
# Gaussian we work in zCDP: allocate `budget_frac` of rho to the domain, split it
# evenly, and invert the conservative pure-eps -> (eps^2 / 2)-zCDP bound
# (eps_q = sqrt(2 * rho_per_query)); this never under-charges the true
# (eps_q^2 / 8)-zCDP cost of a bounded-range exponential mechanism.
dp_quantile_eps <- function(dp, n_queries, budget_frac) {
  n_queries <- as.integer(n_queries)
  if (dp$mechanism == "laplace") {
    (budget_frac * dp$epsilon) / n_queries
  } else {
    rho_dom <- budget_frac * zcdp_rho_for(dp$epsilon, dp$delta)
    sqrt(2 * (rho_dom / n_queries))
  }
}

# Laplace(0, scale) draws built from two exponentials (base R has no rlaplace).
rlaplace <- function(n, scale) {
  u <- stats::runif(n) - 0.5
  -scale * sign(u) * log1p(-2 * abs(u))
}

# Accounting record attached to a DP synth_result$privacy slot. `domain`
# describes how bin edges were chosen and, under DP estimation, what it cost:
# a list with `mode`, the estimated variables, the per-query epsilon, and the
# budget fraction reserved for it (0 when nothing was estimated). `longitudinal`
# is NULL for a flat release, or a list describing the DP Markov model (the
# number of initial-state marginals, the number of transition matrices, and the
# per-person row cap that bounds the transition sensitivity) when within-unit
# temporal structure was modelled.
new_dp_accounting <- function(dp, calib, cap, n_marginals, variables, dropped,
                              domain = list(mode = dp$domain, vars = character(0),
                                            eps_per_query = NA_real_, frac = 0),
                              longitudinal = NULL, linked = NULL, learn = NULL,
                              adaptive = NULL, bayes = NULL, aim = NULL,
                              estimator = "local") {
  structure(
    list(
      epsilon = dp$epsilon,
      delta = dp$delta,
      unit = dp$unit,
      mechanism = calib$mechanism,
      dependence = dp$dependence,
      cap = cap,
      n_marginals = n_marginals,
      noise = if (calib$mechanism == "laplace") calib$scale else calib$sigma,
      rho = calib$rho,
      variables = variables,
      rows_dropped = dropped,
      domain = domain,
      longitudinal = longitudinal,
      linked = linked,
      learn = learn,
      adaptive = adaptive,
      bayes = bayes,
      aim = aim,
      estimator = estimator
    ),
    class = "dp_accounting"
  )
}

#' @export
print.dp_accounting <- function(x, ...) {
  cat("<dp_accounting> Track B differential privacy\n")
  cat("  guarantee : (epsilon =", x$epsilon,
      if (x$delta > 0) paste0(", delta = ", format(x$delta, scientific = TRUE)) else ", pure eps",
      ")-DP,", x$unit, "level\n")
  cat("  mechanism :", x$mechanism,
      if (x$mechanism == "gaussian") paste0("(rho = ", signif(x$rho, 3), " zCDP)") else "", "\n")
  model_kind <- if (!is.null(x$aim))
    paste0("Full AIM graphical model (treewidth ", x$aim$treewidth,
           ", Private-PGM)")
    else if (!is.null(x$adaptive))
      paste0("adaptive junction tree (treewidth ", x$adaptive$treewidth, ")")
    else if (!is.null(x$bayes))
      paste0("degree ", x$bayes$degree, " Bayesian network (GreedyBayes)")
    else x$dependence
  cat("  model     :", model_kind, "over", length(x$variables), "variables",
      if (!is.null(x$longitudinal)) "(DP Markov: initial-state + transitions)"
      else if (!is.null(x$linked)) paste0("(linked: ", x$linked$n_tables, " tables)")
      else "",
      "\n")
  if (!is.null(x$linked)) {
    cat("  histograms:", x$n_marginals,
        "(per-table variable marginals + child count models, composed budget)\n")
  } else if (is.null(x$longitudinal)) {
    if (!is.null(x$aim)) {
      am <- x$aim
      if (isTRUE(am$anneal)) {
        cat("  marginals :", x$n_marginals,
            paste0("(", length(x$variables), " one-way + ", am$n_new,
                   " loopy + ", am$n_refine,
                   " refinement pair(s), composed budget)\n"))
      } else {
        cat("  marginals :", x$n_marginals,
            paste0("(", length(x$variables), " one-way + ", am$n_rounds,
                   " adaptively selected pairs, loopy, composed budget)\n"))
      }
    } else if (!is.null(x$adaptive)) {
      ad <- x$adaptive
      if (isTRUE(ad$anneal)) {
        cat("  marginals :", x$n_marginals,
            paste0("(", length(x$variables), " one-way + ", ad$n_cliques,
                   " spanning + ", ad$n_refine,
                   " refinement clique(s), composed budget)\n"))
      } else {
        cat("  marginals :", x$n_marginals,
            paste0("(", length(x$variables), " one-way + ", ad$n_cliques,
                   " adaptively selected cliques, composed budget)\n"))
      }
    } else if (!is.null(x$bayes)) {
      by <- x$bayes
      cat("  marginals :", x$n_marginals,
          paste0("(", by$n_nodes, " one-way + ", by$n_families,
                 " family joints, composed budget)\n"))
    } else if (!is.null(x$learn)) {
      lr <- x$learn
      cat("  marginals :", x$n_marginals,
          paste0("(", lr$n_struct, " pairwise scans + ", lr$n_param,
                 " parameter marginals, composed budget)\n"))
    } else {
      cat("  marginals :", x$n_marginals,
          "(measured under composed budget)\n")
    }
  } else {
    lg <- x$longitudinal
    cat("  histograms:", x$n_marginals,
        paste0("(1 length + ", lg$n_init_marg, " initial + ",
               lg$n_transitions, " transition, composed budget)\n"))
    if (!is.null(lg$baseline) && length(lg$baseline))
      cat("  baseline  : held constant within unit:",
          paste(lg$baseline, collapse = ", "), "\n")
    ord <- if (is.null(lg$order)) 1L else lg$order
    crs <- if (is.null(lg$cross)) 0L else lg$cross
    if (ord > 1L || crs > 0L) {
      cat("  transitions: order ", ord, " + ", crs, " cross-parent(s)",
          " (sensitivity cap - order)\n", sep = "")
      if (!is.null(lg$cross_parents)) {
        pairs <- vapply(names(lg$cross_parents), function(v) {
          cp <- lg$cross_parents[[v]]
          if (length(cp)) paste0(v, " ~ ", paste(cp, collapse = "+")) else NA_character_
        }, character(1))
        pairs <- pairs[!is.na(pairs)]
        if (length(pairs))
          cat("               cross-parents:", paste(pairs, collapse = "; "), "\n")
      }
    }
  }
  anneal_rec <- if (!is.null(x$adaptive) && isTRUE(x$adaptive$anneal)) x$adaptive
                else if (!is.null(x$aim) && isTRUE(x$aim$anneal)) x$aim
                else NULL
  if (!is.null(anneal_rec)) {
    cat("  noise     :",
        if (x$mechanism == "laplace") "Laplace scale" else "Gaussian sd",
        signif(anneal_rec$noise_min, 4), "-", signif(anneal_rec$noise_max, 4),
        "per cell (annealed range)\n")
  } else {
    cat("  noise     :",
        if (x$mechanism == "laplace") "Laplace scale" else "Gaussian sd",
        signif(x$noise, 4), "per cell",
        if (!is.null(x$learn)) "(parameters)" else "", "\n")
  }
  if (!is.null(x$learn)) {
    lr <- x$learn
    cat("  structure : budget-efficient,", signif(lr$frac, 3),
        "of budget selects the tree;",
        if (x$mechanism == "laplace") "scan Laplace scale" else "scan Gaussian sd",
        signif(lr$struct_noise, 4), "per cell\n")
  }
  if (!is.null(x$adaptive)) {
    ad <- x$adaptive
    if (isTRUE(ad$anneal)) {
      cat("  selection : annealed (AIM-style),", signif(ad$select_frac, 3),
          "of budget;", ad$n_rounds, "round(s) =", ad$n_cliques, "spanning +",
          ad$n_refine, "refinement;",
          if (x$mechanism == "laplace") "noise halved" else "sigma halved",
          paste0(ad$n_anneal_steps, "x"), "\n")
    } else {
      cat("  selection : adaptive (AIM-style),", signif(ad$select_frac, 3),
          "of budget over", ad$n_cliques, "exponential-mechanism round(s);",
          "per-round eps", signif(ad$select_eps, 4), "\n")
    }
  }
  if (!is.null(x$bayes)) {
    by <- x$bayes
    cat("  selection : GreedyBayes,", signif(by$select_frac, 3),
        "of budget over", by$n_select, "exponential-mechanism round(s);",
        "per-round eps", signif(by$select_eps, 4), "\n")
  }
  if (!is.null(x$aim)) {
    am <- x$aim
    if (isTRUE(am$anneal)) {
      cat("  selection : annealed Full AIM (loopy marginals),",
          signif(am$select_frac, 3), "of budget;", am$n_rounds,
          "round(s) =", am$n_new, "loopy +", am$n_refine, "refinement;",
          if (x$mechanism == "laplace") "noise halved" else "sigma halved",
          paste0(am$n_anneal_steps, "x"), "\n")
    } else {
      cat("  selection : Full AIM (loopy marginals),", signif(am$select_frac, 3),
          "of budget over", am$n_rounds, "exponential-mechanism round(s);",
          "per-round eps", signif(am$select_eps, 4), "\n")
    }
    cat("  estimator : Private-PGM reconciliation over the triangulated model\n")
    cat("              (belief propagation + mirror descent; post-processing,",
        "no extra budget)\n")
    if (identical(am$scoring, "model"))
      cat("  scoring   : model-projection (candidates scored against the",
          "reconciled model; post-processing, default)\n")
    else
      cat("  scoring   : independence (one-way product;",
          "model-projection disabled)\n")
  }
  if (identical(x$estimator, "pgm")) {
    cat("  estimator : Private-PGM reconciliation of the measured marginals\n")
    cat("              (belief propagation + mirror descent; post-processing,",
        "no extra budget)\n")
  }
  if (!is.null(x$linked)) {
    for (ti in x$linked$tables) {
      role <- if (ti$role == "root") "root, cap 1"
              else paste0("child of ", ti$parent, ", <=", ti$local_cap,
                          "/parent, path cap ", ti$path_cap,
                          if (isTRUE(ti$cross)) ", cond. on parent" else "",
                          if (isTRUE(ti$longitudinal)) {
                            if (isTRUE(ti$cross_init))
                              ", DP Markov over rows (initial state cond. on parent)"
                            else ", DP Markov over rows"
                          } else "")
      cat(sprintf("    - %-14s %s%s\n", ti$name, role,
                  if (ti$rows_dropped > 0)
                    paste0(" (", ti$rows_dropped, " rows dropped)") else ""))
      # Per-child within-unit transition controls (longitudinal children only).
      if (isTRUE(ti$longitudinal)) {
        if (!is.null(ti$baseline) && length(ti$baseline))
          cat("        baseline held:", paste(ti$baseline, collapse = ", "), "\n")
        to <- if (is.null(ti$tran_order)) 1L else ti$tran_order
        tc <- if (is.null(ti$tran_cross)) 0L else ti$tran_cross
        tp <- if (is.null(ti$tran_parent)) 0L else ti$tran_parent
        if (to > 1L || tc > 0L || tp > 0L) {
          cat("        transitions: order ", to, " + ", tc, " cross-parent(s)",
              if (tp > 0L) paste0(" + ", tp, " parent-attr(s)") else "",
              " (sensitivity cap - order)\n", sep = "")
          if (!is.null(ti$tran_cross_parents)) {
            pairs <- vapply(names(ti$tran_cross_parents), function(v) {
              cp <- ti$tran_cross_parents[[v]]
              if (length(cp)) paste0(v, " ~ ", paste(cp, collapse = "+"))
              else NA_character_
            }, character(1))
            pairs <- pairs[!is.na(pairs)]
            if (length(pairs))
              cat("        cross-parents:", paste(pairs, collapse = "; "), "\n")
          }
          if (!is.null(ti$tran_parent_parents)) {
            pairs <- vapply(names(ti$tran_parent_parents), function(v) {
              pp <- ti$tran_parent_parents[[v]]
              if (length(pp)) paste0(v, " ~ ", paste(pp, collapse = "+"))
              else NA_character_
            }, character(1))
            pairs <- pairs[!is.na(pairs)]
            if (length(pairs))
              cat("        parent-attrs:", paste(pairs, collapse = "; "), "\n")
          }
        }
      }
    }
  } else {
    cat("  row cap   :", x$cap, "per", x$unit,
        if (x$rows_dropped > 0) paste0("(", x$rows_dropped, " rows dropped by capping)") else "",
        "\n")
  }
  dm <- x$domain
  if (!is.null(dm)) {
    if (length(dm$vars) > 0) {
      cat("  bin edges :", dm$mode, "- estimated under DP for",
          paste(dm$vars, collapse = ", "),
          paste0("(", signif(dm$frac, 3), " of budget)"), "\n")
    } else if (identical(dm$mode, "data")) {
      cat("  bin edges :", dm$mode,
          "(from data range; NOT in the accounting)\n")
    } else {
      cat("  bin edges :", dm$mode, "(public bounds; no budget spent)\n")
    }
    ct <- dm$categorical
    if (!is.null(ct)) {
      cat("  categories: DP set-union discovered",
          paste(paste0(ct$vars, " (", ct$n_kept, " kept)"), collapse = ", "),
          paste0("(per-op eps ", signif(ct$eps_op, 3), ", delta ",
                 format(ct$delta_cat, scientific = TRUE), ")\n"))
    }
  }
  invisible(x)
}
