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

# Calibrate the noise for a release of `n_marginals` marginals at sensitivity
# `cap`. Returns a list describing the mechanism and the per-cell noise, plus a
# closure `add_noise(counts)` that draws the calibrated noise for one histogram.
dp_calibrate <- function(dp, n_marginals, cap) {
  n_marginals <- as.integer(n_marginals)
  cap <- as.numeric(cap)
  eps <- dp$epsilon
  if (dp$mechanism == "laplace") {
    # Total L1 sensitivity of the concatenated release is n_marginals * cap.
    scale <- (n_marginals * cap) / eps
    add_noise <- function(counts)
      counts + rlaplace(length(counts), scale = scale)
    list(mechanism = "laplace", scale = scale, epsilon = eps, delta = 0,
         rho = NA_real_, n_marginals = n_marginals, cap = cap,
         add_noise = add_noise)
  } else {
    rho <- zcdp_rho_for(eps, dp$delta)
    sigma <- cap * sqrt(n_marginals / (2 * rho))
    add_noise <- function(counts)
      counts + stats::rnorm(length(counts), sd = sigma)
    list(mechanism = "gaussian", sigma = sigma, epsilon = eps, delta = dp$delta,
         rho = rho, n_marginals = n_marginals, cap = cap,
         add_noise = add_noise)
  }
}

# Laplace(0, scale) draws built from two exponentials (base R has no rlaplace).
rlaplace <- function(n, scale) {
  u <- stats::runif(n) - 0.5
  -scale * sign(u) * log1p(-2 * abs(u))
}

# Accounting record attached to a DP synth_result$privacy slot.
new_dp_accounting <- function(dp, calib, cap, n_marginals, variables, dropped) {
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
      rows_dropped = dropped
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
  cat("  model     :", x$dependence, "over", length(x$variables), "variables\n")
  cat("  marginals :", x$n_marginals,
      "(measured under composed budget)\n")
  cat("  noise     :",
      if (x$mechanism == "laplace") "Laplace scale" else "Gaussian sd",
      signif(x$noise, 4), "per cell\n")
  cat("  row cap   :", x$cap, "per", x$unit,
      if (x$rows_dropped > 0) paste0("(", x$rows_dropped, " rows dropped by capping)") else "",
      "\n")
  invisible(x)
}
