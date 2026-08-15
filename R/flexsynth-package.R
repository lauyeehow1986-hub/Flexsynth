#' flexsynth: flexible synthetic data for nested, longitudinal and linked data
#'
#' Generate high-quality synthetic data from real datasets of any structure,
#' natively in long format, with first-class support for nested / longitudinal
#' and multi-table *linked* data.
#'
#' @section Two privacy tracks:
#' - **Track A (default):** high-utility sequential synthesis. No formal
#'   guarantee; ships empirical disclosure-risk diagnostics.
#' - **Track B (opt-in):** differentially private synthesis via
#'   [dp_control()], with a person-level (\eqn{\epsilon}, \eqn{\delta})
#'   guarantee. See `vignette("differential-privacy")` for scope and accounting.
#'
#' Synthetic data is not anonymisation. Track A output must never be described
#' as differentially private.
#'
#' @keywords internal
#' @importFrom stats ave predict setNames
#' @importFrom utils modifyList
"_PACKAGE"
