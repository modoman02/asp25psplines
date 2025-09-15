#' Base B-spline construction utility
#'
#' Low-level helper to construct B-spline basis functions.
#'
#' @param x Numeric vector of covariates.
#' @param knots List with elements `knots` (numeric knot vector), `tp` (logical thin-plate flag), and `l` (integer degree/order).
#' @param m Integer penalty/order parameter (difference order).
#' @param j Integer basis index.
#'
#' @return Numeric vector with basis evaluations for index `j`.
#' @examples
#' \donttest{
#' x <- seq(0, 1, length.out = 10)
#' K <- base_fun(
#'   x,
#'   knots = list(knots = seq(0, 1, length.out = 5), tp = FALSE, l = 2),
#'   m = 2, j = 1
#' )
#' }
#' @export
base_fun <- function(x, knots, m, j) {
  l  <- knots$l
  tp <- knots$tp
  kv <- knots$knots
  if (tp) {
    kappa_j <- kv[j]
    return(pmax(0, x - kappa_j)^l)
  } else {
    if (l == 0) {
      if (j == length(kv) - 1) {
        return(ifelse(x >= kv[j] & x <= kv[j+1], 1, 0))
      } else {
        return(ifelse(x >= kv[j] & x < kv[j+1], 1, 0))
      }
    }
    denom1 <- kv[j + l] - kv[j]
    denom2 <- kv[j + l + 1] - kv[j + 1]
    term1 <- if (denom1 == 0) 0 else ((x - kv[j]) / denom1) *
      base_fun(x = x, knots = list(knots = kv, tp = tp, l = l - 1), m = m, j = j)
    term2 <- if (denom2 == 0) 0 else ((kv[j + l + 1] - x) / denom2) *
      base_fun(x = x, knots = list(knots = kv, tp = tp, l = l - 1), m = m, j = j + 1)
    term1 + term2
  }
}
