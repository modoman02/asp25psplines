#' Fit a penalised B-spline smoother
#'
#' @param x Numeric covariate vector or design input.
#' @param knots Numeric vector of knot locations.
#' @return A list with basis/design matrix and fitted coefficients.
#' @examples
#' @examples
#' \donttest{
#' x <- seq(0, 1, length.out = 50)
#' k <- list(
#'   knots = seq(0, 1, length.out = 10),
#'   tp    = FALSE,
#'   l     = 2
#' )
#' fs <- fit_spline(x = x, knots = k)
#' }

#' @export

fit_spline <- function(x, knots) {
  tp <- knots$tp
  l  <- knots$l
  kv <- knots$knots
  n  <- length(x)
  if (tp) {
    m <- length(kv)
    d <- (l + 1) + m
    Z <- matrix(NA_real_, nrow = n, ncol = d)
    # polynomial part
    for (j in 1:(l+1)) Z[, j] <- x^(j-1)
    # truncated parts
    for (j in seq_len(m)) Z[, l + 1 + j] <- base_fun(x, knots, m = m, j = j)
  } else {
    d <- length(kv) - l - 1
    Z <- matrix(NA_real_, nrow = n, ncol = d)
    for (j in seq_len(d)) Z[, j] <- base_fun(x, knots, m = NA, j = j)
  }
  Z
}
