#' Create knot locations for P-spline basis construction
#'
#' @param x Numeric covariate vector.
#' @param a Numeric lower boundary of the covariate range (left endpoint of \eqn{[a,b]}).
#' @param b Numeric upper boundary of the covariate range (right endpoint of \eqn{[a,b]}).
#' @param m Integer penalty/order parameter for differences.
#' @param l Integer number of interior knots or basis size control.
#' @param equi Logical; if TRUE, equally spaced knots, otherwise quantile-based.
#' @param tp Logical; if TRUE, use thin-plate style option if implemented.
#' @param buffer Numeric; extra range added beyond \eqn{[a,b]} for boundary knots.
#' @return A numeric vector of knot locations.
#' @examples
#' \donttest{
#' create_knots(x = runif(50), a = 0, b = 1, m = 2, l = 10, equi = TRUE, tp = FALSE, buffer = 0.05)
#' }
#' @export
create_knots <- function(x, a, b, m, l, equi = TRUE, tp = FALSE, buffer = 0) {
  ## ... dein Code ...
}


create_knots <- function(x, a = min(x), b = max(x), m = 20, l = 3,
                         equi = TRUE, tp = FALSE, buffer = 0.05) {

  # extend boundaries
  range <- max(x) - min(x)
  a <- min(x) - buffer * range
  b <- max(x) + buffer * range

  if (equi) {
    inner_knots <- seq(a, b, length.out = m)
  } else {
    probs <- seq_len(m) / (m + 1)
    inner_knots <- as.numeric(quantile(x, probs = probs, names = FALSE))
  }

  if (!tp) {
    # add l boundary knots on each side, if tp=F is selected
    left  <- rep(a, l)
    right <- rep(b, l)
    knots <- c(left, inner_knots, right)
  } else {
    # if tp=T, no boundary knots are needed
    knots <- inner_knots
  }
  list(knots = knots, tp = tp, l = l)
}
