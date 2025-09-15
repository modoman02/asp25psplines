
#' Create Penalty Matrix for respective spline method 
#'
#' For B-splines a discrete difference penalty of order \code{r} is used.
#' For truncated power splines the unpenalised part is the polynomial of
#' degree 0..l; the truncated parts (one per inner knot) are penalised with
#' an identity penalty.
#'
#' @param knots A list as returned by \code{create\_knots()} with entries
#'   \code{knots}, \code{tp}, and \code{l}.
#' @param r Integer order of differences for P-splines (default 2).
#'
#' @return A symmetric positive semi-definite penalty matrix.
#' @export
get_pen_mat <- function(knots, r = 2) {
  tp <- knots$tp
  l  <- knots$l
  if (tp) {
    m <- length(knots$knots)
    d <- (l + 1) + m
    K <- diag(c(rep(0, l + 1), rep(1, m)), nrow = d, ncol = d)
  } else {
    d <- length(knots$knots) - l - 1
    D <- diag(d)
    for (i in seq_len(r)) {
      D <- D[-1, , drop = FALSE] - D[-nrow(D), , drop = FALSE]
    }
    K <- t(D) %*% D
  }
  K
}
