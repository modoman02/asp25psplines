#' High-level fit for Gaussian location–scale P-splines
#'
#' Fits the conditional mean and log-scale by alternating penalised weighted
#' least squares. Smoothing parameters are selected per component; an optional
#' refinement step can calibrate interval coverage.
#'
#' @param x Numeric covariate vector.
#' @param y Numeric response vector of the same length as `x`.
#' @param m Integer difference order (penalty).
#' @param l Integer basis size/degree controller (as used in your implementation).
#' @param equi Logical; if TRUE, uses equally spaced knots, otherwise quantile-based.
#' @param tp Logical; thin-plate style option if implemented.
#' @param buffer Numeric boundary buffer beyond \eqn{[a,b]} for boundary knots.
#' @param r_pen Numeric or matrix; additional penalty specification if applicable.
#' @param lambda_grid_mu Numeric vector of candidate smoothing values for the mean.
#' @param lambda_grid_sigma Numeric vector of candidate smoothing values for the scale.
#' @param method_mu Character; selection criterion for `lambda_mu` (e.g., "GCV", "BIC").
#' @param tolerance Numeric convergence tolerance for the alternating updates.
#' @param method_sigma Character; selection criterion for `lambda_sigma`.
#' @param max_iterations Integer; maximum number of alternating iterations.
#' @return An object of class `asp25psplines_fit` with elements such as
#'   `mu_hat`, `sigma_hat`, `coef_mu`, `coef_sigma`, `lambda_mu`, `lambda_sigma`,
#'   and `diagnostics`.
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- sort(runif(80)); y <- sin(2*pi*x) + rnorm(80, 0.25)
#' res <- run_all_funcs(x, y, m = 2, l = 10, equi = TRUE, tp = FALSE,
#' buffer = 0.05, lambda_grid_mu = NULL, lambda_grid_sigma = NULL,
#' method_mu = "GCV", method_sigma = "GCV", max_iterations = 20)
#' }
#' @export



run_all_funcs <- function(x, y,
                          m = 20, l = 3,
                          equi = TRUE, tp = FALSE,
                          buffer = 0.05,
                          r_pen = 2,
                          lambda_grid_mu = exp(seq(log(1e-4), log(1e4), length.out = 40)),
                          lambda_grid_sigma = exp(seq(log(1e-4), log(1e4), length.out = 40)),
                          method_mu = "BIC",
                          method_sigma = "BIC",
                          max_iterations = 50,
                          tolerance = 1e-5) {

  stopifnot(length(x) == length(y))

  pre <- prep_xy(x, y)
  x2 <- pre$x; y2 <- pre$y

  knots <- create_knots(x2, m = m, l = l, equi = equi, tp = tp, buffer = buffer)
  X <- fit_spline(x2, knots)
  Z <- X

  K <- get_pen_mat(knots, r = r_pen)
  K_mu <- make_spd(K); K_sigma <- make_spd(K)

  res <- update_parameters(X, Z, y2,
                           K_mu = K_mu, K_sigma = K_sigma,
                           max_iterations = max_iterations,
                           lambda_grid_mu = lambda_grid_mu,
                           lambda_grid_sigma = lambda_grid_sigma,
                           method_mu = method_mu,
                           method_sigma = method_sigma,
                           tolerance = tolerance)

  res$x <- x2; res$y <- y2
  res$X <- X; res$Z <- Z
  res$knots <- knots
  class(res) <- "asppsplines_fit"
  res
}
