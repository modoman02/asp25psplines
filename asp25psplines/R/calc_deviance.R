#' Calculate deviance for Gaussian location-scale model
#'
#' @param y Numeric response vector.
#' @param mu\_hat Numeric vector of fitted means.
#' @param sigma\_hat Numeric vector of fitted standard deviations.
#' @return Numeric deviance-like scalar (or vector) for the fit.
#' @examples
#' \donttest{
#' y <- rnorm(10); mu <- rep(0,10); sig <- rep(1,10)
#' calc_deviance(y, mu, sig)
#' }
#' @export

calc_deviance <- function(y, mu_hat, sigma_hat) {
  sum(log(sigma_hat^2) + (y - mu_hat)^2 / (sigma_hat^2))
}
