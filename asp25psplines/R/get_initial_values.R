#' Initial values for location-scale fitting
#'
#' @param X Design matrix for the mean.
#' @param Z Design matrix for the log-scale.
#' @param y Numeric response vector.
#' @return A list with initial mean and scale (or their coefficients).
#' @examples
#' \donttest{
#' # init <- get_initial_values(X, Z, y)
#' }
#' @export

get_initial_values <- function(X, Z, y) {
  # OLS for mu with small ridge
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  beta_ols <- .solve_spd(XtX, Xty, ridge = 1e-8)
  mu_hat   <- as.numeric(X %*% beta_ols)
  res      <- y - mu_hat
  s_hat    <- log(pmax(abs(res), .Machine$double.eps)) + 0.635
  ZtZ <- crossprod(Z); Zts <- crossprod(Z, s_hat)
  gamma_ols <- .solve_spd(ZtZ, Zts, ridge = 1e-8)
  sigma_hat <- as.numeric(exp(Z %*% gamma_ols))

  # one weighted step
  w <- 1 / pmax(sigma_hat^2, 1e-12)
  XtWX <- crossprod(X * sqrt(w))
  XtWy <- crossprod(X * sqrt(w), y * sqrt(w))
  beta_w <- .solve_spd(XtWX, XtWy, ridge = 1e-8)
  mu_hat <- as.numeric(X %*% beta_w)
  list(mu_hat = mu_hat, sigma_hat = sigma_hat)
}
