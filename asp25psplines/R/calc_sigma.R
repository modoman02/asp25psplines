#' Update step for the log-scale
#'
#' Penalised weighted least squares update for the log-scale.
#'
#' @param Z Design matrix for the log-scale.
#' @param y Numeric response vector.
#' @param K\_sigma Penalty matrix for the scale.
#' @param mu\_hat Current estimate of the mean (vector).
#' @param lambda\_grid Numeric vector or scalar of smoothing values.
#' @param method Character; selection method for smoothing (e.g., "GCV","BIC").
#' @param max\_iterations\_sigma Integer; max iterations for the scale step.
#' @param tolerance Numeric convergence tolerance.
#' @return A list with updated coefficients, fitted log-scale and selected lambda.
#' @examples
#' \donttest{
#' # pseudo:
#' # out <- calc_sigma(Z, y, K_sigma, mu_hat, 10^seq(-3,3), "GCV", 50, 1e-6)
#' }
#' @export

calc_sigma <- function(Z, y, K_sigma, mu_hat,
                       lambda_grid, method = "AIC",
                       max_iterations_sigma = 50,
                       tolerance = 1e-6) {
  n <- NROW(Z); d <- NCOL(Z)
  r2 <- (y - mu_hat)^2
  sel <- select_lambda(y = log(pmax(r2, .Machine$double.eps))/2,
                              X = Z, K = K_sigma,
                              sigma_hat = rep(1, n),
                              lambda_grid = lambda_grid,
                              method = method)
  lam <- sel$lambda_best
  W <- diag(1, n)
  A <- t(Z) %*% W %*% Z + lam * K_sigma + 1e-8 * diag(d)
  gamma_hat <- as.numeric(.solve_spd(A, t(Z) %*% (log(pmax(r2, .Machine$double.eps))/2), ridge = 1e-8))
  sigma_hat <- as.numeric(exp(Z %*% gamma_hat))
  list(sigma_new = sigma_hat,
       gamma = gamma_hat,
       lambda_hat = lam,
       edf = sel$edf[which.min(abs(lambda_grid - lam))],
       crit_curve = sel$crit_values)
}
