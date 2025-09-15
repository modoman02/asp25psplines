#' Update step for the mean
#'
#' Penalised weighted least squares update for the mean curve.
#'
#' @param X Design matrix for the mean.
#' @param y Numeric response vector.
#' @param K\_mu Penalty matrix for the mean.
#' @param sigma\_hat Current estimate of the standard deviation (vector).
#' @param lambda\_grid Numeric vector or scalar of smoothing values.
#' @param method Character; selection method for smoothing (e.g., "GCV","BIC").
#' @param max\_iterations\_mu Integer; max iterations for the mean step.
#' @param tolerance Numeric convergence tolerance.
#' @return A list with updated coefficients, fitted values and selected lambda.
#' @examples
#' \donttest{
#' # pseudo:
#' # out <- calc_mu(X, y, K_mu, sigma_hat, lambda_grid=10^seq(-3,3), method="GCV",
#' #                max_iterations_mu=50, tolerance=1e-6)
#' }
#' @export

calc_mu <- function(X, y, K_mu, sigma_hat,
                    lambda_grid, method = "AIC",
                    max_iterations_mu = 50,
                    tolerance = 1e-6) {
  n <- NROW(X); d <- NCOL(X)
  sel <- select_lambda(y = y, X = X, K = K_mu,
                              sigma_hat = sigma_hat,
                              lambda_grid = lambda_grid,
                              method = method)
  lam <- sel$lambda_best
  w <- 1 / (sigma_hat^2)
  W <- diag(w, n, n)
  A <- t(X) %*% W %*% X + lam * K_mu + 1e-8 * diag(d)
  beta <- as.numeric(solve(A, t(X) %*% W %*% y))
  mu_hat <- as.numeric(X %*% beta)
  list(mu_new = mu_hat,
       beta = beta,
       lambda_hat = lam,
       edf = sel$edf[which(lambda_grid == lam)][1],
       crit_curve = sel$crit_values)
}
