
# Wrapper, that makes alternating updates for mu and sigma with lambda selection
#' Update parameters for alternating location-scale fitting
#'
#' Internal update engine used by high-level wrappers.
#'
#' @param X TODO: describe parameter 'X'.
#' @param Z TODO: describe parameter 'Z'.
#' @param y TODO: describe parameter 'y'.
#' @param K\_mu TODO: describe parameter 'K\_mu'.
#' @param K\_sigma TODO: describe parameter 'K\_sigma'.
#' @param max\_iterations TODO: describe parameter 'max\_iterations'.
#' @param lambda\_grid\_mu TODO: describe parameter 'lambda\_grid\_mu'.
#' @param lambda\_grid\_sigma TODO: describe parameter 'lambda\_grid\_sigma'.
#' @param method\_mu TODO: describe parameter 'method\_mu'.
#' @param method\_sigma TODO: describe parameter 'method\_sigma'.
#' @param tolerance TODO: describe parameter 'tolerance'.
#' @return TODO: describe return value.
#' @examples
#' \donttest{
#' # minimal example (pseudo)
#' # res <- update_parameters()
#' }
#' @keywords internal
update_parameters <- function(X, Z, y,
                              K_mu, K_sigma,
                              max_iterations = 50,
                              lambda_grid_mu = exp(seq(log(1e-4), log(1e4), length.out = 40)),
                              lambda_grid_sigma = exp(seq(log(1e-4), log(1e4), length.out = 40)),
                              method_mu = "BIC",
                              method_sigma = "BIC",
                              tolerance = 1e-5) {
  
  n <- NROW(X)
  
  init <- get_initial_values(X, Z, y)
  mu_hat    <- as.numeric(init$mu_hat)
  sigma_hat <- pmax(as.numeric(init$sigma_hat), 1e-6)
  
  lam_mu_seq    <- rep(NA_real_, max_iterations)
  lam_sigma_seq <- rep(NA_real_, max_iterations)
  edf_mu_seq    <- rep(NA_real_, max_iterations)
  edf_si_seq    <- rep(NA_real_, max_iterations)
  
  beta_last  <- rep(0, ncol(X))
  gamma_last <- rep(0, ncol(Z))
  
  dev_prev <- Inf
  iters <- 0
  
  for (it in seq_len(max_iterations)) {
    # update mu
    sel_mu <- select_lambda(y = y, X = X, K = K_mu,
                                   sigma_hat = sigma_hat,
                                   lambda_grid = lambda_grid_mu,
                                   method = method_mu)
    lam_mu <- sel_mu$lambda_best
    w <- 1 / pmax(sigma_hat^2, 1e-12)
    Xw <- X * sqrt(w)
    XtWX <- crossprod(Xw) + lam_mu * K_mu
    XtWy <- crossprod(Xw, y * sqrt(w))
    ch <- tryCatch(chol(XtWX), error=function(e) NULL)
    if (!is.null(ch)) {
      beta  <- backsolve(ch, forwardsolve(t(ch), XtWy))
      edf_mu <- edf_exact_from_chol(ch, crossprod(Xw))
    } else {
      beta  <- .solve_spd(XtWX, XtWy, ridge = 1e-8)
      edf_mu <- suppressWarnings(edf_hutchinson(X, Wsqrt = sqrt(w), K = K_mu, lambda = lam_mu, n_probe = 15))
    }
    mu_hat <- as.numeric(X %*% beta)
    
    lam_mu_seq[it] <- lam_mu
    edf_mu_seq[it] <- edf_mu
    
    # update sigma 
    r2 <- pmax((y - mu_hat)^2, .Machine$double.eps)
    y_sigma <- 0.5 * log(r2)
    sel_si <- select_lambda(y = y_sigma, X = Z, K = K_sigma,
                                   sigma_hat = rep(1, n),
                                   lambda_grid = lambda_grid_sigma,
                                   method = method_sigma)
    lam_si <- sel_si$lambda_best
    ZtZ <- crossprod(Z) + lam_si * K_sigma
    Zts <- crossprod(Z, y_sigma)
    ch2 <- tryCatch(chol(ZtZ), error=function(e) NULL)
    if (!is.null(ch2)) {
      gamma <- backsolve(ch2, forwardsolve(t(ch2), Zts))
      edf_si <- edf_exact_from_chol(ch2, crossprod(Z))
    } else {
      gamma <- .solve_spd(ZtZ, Zts, ridge = 1e-8)
      edf_si <- suppressWarnings(edf_hutchinson(Z, Wsqrt = rep(1, n), K = K_sigma, lambda = lam_si, n_probe = 15))
    }
    sigma_hat <- link_sigma(as.numeric(Z %*% gamma))
    
    lam_sigma_seq[it] <- lam_si
    edf_si_seq[it] <- edf_si
    
    dev_now <- calc_deviance(y, mu_hat, sigma_hat)
    if (is.finite(dev_prev) && is.finite(dev_now)) {
      rel <- abs(dev_prev - dev_now) / (abs(dev_prev) + 1e-8)
      if (rel < tolerance) { iters <- it; break }
    }
    dev_prev <- dev_now
    beta_last <- beta; gamma_last <- gamma
    iters <- it
  }
  
  list(mu_hat = mu_hat, sigma_hat = sigma_hat,
       lambda_mu = lam_mu_seq[iters], lambda_sigma = lam_sigma_seq[iters],
       edf_mu = edf_mu_seq[iters], edf_sigma = edf_si_seq[iters],
       beta = beta_last, gamma = gamma_last,
       iterations = iters,
       y = y)
}
