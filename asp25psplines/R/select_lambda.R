
# robust lambda selection with fallback
#' select\_lambda
#'
#' Internal helper function.
#'
#' @param y TODO: describe parameter 'y'.
#' @param X TODO: describe parameter 'X'.
#' @param K TODO: describe parameter 'K'.
#' @param sigma\_hat TODO: describe parameter 'sigma\_hat'.
#' @param lambda\_grid TODO: describe parameter 'lambda\_grid'.
#' @param method TODO: describe parameter 'method'.
#' @param ridge TODO: describe parameter 'ridge'.
#' @return TODO: describe return value.
#' @examples
#' \donttest{
#' # minimal example (pseudo)
#' # res <- select_lambda()
#' }
#' @keywords internal
select_lambda <- function(y, X, K, sigma_hat, lambda_grid,
                                 method = c("AIC","BIC","ML","GCV"),
                                 ridge = 1e-8) {
  method <- match.arg(method)
  n <- NROW(X); p <- NCOL(X)
  lambda_grid <- .sanitize_lambda_grid(lambda_grid)
  
  w <- 1 / pmax(sigma_hat^2, 1e-12)
  Wsqrt <- sqrt(w)
  Xw <- X * Wsqrt
  yw <- y * Wsqrt
  
  XtWX0 <- crossprod(Xw)            
  XtWy  <- crossprod(Xw, yw)        
  
  crit_values <- rep(NA_real_, length(lambda_grid))
  edf_vec     <- rep(NA_real_, length(lambda_grid))
  beta_list   <- vector("list", length(lambda_grid))
  
  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]
    A <- XtWX0 + lam * K
    
    # Try Cholesky first
    ch <- tryCatch(chol(A), error = function(e) NULL)
    if (!is.null(ch)) {
      beta <- backsolve(ch, forwardsolve(t(ch), XtWy))
      mu  <- as.numeric(X %*% beta)
      rss <- sum(w * (y - mu)^2)
      
      # Exact ED via helper function
      edf <- edf_exact_from_chol(ch, XtWX0)
      
    } else {
      # Fallback, in case method before didnt work 
      beta <- tryCatch(.solve_spd(A, XtWy, ridge = ridge), error = function(e) rep(NA_real_, p))
      if (any(!is.finite(beta))) next
      mu  <- as.numeric(X %*% beta)
      rss <- sum(w * (y - mu)^2)
      
      # Hutchinson ED 
      edf <- suppressWarnings(edf_hutchinson(X, Wsqrt = Wsqrt, K = K, lambda = lam, n_probe = 15))
    }
    
    if (!is.finite(edf)) edf <- min(p, n - 1)
    
    if (method %in% c("AIC","BIC")) {
      dev <- n * log(rss / n)
      if (method == "AIC") crit <- dev + 2 * edf else crit <- dev + log(n) * edf
    } else if (method == "ML") {
      crit <- rss
    } else {
      crit <- rss / (n - edf)^2
    }
    
    crit_values[j] <- crit
    edf_vec[j]     <- edf
    beta_list[[j]] <- beta
  }
  
  # Fallback, in case no real lambda could be selected 
  if (all(!is.finite(crit_values))) {
    mid <- lambda_grid[ceiling(length(lambda_grid)/2)]
    return(list(lambda_best = mid, crit_values = crit_values, edf = edf_vec, beta = rep(0, p)))
  }
  
  idx <- which.min(crit_values)
  list(lambda_best = lambda_grid[idx],
       crit_values = setNames(crit_values, paste0("lam=", format(lambda_grid, digits = 4))),
       edf = edf_vec,
       beta = beta_list[[idx]])
}
