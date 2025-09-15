#' Compact console summary for a fitted asp25psplines model
#'
#' @param result Fitted object returned by `run_all_funcs()`.
#' @return A list with key summary elements; prints a concise summary.
#' @examples
#' \donttest{
#' # res <- run_all_funcs(x, y)
#' # summary_all(res)
#' }
#' @export



summary_all <- function(result) {
  cat("=== penalised spline fit ===\n")
  cat("Iterations:", result$iterations, "\n")
  cat("lambda_mu:", format(result$lambda_mu, digits = 6),
      " | lambda_sigma:", format(result$lambda_sigma, digits = 6), "\n")
  mu <- result$mu_hat; sg <- result$sigma_hat
  if (all(is.finite(mu))) cat("mu:    mean =", round(mean(mu),4), "range = [", round(min(mu),4), ",", round(max(mu),4), "]\n")
  if (all(is.finite(sg))) cat("sigma: mean =", round(mean(sg),4), "range = [", round(min(sg),4), ",", round(max(sg),4), "]\n")
}
