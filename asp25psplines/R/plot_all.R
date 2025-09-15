#' Quick diagnostic plot for fitted mean and (optionally) scale
#'
#' @param result Fitted object returned by `run_all_funcs()`.
#' @param x Optional numeric x-values for plotting on a specific grid.
#' @param y Optional observed responses for residual diagnostics.
#' @param show_sigma Logical; if TRUE, draw a band derived from the estimated scale.
#' @param alpha Numeric alpha level for bands (0–1). Default: 0.1.
#' @return Invisibly returns `NULL`. Called for its side effects (plot).
#' @examples
#' \donttest{
#' # res <- run_all_funcs(x, y)
#' # plot_all(res, show_sigma = TRUE)
#' }
#' @export



plot_all <- function(result, x, y, alpha = 0.05, show_sigma = TRUE) {
  if (missing(x) || missing(y)) { x <- result$x; y <- result$y }
  ord <- order(x); x <- x[ord]; y <- y[ord]
  mu <- finite_guard(result$mu_hat)[ord]
  sg <- finite_guard(result$sigma_hat)[ord]

  plot(x, y, pch = 16, cex = 0.6, xlab = "x", ylab = "y",
       main = "Penalised spline fit")

  if (all(is.finite(mu))) lines(x, mu, lwd = 2) else message("mu contains non-real numbers. No line printed")
  if (show_sigma && all(is.finite(mu)) && all(is.finite(sg))) {
    z <- 1.96
    up <- mu + z * sg; lo <- mu - z * sg
    polygon(c(x, rev(x)), c(up, rev(lo)), border = NA, col = rgb(0, 0, 0, 0.1))
    lines(x, up, lty = 2); lines(x, lo, lty = 2)
  }

  legend("topleft", bty = "n",
         legend = c(
           paste0("lambda_mu=", format(result$lambda_mu, digits = 4)),
           paste0("lambda_sigma=", format(result$lambda_sigma, digits = 4)),
           paste0("edf_mu=", format(suppressWarnings(result$edf_mu), digits = 4)),
           paste0("edf_sigma=", format(suppressWarnings(result$edf_sigma), digits = 4)),
           paste0("iter=", result$iterations)
         ))
}
