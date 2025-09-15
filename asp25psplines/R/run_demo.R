#' Minimal demo for asp25psplines
#'
#' @param seed Integer random seed for reproducibility.
#' @return A fitted object (class \code{asp25psplines\_fit}) from a tiny example.
#' @examples
#' \donttest{
#' run_demo(seed = 1)
#' }
#' @export

run_demo <- function(seed = 1) {
  set.seed(seed)
  n <- 1000
  x <- runif(n, -1, 1)
  y <- sin(3 * x) + stats::rnorm(n, 0.25)

  d <- prep_xy_local(x, y)   # lokale Helper-Funktion (s.u.)
  x <- d$x; y <- d$y
  # in run_demo():
  set.seed(seed)
  x <- sort(runif(100))
  mu <- sin(2*pi*x)
  y  <- mu + rnorm(100, sd = 0.2)
  # ... fit auf (x,y) ...
  fit <- run_all_funcs(
    x, y,
    tp = TRUE, m = 20, l = 3,
    method_mu = "BIC",
    method_sigma = "BIC",
    lambda_grid_mu    = exp(seq(log(1e-3), log(1e3), length.out = 30)),
    lambda_grid_sigma = exp(seq(log(1e-3), log(1e3), length.out = 30)),
    max_iterations = 50,
    tolerance = 1e-5
  )

  # Zusammenfassung + Plot (falls vorhanden)
  if (exists("summary_all")) try(summary_all(fit), silent = TRUE)
  if (exists("plot_all"))    try(plot_all(fit, x, y), silent = TRUE)

  cat("Demo finished.\n")
  invisible(fit)
}

# Lokale, robuste Vorverarbeitung; vermeidet Abhängigkeit auf externes prep_xy()
#' prep\_xy\_local
#'
#' Internal helper function.
#'
#' @param x TODO: describe parameter 'x'.
#' @param y TODO: describe parameter 'y'.
#' @param sort\_x TODO: describe parameter 'sort\_x'.
#' @return TODO: describe return value.
#' @examples
#' \donttest{
#' # minimal example (pseudo)
#' # res <- prep_xy_local()
#' }
#' @keywords internal
prep_xy_local <- function(x, y, sort_x = TRUE) {
  stopifnot(length(x) == length(y))
  ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x <- as.numeric(x[ok]); y <- as.numeric(y[ok])
  if (sort_x) {
    o <- order(x); x <- x[o]; y <- y[o]
  }
  list(x = x, y = y)
}
