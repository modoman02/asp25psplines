
# Base-R utilities: stability, solvers, ED, preprocessing (fast variants)

#' @keywords internal
#' @noRd

.sanitize_lambda_grid <- function(grid, lo=1e-6, hi=1e6, fallback_len=25L) {
  g <- as.numeric(grid)
  g <- g[is.finite(g) & g > 0]
  if (!length(g)) g <- exp(seq(log(1e-3), log(1e3), length.out = fallback_len))
  g[g < lo] <- lo
  g[g > hi] <- hi
  unique(sort(g))
}

#' @keywords internal
#' @noRd

finite_guard <- function(x, val = NA_real_) {
  x[!is.finite(x)] <- val
  x
}

#' @keywords internal
#' @noRd

link_sigma <- function(eta) {
  pmax(exp(eta), 1e-8)
}

#' @keywords internal
#' @noRd

make_spd <- function(K, eps = 1e-8) {
  K <- 0.5*(K + t(K))
  eva <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  if (min(eva) < eps) K <- K + (eps - min(eva)) * diag(nrow(K))
  K
}

#' @keywords internal
#' @noRd

.solve_spd <- function(A, b = NULL, ridge = 1e-8, ridge_max = 1e2, qr_tol = 1e-10) {
  if (any(!is.finite(A))) stop(".solve_spd: A contains non-numbers")
  if (!is.null(b) && any(!is.finite(b))) stop(".solve_spd: b contains non-numbers")
  
  A <- 0.5 * (A + t(A))
  dmean <- mean(abs(diag(A)))
  if (!is.finite(dmean) || dmean <= 0) dmean <- 1
  base <- dmean * max(1e-12, ridge)
  
  lambda <- base
  for (i in 0:12) {
    Ai <- A + (lambda * diag(nrow(A)))
    ch <- tryCatch(chol(Ai), error = function(e) NULL)
    if (!is.null(ch)) {
      if (is.null(b)) return(ch)
      y <- backsolve(ch, forwardsolve(t(ch), b))
      return(as.numeric(y))
    }
    lambda <- lambda * 10
    if (lambda > ridge_max) break
  }
  
  if (is.null(b)) {
    q <- qr(A + base * diag(nrow(A)))
    if (q$rank < nrow(A)) warning(".solve_spd: QR rank-deficient")
    return(q)
  } else {
    q <- qr(A + base * diag(nrow(A)), tol = qr_tol, LAPACK = TRUE)
    if (q$rank == 0) stop(".solve_spd: QR rank = 0")
    return(as.numeric(qr.coef(q, b)))
  }
}

# Exact ED if Cholesky available
#' @keywords internal
#' @noRd

edf_exact_from_chol <- function(ch, XtWX0) {
  S <- backsolve(ch, forwardsolve(t(ch), XtWX0))
  sum(diag(S))
}

# Hutchinson ED so no big diag of W is neccessary
#' @keywords internal
#' @noRd

edf_hutchinson <- function(X, Wsqrt, K, lambda, n_probe = 15) {
  n <- NROW(X)
  acc <- 0
  XtWX0 <- crossprod(X * Wsqrt)
  for (i in 1:n_probe) {
    z <- sample(c(-1,1), n, replace = TRUE)
    v  <- crossprod(X * Wsqrt, z)
    A  <- XtWX0 + lambda * K
    u  <- .solve_spd(A, v, ridge = 1e-8)
    Sz <- (X %*% u)
    acc <- acc + sum(z * Sz)
  }
  acc / n_probe
}

#' @keywords internal
#' @noRd

nr_update_damped <- function(theta, score, hess, max_halves = 8) {
  step <- tryCatch(.solve_spd(hess, score), error = function(e) rep(NA_real_, length(theta)))
  if (any(!is.finite(step))) return(list(theta = theta, ok = FALSE, reason = "solve_failed"))
  alpha <- 1
  for (k in 0:max_halves) {
    cand <- theta + alpha * step
    if (all(is.finite(cand)) && all(abs(cand) < 1e6)) return(list(theta = cand, ok = TRUE, alpha = alpha))
    alpha <- alpha * 0.5
  }
  list(theta = theta, ok = FALSE, reason = "no_finite_step")
}

#' @keywords internal
#' @noRd

prep_xy <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  o <- order(x); x <- x[o]; y <- y[o]
  if (any(duplicated(x))) {
    dup <- duplicated(x)
    x[dup] <- x[dup] + rnorm(sum(dup), sd = 1e-8)
    o <- order(x); x <- x[o]; y <- y[o]
  }
  xm <- mean(x); xs <- sd(x); if (!is.finite(xs) || xs == 0) xs <- 1
  list(x = (x - xm)/xs, y = y, center = xm, scale = xs)
}
