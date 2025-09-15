#!/usr/bin/env Rscript
# ================================================================
# load_and_run_big_sim_fast.R
# 1) Load all .R functions from a folder (default: ./R) into workspace
# 2) Run BIG SIM (FAST) using those functions (no predict(), package optional)
# ------------------------------------------------
# Usage:
#   Rscript load_and_run_big_sim_fast.R                # assumes ./R exists or auto-detect
#   Rscript load_and_run_big_sim_fast.R /path/to/R     # explicit R folder
# ================================================================

options(stringsAsFactors = FALSE, warn = 1)

message("=== Step 1: Resolve R/ source folder ===")
args <- commandArgs(trailingOnly = TRUE)

# Try to resolve the R source folder robustly
find_r_dir <- function(start_dir = getwd()) {
  # 1) If ./R exists and has run_all_funcs.R, use it
  cand <- file.path(start_dir, "R")
  if (dir.exists(cand) && file.exists(file.path(cand, "run_all_funcs.R"))) return(normalizePath(cand))
  # 2) Search recursively for a folder named R that contains run_all_funcs.R
  dirs <- list.dirs(start_dir, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    if (basename(d) == "R" && file.exists(file.path(d, "run_all_funcs.R"))) {
      return(normalizePath(d))
    }
  }
  # 3) Fallback: if user provided a zip-extracted path via args, we'll try args[1]
  return(NA_character_)
}

R_DIR <- if (length(args) >= 1 && dir.exists(args[1])) normalizePath(args[1]) else find_r_dir()
if (is.na(R_DIR)) {
  stop("Could not find an 'R' folder with run_all_funcs.R. ",
       "Run as: Rscript load_and_run_big_sim_fast.R /path/to/extracted/R")
}
message("R_DIR = ", R_DIR)

message("=== Step 2: Ensure dependencies ===")
need <- function(pkgs) for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) {
  install.packages(p, repos="https://cloud.r-project.org")
}
need(c("ggplot2","splines"))

# Optional: load your package if available (not required)
if (requireNamespace("asp25psplines", quietly = TRUE)) {
  suppressPackageStartupMessages(library(asp25psplines))
  message("Loaded installed package 'asp25psplines'.")
} else {
  message("Package 'asp25psplines' not installed; using sourced functions.")
}

message("=== Step 3: Source all functions from R_DIR into workspace ===")

# Priority order to avoid dependency hiccups; remaining files source alphabetically
priority <- c(
  "zzz_utils_base.R",
  "helpers.R",
  "create_knots.R",
  "get_pen_mat.R",
  "fit_spline.R",
  "select_lambda.R",
  "get_initial_values.R",
  "update_parameters.R",
  "calc_mu.R",
  "calc_sigma.R",
  "calc_deviance.R",
  "base_fun.R",
  "run_all_funcs.R",
  "plot_all.R",
  "summary_all.R",
  "run_demo.R",
  "asppsplines-package.R"
)

all_files <- list.files(R_DIR, pattern = "\\.[rR]$", full.names = TRUE)
# Drop macOS resource forks if present
all_files <- all_files[!grepl("__MACOSX", all_files, fixed = TRUE)]
# Order: priority first if present, then remaining alphabetically
order_idx <- integer(0)
for (p in priority) {
  hit <- which(basename(all_files) == p)
  if (length(hit)) order_idx <- c(order_idx, hit[1])
}
rest <- setdiff(seq_along(all_files), order_idx)
ordered_files <- c(all_files[order_idx], all_files[rest][order(basename(all_files[rest]))])

loaded <- character(0)
for (f in ordered_files) {
  tryCatch({
    sys.source(f, envir = .GlobalEnv)
    loaded <- c(loaded, basename(f))
    message("sourced: ", basename(f))
  }, error = function(e){
    stop("Failed to source ", f, " -> ", conditionMessage(e))
  })
}

message("Loaded ", length(loaded), " files into workspace.")
if (!exists("run_all_funcs", envir = .GlobalEnv) && !("asp25psplines" %in% loadedNamespaces())) {
  stop("Function 'run_all_funcs' not found in workspace and package not loaded. Check your R folder.")
}

# Helper to resolve functions either from workspace or package namespace
resolve_fun <- function(name) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) return(get(name, envir = .GlobalEnv))
  if ("asp25psplines" %in% loadedNamespaces()) {
    env <- asNamespace("asp25psplines")
    if (exists(name, envir = env, inherits = FALSE)) return(get(name, envir = env))
  }
  stop("Function '", name, "' not found in workspace or namespace.")
}

run_all_funcs <- resolve_fun("run_all_funcs")
plot_all_safe <- function(fit) {
  fn <- NULL
  if (exists("plot_all", envir = .GlobalEnv, inherits = FALSE)) fn <- get("plot_all", envir=.GlobalEnv)
  else if ("asp25psplines" %in% loadedNamespaces() && exists("plot_all", envir=asNamespace("asp25psplines"), inherits=FALSE)) fn <- get("plot_all", envir=asNamespace("asp25psplines"))
  if (!is.null(fn)) try(fn(fit), silent=TRUE)
}

message("=== Step 4: BIG SIM configuration ===")



# ---------- Config ----------
CFG <- list(
  # Grid
  n_grid        = c(200, 500, 1000),
  tp_grid       = c("smooth","piecewise","sine"), # Wahrheits-Szenarien (Label)
  tp_flag_grid  = c(FALSE, TRUE),                 # Paket-Boolean-Schalter (50/50)
  rep_per_cell  = 80,                             # anpassen für Laufzeit
  seed          = 123,
  
  # Methoden – nur erlaubte!
  crit_mu_grid    = c("BIC","AIC","GCV","ML"),
  crit_sigma_grid = c("GCV","ML"),
  
  # Spline-Setup fürs Paket
  m_mu      = 20,   # m (Basisanzahl) – zu deinem Paket passend
  l_order   = 3,    # l (Spline-Ordnung) IMMER 3 (kubisch) – WICHTIG
  equi      = TRUE,
  buffer    = 0.05,
  r_pen     = 2,
  
  # Lambda-Grid (Basis) – dicht genug
  lambda_grid_base = exp(seq(log(1e-4), log(1e3), length.out=60)),
  
  # Refinement: dichtere Grids in Runde 2 (ohne Kenntnis “gewählter λ”)
  refine_enable = c(FALSE, TRUE),
  lambda_grid_refined = exp(seq(log(5e-5), log(5e3), length.out=120)),
  
  # Intervalle
  coverage_target = 0.95,
  leverage_inflation_cap = 1.5,
  leverage_eps = 1e-6,
  
  # Ordner
  results_dir = "results",
  fig_dir     = "figures"
)

# ---------- RUN folders ----------
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
CFG$results_dir <- file.path(CFG$results_dir, RUN_ID)
CFG$fig_dir     <- file.path(CFG$fig_dir, RUN_ID)
dir.create(CFG$results_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(CFG$fig_dir,     recursive=TRUE, showWarnings=FALSE)

# ---------- Truths & data gen ----------
truth_mu <- function(x, tp="smooth") {
  switch(tp,
         "smooth"    = 2*(x-0.5)^2 - 0.1,
         "piecewise" = ifelse(x < 0.4, -0.2 + 0.5*x,
                              ifelse(x < 0.7, 0.1, 0.1 + 0.8*(x-0.7))),
         "sine"      = 0.3*sin(2*pi*x),
         2*(x-0.5)^2 - 0.1
  )
}
truth_sigma <- function(x, tp="smooth") {
  base <- 0.3 + 0.2*(x-0.5)^2
  switch(tp,
         "smooth"    = base,
         "piecewise" = base + 0.2*(x>0.6),
         "sine"      = base + 0.1*sin(2*pi*x),
         base
  )
}
gen_data <- function(n, tp) {
  x <- sort(runif(n))
  mu <- truth_mu(x, tp)
  sg <- pmax(1e-6, truth_sigma(x, tp))
  y <- rnorm(n, mu, sg)
  list(x=x, y=y, mu=mu, sigma=sg)
}

# ---------- Metrics & intervals ----------
rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
neg_loglik <- function(y, mu, sigma) { sigma <- pmax(1e-10, sigma); -mean(dnorm(y, mean=mu, sd=sigma, log=TRUE)) }
coverage_pi <- function(y, mu, sigma, level=0.95) {
  z <- qnorm((1+level)/2); mean(y >= mu - z*sigma & y <= mu + z*sigma, na.rm=TRUE)
}
ci_length <- function(mu, sigma, level=0.95) {
  z <- qnorm((1+level)/2); mean(2*z*sigma, na.rm=TRUE)
}
inflate_sigma <- function(sigma_hat, df_mu, n) {
  fac <- sqrt(n / max(1, n - df_mu + CFG$leverage_eps))
  fac <- pmin(CFG$leverage_inflation_cap, pmax(1.0, fac))
  sigma_hat * fac
}
calibrate_sigma_factor <- function(y, mu, sigma, level=0.95) {
  z <- qnorm((1+level)/2)
  grid <- seq(0.5, 3.0, length.out=80)
  vals <- sapply(grid, function(c){ l<-mu - z*sigma*c; u<-mu + z*sigma*c; mean(y>=l & y<=u) - level })
  grid[which.min(abs(vals))]
}

# ---------- Helpers ----------
normalize_method <- function(s) {
  s <- toupper(as.character(s))
  ok <- c("AIC","BIC","ML","GCV")
  if (!s %in% ok) stop("Unbekannte Methode: ", s, " (erlaubt: ", paste(ok, collapse=", "), ")")
  s
}
is_bad_vec <- function(v) { any(!is.finite(v)) || any(is.na(v)) }
get_num <- function(obj, name, default=NA_real_) if (!is.null(obj[[name]]) && is.numeric(obj[[name]]) && length(obj[[name]])==1) obj[[name]] else default

# ---------- (Very small) internal fallback – only if pkg-fit fails/NA ----------
# Keep minimal to avoid interference; we prefer 100% package usage.
bs_basis <- function(x, k, degree=3) {
  knots <- if (k > degree+1) seq(min(x), max(x), length.out=k - degree + 1L)[-c(1,(k-degree+1L))] else NULL
  splines::bs(x, degree=degree, knots=knots, intercept=TRUE)
}
diff_matrix <- function(k, m) { D <- diag(k); for (i in seq_len(m)) D <- diff(D); D }
pspline_fit <- function(x, y, k, m, lambda) {
  X <- bs_basis(x, k=k, degree=3); k_eff <- ncol(X)
  D <- diff_matrix(k_eff, m); P <- crossprod(D)
  XtX <- crossprod(X); XtY <- crossprod(X, y)
  A <- XtX + lambda*P + 1e-8*diag(k_eff)
  beta <- solve(A, XtY); yhat <- as.vector(X%*%beta)
  df <- sum(diag(solve(A, XtX)))
  list(beta=beta, X=X, df=df, yhat=yhat)
}
fallback_once <- function(x,y) {
  # very small stable fallback: one lambda for mu via GCV, then log-σ spline via GCV
  k_mu <- 20; k_sg <- 7; m <- 2
  lam_grid <- exp(seq(log(1e-4), log(1e3), length.out=40))
  best <- list(score=Inf, fm=NULL, lam=NA_real_)
  for (lam in lam_grid) {
    fm <- pspline_fit(x,y,k=k_mu,m=m,lambda=lam)
    gcv <- sum((y-fm$yhat)^2) / (max(1, length(y)-fm$df))^2
    if (gcv < best$score) best <- list(score=gcv,fm=fm,lam=lam)
  }
  mu_hat <- best$fm$yhat
  r2 <- pmax(1e-12,(y-mu_hat)^2); z <- log(r2)
  bests <- list(score=Inf, fs=NULL, lam=NA_real_)
  for (lam in lam_grid) {
    fs <- pspline_fit(x,z,k=k_sg,m=m,lambda=lam)
    gcv <- sum((z-fs$yhat)^2) / (max(1, length(z)-fs$df))^2
    if (gcv < bests$score) bests <- list(score=gcv,fs=fs,lam=lam)
  }
  list(mu_hat=mu_hat, sigma_hat=pmax(1e-8, sqrt(exp(bests$fs$yhat))),
       edf_mu=best$fm$df, edf_sigma=bests$fs$df,
       lambda_mu=best$lam, lambda_sigma=bests$lam, backend="fallback")
}

# ---------- Core fit (package first; optional refined λ-grid; NA-guard) ----------
fit_pkg_round <- function(x,y, tp_flag, method_mu, method_sigma, lambda_grid) {
  run_all_funcs(
    x=x, y=y,
    m=CFG$m_mu, l=CFG$l_order,
    equi=CFG$equi, tp=tp_flag, buffer=CFG$buffer, r_pen=CFG$r_pen,
    lambda_grid_mu=lambda_grid, lambda_grid_sigma=lambda_grid,
    method_mu=method_mu, method_sigma=method_sigma,
    max_iterations=50, tolerance=1e-5
  )
}

fit_with_pkg_then_maybe_refine <- function(x,y,tp_flag, method_mu,method_sigma, refine=FALSE) {
  mu_meth <- normalize_method(method_mu)
  sg_meth <- normalize_method(method_sigma)
  
  fit1 <- try(fit_pkg_round(x,y,tp_flag, mu_meth, sg_meth, CFG$lambda_grid_base), silent=TRUE)
  if (inherits(fit1,"try-error")) return(fit1)
  
  # NA guard
  mu_hat <- fit1$mu_hat
  sigma_hat <- if (!is.null(fit1$sigma_hat)) fit1$sigma_hat else rep(sd(y - mu_hat), length(mu_hat))
  if (is_bad_vec(mu_hat) || is_bad_vec(sigma_hat)) {
    return(structure(simpleError("NA/Inf in mu_hat/sigma_hat after pkg round 1"), class="try-error"))
  }
  
  if (!refine) { fit1$backend <- "pkg"; return(fit1) }
  
  # Refinement: zweite Runde mit dichterem/breiterem Grid
  fit2 <- try(fit_pkg_round(x,y,tp_flag, mu_meth, sg_meth, CFG$lambda_grid_refined), silent=TRUE)
  if (inherits(fit2,"try-error")) { fit1$backend <- "pkg"; return(fit1) }
  
  mu2 <- fit2$mu_hat
  sg2 <- if (!is.null(fit2$sigma_hat)) fit2$sigma_hat else rep(sd(y - mu2), length(mu2))
  if (is_bad_vec(mu2) || is_bad_vec(sg2)) { fit1$backend <- "pkg"; return(fit1) }
  
  fit2$backend <- "pkg"
  fit2
}

# ---------- One replicate ----------
run_one <- function(n, tp_lab, tp_flag, method_mu, method_sigma, refine, seed) {
  set.seed(seed)
  dat <- gen_data(n, tp_lab)
  
  t0 <- proc.time()[[3]]
  pkg_fit <- try(fit_with_pkg_then_maybe_refine(dat$x, dat$y, tp_flag, method_mu, method_sigma, refine=refine), silent=TRUE)
  
  used_backend <- "pkg"
  used_stage   <- if (refine) "pkg_refined" else "pkg_base"
  err_msg <- ""
  
  if (inherits(pkg_fit,"try-error")) {
    err_msg <- conditionMessage(attr(pkg_fit, "condition"))
    fb <- fallback_once(dat$x, dat$y)
    pkg_fit <- fb
    used_backend <- "fallback"
    used_stage   <- "fallback"
  } else {
    # zusätzliche NA-Guard direkt nach pkg
    mu_test <- pkg_fit$mu_hat
    sg_test <- if (!is.null(pkg_fit$sigma_hat)) pkg_fit$sigma_hat else rep(sd(dat$y - mu_test), length(mu_test))
    if (is_bad_vec(mu_test) || is_bad_vec(sg_test)) {
      err_msg <- "NA/Inf after pkg fit; using fallback."
      fb <- fallback_once(dat$x, dat$y)
      pkg_fit <- fb
      used_backend <- "fallback"
      used_stage   <- "fallback"
    }
  }
  
  elapsed <- proc.time()[[3]] - t0
  
  # Extract hats (robust)
  mu_hat <- pkg_fit$mu_hat
  sigma_hat <- if (!is.null(pkg_fit$sigma_hat)) pkg_fit$sigma_hat else rep(sd(dat$y - mu_hat), length(mu_hat))
  
  # EDF / lambdas if present
  edf_mu    <- get_num(pkg_fit,"edf_mu",NA_real_)
  edf_sigma <- get_num(pkg_fit,"edf_sigma",NA_real_)
  lam_mu    <- get_num(pkg_fit,"lambda_mu",NA_real_)
  lam_sigma <- get_num(pkg_fit,"lambda_sigma",NA_real_)
  
  # Intervals: raw, leverage-inflated, calibrated
  z <- qnorm((1+CFG$coverage_target)/2)
  
  cov_raw <- coverage_pi(dat$y, mu_hat, sigma_hat, CFG$coverage_target)
  ci_raw  <- ci_length(mu_hat, sigma_hat, CFG$coverage_target)
  
  sigma_infl <- if (is.finite(edf_mu)) inflate_sigma(sigma_hat, edf_mu, n) else sigma_hat
  cov_infl <- coverage_pi(dat$y, mu_hat, sigma_infl, CFG$coverage_target)
  ci_infl  <- ci_length(mu_hat, sigma_infl, CFG$coverage_target)
  
  cfac <- calibrate_sigma_factor(dat$y, mu_hat, sigma_hat, CFG$coverage_target)
  cov_cal <- coverage_pi(dat$y, mu_hat, sigma_hat*cfac, CFG$coverage_target)
  ci_cal  <- ci_length(mu_hat, sigma_hat*cfac, CFG$coverage_target)
  
  # Print progress line
  cat(sprintf("n=%-4d tp_lab=%-9s tp=%-5s mu=%-3s sg=%-3s refine=%-5s | backend=%-8s stage=%-11s | cov_raw=%.3f cov_cal=%.3f | %s\n",
              n, tp_lab, as.character(tp_flag), toupper(method_mu), toupper(method_sigma),
              as.character(refine), used_backend, used_stage, cov_raw, cov_cal, err_msg))
  
  data.frame(
    n=n, tp=tp_lab, tp_flag=tp_flag, seed=seed,
    criterion_mu=toupper(method_mu), criterion_sigma=toupper(method_sigma),
    refine=refine, backend=used_backend, stage=used_stage,
    time_sec=elapsed,
    rmse_mu = rmse(mu_hat, dat$mu),
    nll     = neg_loglik(dat$y, mu_hat, sigma_hat),
    
    coverage_raw = cov_raw,  ci_len_raw  = ci_raw,
    coverage_infl= cov_infl, ci_len_infl = ci_infl,
    coverage_cal = cov_cal,  ci_len_cal  = ci_cal,
    
    edf_mu=edf_mu, edf_sigma=edf_sigma,
    lambda_mu=lam_mu, lambda_sigma=lam_sigma,
    stringsAsFactors=FALSE
  )
}

# ---------- Grid + Chunks + Resume ----------
set.seed(CFG$seed)
GRID <- expand.grid(
  n = CFG$n_grid,
  tp = CFG$tp_grid,
  tp_flag = CFG$tp_flag_grid,
  criterion_mu = CFG$crit_mu_grid,
  criterion_sigma = CFG$crit_sigma_grid,
  refine = CFG$refine_enable,
  rep = seq_len(CFG$rep_per_cell),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

done_file <- file.path(CFG$results_dir, "progress_done.txt")
if (!file.exists(done_file)) file.create(done_file)
done <- tryCatch(scan(done_file, what=integer(), quiet=TRUE), error=function(e) integer())

chunk_size <- 500L
rows <- nrow(GRID); all_chunks <- ceiling(rows/chunk_size)
cat("RUN_ID: ", RUN_ID, "\n", sep="")
cat("results_dir: ", CFG$results_dir, "\n", sep="")
cat("fig_dir    : ", CFG$fig_dir, "\n", sep="")
cat("Total rows: ", rows, " in ", all_chunks, " chunks\n", sep="")

pkg_err_log <- file.path(CFG$results_dir, "pkg_errors.log")

for (chunk in seq_len(all_chunks)) {
  idx_start <- (chunk-1L)*chunk_size + 1L
  idx_end   <- min(rows, chunk*chunk_size)
  sel <- idx_start:idx_end
  
  res_list <- list()
  for (i in sel) {
    if (i %in% done) next
    g <- GRID[i, ]
    seed_i <- CFG$seed + i
    rr <- tryCatch({
      row <- run_one(n=g$n, tp_lab=g$tp, tp_flag=g$tp_flag,
                     method_mu=g$criterion_mu, method_sigma=g$criterion_sigma,
                     refine=g$refine, seed=seed_i)
      row$id <- i; row$rep <- g$rep
      row
    }, error=function(e){
      msg <- conditionMessage(e)
      cat(sprintf("[FATAL] row %d error: %s\n", i, msg))
      cat(sprintf("row=%d n=%d tp=%s tp_flag=%s mu=%s sg=%s refine=%s | %s\n",
                  i, g$n, g$tp, as.character(g$tp_flag), g$criterion_mu, g$criterion_sigma,
                  as.character(g$refine), msg),
          file=pkg_err_log, append=TRUE)
      NULL
    })
    if (!is.null(rr)) {
      res_list[[length(res_list)+1L]] <- rr
      cat(sprintf("Row %d/%d written (backend=%s)\n", i, rows, rr$backend[1]))
      cat(i, file=done_file, append=TRUE); cat("\n", file=done_file, append=TRUE)
    }
  }
  if (length(res_list)) {
    res <- do.call(rbind, res_list)
    out_path <- file.path(CFG$results_dir, sprintf("results_chunk_%03d.csv", chunk))
    utils::write.csv(res, out_path, row.names=FALSE)
    cat("Wrote ", out_path, "\n", sep="")
  }
}

# Merge
csvs <- list.files(CFG$results_dir, pattern="^results_chunk_\\d+\\.csv$", full.names=TRUE)
if (length(csvs)) {
  all <- do.call(rbind, lapply(csvs, utils::read.csv, stringsAsFactors=FALSE))
  out_merged <- file.path(CFG$results_dir, "sim_results_merged.csv")
  utils::write.csv(all, out_merged, row.names=FALSE)
  cat("Merged -> ", out_merged, "\n", sep="")
} else {
  cat("No chunk CSVs found; nothing to merge.\n")
}

# ---------- Plots ----------
suppressPackageStartupMessages(library(ggplot2))
in_file <- file.path(CFG$results_dir, "sim_results_merged.csv")
if (!file.exists(in_file)) stop("Keine Ergebnisse gefunden: ", in_file)
DT <- utils::read.csv(in_file, stringsAsFactors=FALSE)
DT$refine <- as.logical(DT$refine)
DT$criterion_mu    <- factor(DT$criterion_mu, levels=CFG$crit_mu_grid)
DT$criterion_sigma <- factor(DT$criterion_sigma, levels=CFG$crit_sigma_grid)
DT$tp              <- factor(DT$tp, levels=CFG$tp_grid)
DT$tp_flag         <- as.logical(DT$tp_flag)

save_plot <- function(p, name, w=10, h=6) {
  fp <- file.path(CFG$fig_dir, paste0(name, ".png"))
  ggsave(fp, p, width=w, height=h, dpi=150)
  message("Saved: ", fp)
}

mline <- function() stat_summary(fun=mean, geom="line")
mpoint<- function() stat_summary(fun=mean, geom="point")

# 1) RMSE(mu) vs n by criterion_mu (facet: tp)
p1 <- ggplot(DT, aes(x=factor(n), y=rmse_mu, color=criterion_mu, group=criterion_mu)) +
  mline() + mpoint() + facet_wrap(~tp, scales="free_y") +
  labs(x="n", y="RMSE(mu)", title="RMSE(mu) vs n by criterion_mu") +
  theme_bw()
save_plot(p1, "rmse_mu_vs_n_by_criterion_mu")

# 2) Coverage raw/infl/cal vs n by criterion_sigma (facet tp)
DL <- rbind(
  transform(DT[,c("n","tp","criterion_sigma","coverage_raw")],  kind="raw",  val=coverage_raw)[,c("n","tp","criterion_sigma","kind","val")],
  transform(DT[,c("n","tp","criterion_sigma","coverage_infl")], kind="infl", val=coverage_infl)[,c("n","tp","criterion_sigma","kind","val")],
  transform(DT[,c("n","tp","criterion_sigma","coverage_cal")],  kind="cal",  val=coverage_cal)[,c("n","tp","criterion_sigma","kind","val")]
)
p2 <- ggplot(DL, aes(x=factor(n), y=val, color=criterion_sigma, group=criterion_sigma)) +
  mline() + mpoint() + geom_hline(yintercept=0.95, linetype=2) +
  facet_grid(kind ~ tp) +
  labs(x="n", y="Coverage (95%)", title="Coverage: raw vs leverage-inflated vs calibrated") +
  theme_bw()
save_plot(p2, "coverage_all_vs_n_by_criterion_sigma")

# 3) CI lengths (raw/infl/cal)
CL <- rbind(
  transform(DT[,c("n","tp","criterion_sigma","ci_len_raw")],  kind="raw",  val=ci_len_raw)[,c("n","tp","criterion_sigma","kind","val")],
  transform(DT[,c("n","tp","criterion_sigma","ci_len_infl")], kind="infl", val=ci_len_infl)[,c("n","tp","criterion_sigma","kind","val")],
  transform(DT[,c("n","tp","criterion_sigma","ci_len_cal")],  kind="cal",  val=ci_len_cal)[,c("n","tp","criterion_sigma","kind","val")]
)
p3 <- ggplot(CL, aes(x=factor(n), y=val, color=criterion_sigma, group=criterion_sigma)) +
  mline() + mpoint() + facet_grid(kind ~ tp) +
  labs(x="n", y="Average CI length", title="CI length: raw vs infl vs cal") +
  theme_bw()
save_plot(p3, "ci_length_all_vs_n_by_criterion_sigma")

# 4) EDF, Zeit
p4 <- ggplot(DT, aes(x=factor(n), y=edf_mu, color=criterion_mu, group=criterion_mu)) +
  mline()+mpoint()+facet_wrap(~tp)+labs(x="n", y="edf_mu", title="EDF(mu) vs n")+theme_bw()
save_plot(p4, "edf_mu_vs_n")
p5 <- ggplot(DT, aes(x=factor(n), y=edf_sigma, color=criterion_sigma, group=criterion_sigma)) +
  mline()+mpoint()+facet_wrap(~tp)+labs(x="n", y="edf_sigma", title="EDF(sigma) vs n")+theme_bw()
save_plot(p5, "edf_sigma_vs_n")
p6 <- ggplot(DT, aes(x=factor(n), y=time_sec, color=interaction(criterion_mu,criterion_sigma), group=interaction(criterion_mu,criterion_sigma))) +
  mline()+mpoint()+facet_wrap(~tp)+labs(x="n", y="Seconds", color="mu:sigma", title="Elapsed time")+theme_bw()
save_plot(p6, "time_vs_n_by_criteria")

# 5) tp_flag-Effekt (nur Kontrolle)
p7 <- ggplot(DT, aes(x=factor(n), y=rmse_mu, color=factor(tp_flag), group=tp_flag)) +
  mline()+mpoint()+facet_grid(tp ~ interaction(criterion_mu,criterion_sigma))+
  labs(x="n", y="RMSE(mu)", color="tp_flag", title="RMSE vs n by tp_flag (package boolean)") + theme_bw()
save_plot(p7, "rmse_mu_by_tpflag")

cat("Alle Ergebnisse gespeichert unter:\n  - CSVs:   ", CFG$results_dir,
    "\n  - Plots:  ", CFG$fig_dir, "\n", sep="")