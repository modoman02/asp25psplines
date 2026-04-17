---
output: github_document
---

# asp25psplines
*Gaussian location–scale regression with P-splines*

**asp25psplines** provides a reproducible implementation of semiparametric **location–scale** modeling using **P-splines**.  
It supports automatic smoothing selection, prediction intervals, and diagnostic tooling for simulation and applied work.

> Status: research/teaching package (not on CRAN). API subject to change.

---

## Table of Contents

- [Features](#features)  
- [Installation](#installation)  
- [Quick start](#quick-start)  
- [Model specification](#model-specification)  
- [API](#api)  
- [Diagnostics](#diagnostics)  
- [Simulation workflow](#simulation-workflow)  
- [Reproducibility](#reproducibility)  
- [Citation](#citation)  
- [Contributing](#contributing)  
- [License](#license)  
- [Acknowledgements](#acknowledgements)

---

## Features

- **Gaussian location–scale** models with B-/P-spline smooths (default: cubic, spline order \(\ell = 3\)).
- **Smoothing / model selection**: `AIC`, `BIC`, `ML`/`REML`, `GCV` (configurable).
- **Prediction** with standard errors and (prediction) intervals.
- **Diagnostics**: residual checks, coverage summaries (raw vs. calibrated), and plotting helpers.
- **Simulation support**: scripts to run large experiments, incremental saving, chunked logging, and post-hoc plotting.

---

## Installation

### Install from a local `.tar.gz` using RStudio (GUI) 

This is the easiest route if you download the file `asp25psplines_0.1.0.tar.gz` on your computer.

1. Open **RStudio**.
2. Go to **Tools → Install Packages…**.
3. In **Install from:** choose **Package Archive File (.zip; .tar.gz)**.
4. Click **Browse…** and select your `asp25psplines_0.1.0.tar.gz` file (e.g., from `Downloads`).
5. (Optional but recommended) Tick **Install dependencies**.
6. Click **Install**.
7. After it finishes, load the package:

``` r
  library(asp25psplines)
   packageVersion("asp25psplines")
#> [1] '0.1.0'
```

---

## Quick start


``` r
library(asp25psplines)

# 1) Synthetic data
set.seed(1)
n  <- 400
x  <- sort(runif(n, 0, 1))
mu <- sin(2*pi*x)
sd <- 0.3 + 0.2*cos(2*pi*x)
y  <- rnorm(n, mean = mu, sd = sd)
dat <- data.frame(x, y)

# 2) Fit location–scale model
fit <- fit_asppsplines(
  formula   = y ~ s(x, k = 20, degree = 3),   # location
  sigma     = ~ s(x, k = 20, degree = 3),     # scale
  data      = dat,
  selection = "BIC"                            # "AIC", "GCV", "ML" also supported
)

# 3) Predict with 95% prediction intervals
newx <- data.frame(x = seq(0, 1, length.out = 200))
pred <- predict_asppsplines(
  object   = fit,
  newdata  = newx,
  interval = "prediction",
  level    = 0.95
)

# 4) Simple plot
plot(newx$x, pred$fit, type = "l", xlab = "x", ylab = "y")
lines(newx$x, pred$lwr, lty = 2)
lines(newx$x, pred$upr, lty = 2)
points(dat$x, dat$y, pch = 16, cex = 0.4)
```

---

## Model specification

We model both **mean** and **scale** as smooth functions of covariates via P-splines:

- Location: $ \mu(x) = \mathbf{B}_\mu(x)\,\boldsymbol{\beta}_\mu $
- Scale (log-scale recommended): $ \log \sigma(x) = \mathbf{B}_\sigma(x)\,\boldsymbol{\beta}_\sigma $

with B-spline bases $ \mathbf{B} $ (typically cubic, order 3) and difference-penalties on coefficients.  
Smoothing parameters are chosen by information criteria or risk proxies (e.g., AIC, BIC, ML/REML, GCV).

---

## API

- `fit_asppsplines(formula, sigma, data, selection = c("BIC","AIC","GCV","ML"), ...)`  
  Fit location–scale model with P-splines; returns an object of class `"asppsplines_fit"`.

- `predict_asppsplines(object, newdata, interval = c("none","confidence","prediction"), level = 0.95, ...)`  
  Predictions with optional intervals and standard errors.

- `coef()`, `fitted()`, `residuals()`, `vcov()`  
  Standard extractors for model inspection.

- `plot()` / `autoplot()`  
  Quick diagnostics and smooth visualization.

---

## Diagnostics

- **Residuals**: assess structure/misspecification (e.g., vs. fitted or covariates).
- **Coverage**: compare nominal vs. empirical coverage and calibration summaries.
- **Smoothing check**: review chosen $\lambda$ (smoothing parameter) under different selectors.


``` r
res <- residuals(fit)
plot(fitted(fit), res, xlab = "Fitted", ylab = "Residuals")
abline(h = 0, lty = 2)

# if available:
# cov <- coverage_summary(fit)  # returns raw and calibrated coverage
# print(cov)
```

---

## Simulation workflow

Large-scale experiments are provided under `inst/sim/`:

- `inst/sim/load_and_run_big_sim_fast.R` – load package functions and run the full simulation.
- Chunked saving, resume-on-restart, and progress logging.
- Final scripts create coverage/RMSE summaries and figures.

Usage (terminal):

```bash
Rscript inst/sim/load_and_run_big_sim_fast.R              # uses ./R by default
# or
Rscript inst/sim/load_and_run_big_sim_fast.R /path/to/R   # explicit R folder
```

Outputs (CSV/RDS/plots) are written to `inst/sim/out/` by default (configurable).

---

## Reproducibility

- Set seeds where stochastic components are used.
- Record environment info:


``` r
sessionInfo()
```

- Consider `renv` for dependency management:


``` r
# install.packages("renv")
renv::init()
```

---

## Citation

Please cite as:


``` r
citation("asp25psplines")
```

---

## Contributing

Contributions are welcome:
1. Open an issue (bug/feature).
2. Create a feature branch (`feature/<short-name>`).
3. Add tests under `tests/testthat/` and keep `R CMD check` clean.
4. Submit a pull request.

For development:


``` r
devtools::document()
devtools::test()
devtools::check()
# pkgdown site (optional):
# pkgdown::build_site()
```

---

## License

MIT License – see `LICENSE`.

---

## Acknowledgements

Built for an Advanced Statistical Programming seminar; inspired by standard P-splines methodology and location–scale modeling literature.
