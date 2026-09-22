# PFIM 7.1

**Population Fisher Information Matrix** — evaluate and optimize clinical trial designs for nonlinear mixed-effects (NLME) models using the Fisher information matrix.

## Features

- **Design evaluation**: population, individual, and Bayesian FIM with covariates and IOV
- **Design optimization**: Simplex, Multiplicative, Fedorov-Wynn, PSO, and PGBO (C++ kernels)
- **Model library**: PK, PD, and PK/PD models (analytic and ODE)
- **Reports**: HTML evaluation and optimization reports via `Report()`

## Installation

```r
install.packages("PFIM")
```

Requires **R >= 4.5.0**, S7, and a C++ toolchain for source builds.

## Quick start

```r
library(PFIM)

vignette("Example01")       # discrete optimization
vignette("Example02")       # continuous optimization (PSO, PGBO, Simplex)
vignette("Example03")       # population FIM with covariates
vignette("LibraryOfModels") # PK/PD model library
```

Large discrete constraint grids: `pfim_set_option(constraints.maxTasks = 500)`.
Long sessions: `pfim_set_option(fim.cache.maxEntries = 500)` or `pfim_reset_session()`.

## License

GPL-3. See `NEWS.md` for version history.
