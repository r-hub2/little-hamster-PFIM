# PFIM 7.1

[![R-CMD-check](https://github.com/packagePFIM/PFIM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/packagePFIM/PFIM/actions/workflows/R-CMD-check.yaml)
[![Integration tests](https://github.com/packagePFIM/PFIM/actions/workflows/integration-tests.yaml/badge.svg)](https://github.com/packagePFIM/PFIM/actions/workflows/integration-tests.yaml)

**Population Fisher Information Matrix** — evaluate and optimize clinical trial designs for nonlinear mixed-effects (NLME) models using the Fisher information matrix.

## Features

- **Design evaluation**: population, individual, and Bayesian FIM with covariates and IOV
- **Design optimization**: Simplex, Multiplicative, Fedorov-Wynn, PSO, and PGBO (C++ kernels)
- **PSO batch fitness**: one R callback per swarm iteration
- **Model library**: PK, PD, and PK/PD models (analytic and ODE)
- **Reports**: HTML evaluation and optimization reports via `Report()`

## Installation

```r
devtools::install_github("packagePFIM/PFIM",
  subdir = "PFIM_7_1_beta_version/PFIM_S7_FINAL_52")
```

Requires **R >= 4.4.0** (S7) and a C++ toolchain (`Rcpp`, `RcppArmadillo`).

## Quick start

```r
library(PFIM)

vignette("Example01")             # discrete optimization
vignette("Example02")             # continuous optimization (PSO, PGBO, Simplex)
vignette("Example03_GCSF")        # G-CSF PK/PD
vignette("Example04_Fakinumab")   # 2-cpt + MM clearance
vignette("LibraryOfModels")
```

## Development

```r
Rscript tools/rebuild-rcpp.R   # after editing src/
devtools::test()
devtools::check()
Rscript tools/coverage.R       # >= 80% gate
Rscript tools/lint-package.R
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for module layout and release checklist.

## Documentation site

```r
install.packages("pkgdown")
pkgdown::build_site()
```

Configuration: `_pkgdown.yml`. The site lists exported functions and vignettes; pre-built vignette HTML is also shipped under `inst/doc/`.

## License

GPL-3. Maintained by the PFIM team (Inserm). See `NEWS.md` for version history.
