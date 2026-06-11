# PFIM 7.1

## PFIM 7.1 (2026-06-10)

### Maintenance

- Removed parallel FIM evaluation (`future` / `furrr`, `configureParallel()`, `parallelStatus()`)
- Dropped unused OpenMP flags from `src/Makevars` (no `#pragma omp` in C++ kernels)
- PSO batch fitness: `.pfimMetaheuristicFitnessBatch()` for `pso_optimize_Rcpp` (`eval_fitness_batch`)
- Package `README.md` and `_pkgdown.yml` for local documentation builds
- Progressive `lintr` defaults (assignment, line length, return)
- C++ kernel test: PSO batch callback (Fedorov buffers in test-fedrov-cpp-buffers.R)

## PFIM 7.1 (2026-06-06)

### Release

- CRAN release 7.1 (`DESCRIPTION`, `NEWS.md`, `README.md` aligned)
- Final release maintained by Romain Leroux (original developer and maintainer, 2020–2026).
- Maintenance transferred to the PFIM team (Inserm).

### CRAN polish

- Removed non-standard `objdump` logic from `src/Makevars.win`
- Simplex C++: convergence message only when `showProcess = TRUE` (cleaner test output)
- Removed redundant `importFrom(S7, prop)` from `NAMESPACE` (`import(S7)` suffices)

### Architecture

- S7 object system for `Evaluation`, `Optimization`, `Model`, `Fim`, and related classes

### Performance

- ODE gradient performance caches: administration reuse, FD scheme cache, deSolve time grid (`PFIM.perf.*`, default on)
- FIM design cache during constraint enumeration (`options(PFIM.fim.cache = TRUE)`)
- Batch model rebuild during constraint grid evaluation
- C++ optimizer hot paths: pre-converted Armadillo structures (Multiplicative, Simplex, PSO, PGBO)

### Optimization

- Fixed Fedorov-Wynn R/C++ buffer sizing: `fisher` and population vectors now match FIM dimension (`ndimFim`), not only the number of candidate protocols

### Documentation

- New vignettes: `Example03_GCSF`, `Example04_Fakinumab` (PopED-aligned tutorials)
- Updated vignette: `FIM-types-and-covariateTest`
- Package URL points to GitHub (CRAN incoming check)
- `README.md` and `NEWS.md` at package root

### CRAN

- Removed invalid `pfim.biostat.fr` URL
- `RcppArmadillo` moved to `LinkingTo` only
- Fixed non-ASCII characters, `parallelly` declaration, and `.GlobalEnv` assignments in parallel code
- `codecov.yml` excluded from the source tarball (CI config only)

### Testing

- 215+ passing `testthat` tests including gradient performance regression, parallel (Linux/macOS CI), FIM cache, and optimizer smoke tests
- End-to-end tests for all five optimizers (Multiplicative, Fedorov-Wynn, PSO, PGBO, Simplex)
- Fedorov-Wynn C++ integration tests (buffer sizing, no subscript warnings)
- Parallel optimizer tests run on all platforms including Windows (`future::multisession`)
- Parallel workers prefer `PFIM.devPath` when set; relayed package-version warnings from workers are muffled
- GitHub Actions at repository root: multi-OS `R CMD check` and integration tests

## History

PFIM was originally developed and maintained by **Romain Leroux** (2020–2026)
as part of Inserm research. The package is now maintained by the PFIM team.
