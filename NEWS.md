# PFIM 7.1

### Analytic / ODE / optimization

- Passive (non-administered) analytic equations that depend on `t` use
  `t_<that outcome>` at the sampling grid, not `t_<administered>`.
  Disease-progression / placebo baseline catalogue models evaluate again.
- Library ODE Michaelis-Menten first-order and infusion remap `dose_RespPK`
  / `Tinf_RespPK` onto the administered compartment (`C1` or
  `outputs$RespPK`). A leftover mismatch errors with a `PFIM:` prefix.
- Discrete optimizers (Fedorov-Wynn / Multiplicative) derive a singleton dose
  grid from the arm's current `Administration` when
  `administrationsConstraints` is omitted. Continuous optimizers reject
  discrete-style `SamplingTimeConstraints` without `samplingsWindows`.
- `SamplingTimeConstraints(numberOfSamplingsOptimisable=)` is an alias of
  the historical `numberOfsamplingsOptimisable`.
- A declared covariate effect of 0 stays as a `\beta` column in the population
  FIM (SE under the null). A warning is issued; omit the parameter to drop it.
- `SimplexAlgorithm` builds extra vertices inside sampling windows (and reseeds
  when `initialSamplings` violate window counts), so increasing
  `numberOfTimesByWindows` no longer freezes at the start. An all-infeasible
  amoeba is an error, not a false `converged = TRUE` at iteration 0.
- When `samplingsWindows` is set, `SamplingTimeConstraints()` requires
  `initialSamplings` to match the window counts and `minSampling`.
  `checkValiditySamplingConstraint()` also checks that the arm's current
  sampling times occupy those windows. Occupancy splits the sorted times by
  `numberOfTimesByWindows`, so a shared endpoint is not double-counted
  (contiguous partitions such as `[0,4]+[4,12]` are accepted).
  Windows must be declared in increasing time order.
- Continuous optimizers warn when `Arm` `samplingTimes` differ from
  `SamplingTimeConstraints$initialSamplings` (the actual search start).
- `numberOfsamplingsOptimisable` is the total protocol size including
  fixed times, not the count of free times.
- Design FIM cache keys include ODE `initialConditions`.

### Maintenance / CRAN

- Package maintainer (`cre`) for the CRAN submission: France Mentré
  (`PFIM@inserm.fr`).
- Romain Leroux (`aut`) designed and coded PFIM 7.1 in full (R/S7 package,
  C++ kernels, tests, vignettes, and documentation). He stops development
  after this release; subsequent versions will not be developed by him.

### Individual / Bayesian FIM (subject-level)

- Across arms and covariate strata, subject FIMs are aggregated by averaging
  covariances then inverting
  (\eqn{\bar C=\sum w_s M_s^{-1}}, \eqn{M_{\mathrm{eff}}=\bar C^{-1}}), not by
  mixing Fisher matrices. Singular strata contribute infinite variance on the
  null-space parameters (eigen / Inf), so SE is Inf when any positive weight
  falls on a non-identifiable protocol. Null-space detection uses the
  correlation-scale eigen-decomposition (\eqn{R=D^{-1/2}MD^{-1/2}}, relative
  threshold \eqn{10^{-10}}) so unit-dependent roundoff is not treated as
  information.
- Covariate `\beta` is omitted from individual and Bayesian FIMs (MAP conditions
  on `\beta`); use population FIM for covariate tests.
- Bayesian eta set is `\omega > 0` regardless of `fixedMu` / `fixedOmega`
  (population fix flags do not remove subject random effects). FD gradients
  are computed for those `\mu` even when `fixedMu`, and also when `\gamma > 0`
  with `\omega = 0` (IOV residual still needs `\partial f/\partial\mu`).
- Bayesian SE/RSE: on the `\eta` scale; for LogNormal, displayed
  `\mathrm{SE}_\theta=\mu\cdot\mathrm{SE}_\eta` and `\mathrm{RSE}=100\cdot\mathrm{SE}_\eta`.
- Individual and Bayesian IOV: marginal FO residual
  `V_k=R_k+F_k\,\mathrm{diag}(\gamma^2)\,F_k^\top` for every `\gamma>0`
  (including parameters without IIV); data FIM keeps eta / mu columns only.
- Discrete optimizers (Multiplicative / Fedorov-Wynn) pick the single best
  subject protocol for individual / Bayesian designs (covariance-mixture optimum
  is a vertex).
- Fedorov-Wynn default optimality gap `\delta` is `1e-4` (matches Wynn weight
  precision). Single-observation finite-difference gradients stay matrices.

### Public API

PFIM 7.1 uses the S7 object system. Pipeline helpers are **not** on the
public API. Use `run()`, `Report()`,
`get*` / `plot*` instead. Removed names stay documented as `\keyword{internal}`
(`PFIM:::name` for package developers).

- `Dcriterion()` → `getDcriterion()`
- `plotSEFIM()` / `plotRSEFIM()` → `plotSE()` / `plotRSE()`
- `plotWeightsMultiplicativeAlgorithm()` → `plotWeights()`
- `plotFrequenciesFedorovWynnAlgorithm()` → `plotFrequencies()`
- `plotEvaluationResults()` / `plotEvaluationSI()` → `plotEvaluation()` /
  `plotSensitivityIndices()`
- `generateReportEvaluation()` / `generateReportOptimization()` /
  `tablesForReport()` → `Report()`
- `evaluateDesign()` / `evaluateArm()` / `evaluateModel()` / `evaluateFim()` /
  `setEvaluationFim()` → `run(Evaluation)`
- Individual / Bayesian FIM aggregation is a covariance mixture
  (\eqn{\bar C=\sum w_s M_s^{-1}}), not a sum of Fisher matrices; not comparable
  to PFIM 6 / PopED on that case. Non-identifiable protocols with positive
  weight give \code{Inf} SE.
- `optimizeDesign()` / `generateFimsFromConstraints()` → `run(Optimization)`
- `show()` as a PFIM export → `methods::show()` (do not call `PFIM::show`)
- `FedorovWynnAlgorithm_Rcpp()` / `MultiplicativeAlgorithm_Rcpp()` → internal
- `fun.amoeba` / `fisherSimplex` / `computeVMat` / `plot` → removed
- Other pipeline names (`defineModelAdministration`, `defineModelType`,
  `evaluateModelGradient`, `generateDosesCombination`, `getArmConstraints`,
  `setOptimalArms`, `updateSamplingTimes`, \ldots) are likewise internal.

### User-visible changes

- `ModelError` S7 hierarchy: `varianceForm` (`combined1` / `combined2`) drives
  residual variance; shared `.pfimSigmaIsEstimable()` rule for FIM columns and
  dV/dsigma; subclass validators lock Constant / Proportional / Combined*
  constraints; unused `equation` / `derivatives` slots removed (legacy args
  warn and are ignored).
- `Combined2` residual error: variance
  `sigmaInter^2 + (sigmaSlope * f)^2` (PopED-style additive + proportional),
  with residual diagonals computed in Rcpp (`residualErrorDerivatives_Rcpp`).
- Example `inst/examples/gcsf-poped-evaluation.R` / vignette **Example04**:
  G-CSF / filgrastim PK-PD population FIM matching the PopED workshop Design 1 RSE
  (linear FD stencil, with optional quadratic-FD columns).
- S7 object model (`Evaluation`, `Optimization`, `Model`, `Fim`, covariates, library).
  Requires **R >= 4.5.0** (Deriv 4.3.5+).
- Population, individual, and Bayesian FIM; categorical covariates and IOV
  (`CategoricalCovariate`, `CategoricalCovariateWithIOV`, `gamma` on parameters).
- `covariateTest()` / `tost()` for covariate effect tests.
- Five optimizers: Multiplicative, Fedorov-Wynn, PSO, PGBO, Simplex (C++ kernels).
- Joint discrete optimization over multiple arms in one design.
- Discrete allocation uses Hamilton (largest-remainder) rounding so subject
  counts sum exactly to the study size; console reports match that rule.
- Optimizer outputs share nested `algorithmOutput` (status / converged /
  iterations). Discrete algos also store `optimalWeights` and
  `numberOfSubjects` (Fedorov-Wynn reads legacy `numberOfIndividuals`;
  top-level `frequencies` is still written as a plot alias). Continuous
  algos (PSO / PGBO / Simplex) have no mixture weights.
- Continuous optimizers: when `tolerance <= 0`, stall/relative stop is off and
  `converged` is `NA` (not a false failure). Multiplicative returns certified
  weights (optimality checked before the update). Fedorov-Wynn distinguishes
  stationary `success` from usable `incomplete` (cycle budget) with a warning.
- Multiplicative outputs keep `mixtureDcriterion` and `realisedDcriterion`, with
  accessors `getMixtureDcriterion()` / `getRealisedDcriterion()`. Joint
  multi-outcome cells keep the best-weighted protocol only.
- HTML report figures: 16 pt grey-theme plots, 10×5 in
  chunks (not 14 in / tiny SE-RSE fonts), nested ggplots printed one-by-one.
- `Report()`, vignettes (`Example01`, `Example02`, `Example03`, `Example04`, `LibraryOfModels`).
- Session options: `pfim_get_option()`, `pfim_set_option()`, `pfim_reset_session()`.
- `numberOfOccasions` on `Evaluation` and `Optimization` (`NA` infers; explicit
  values are checked against `inferNumberOfOccasions()`).

### Performance

- Evaluation / SI HTML figures use a dense time grid
  (\code{seq(0, tmax, 0.05)}, capped at 400 points) with design samples as red
  markers. FIM evaluation itself stays at the design sampling times.
- Evaluation / continuous optimization hot paths: reuse FD nominal on the flat
  `evaluateArm` path (one fewer model solve per arm); simple population FIM via
  `computePopFimCombo_Rcpp`; leaner `.applyFlatToArms` / FD evaluation loop on
  the fitness path.
- Residual variance assembly keeps `Matrix::Diagonal` / sparse embedded sigma
  blocks until FIM kernels densify (avoids dense zero fills on large designs).
- Continuous optimizers (Simplex / PSO / PGBO) reuse one evaluation context per
  design search instead of cloning eval+design on every fitness call.
- Fedorov-Wynn C++ reuses packed/dense FIM workspaces across Wynn exchanges.
- FIM cache (`fim.cache`), batch rebuild during constraint grids (`eval.batch`).
- ODE caches: `perf.adminCache`, `perf.fdCache` (FD `XcolsInv` keyed on
  `nFree`, not mu), `perf.odeTimesCache` (sim-time grid in C++).
  Constraint-grid FIM evals reuse the prepared model (cloned per arm inside
  `evaluateDesign()`).
- `mfvar` / `mfvar_mixed` keep `T_k = V^{-1} dV_k V^{-1}` and stream `dV`
  from the R list (one conversion per column: n reads, not n²/2 copies).
- PSO: one R fitness callback per swarm iteration (`eval_fitness_batch`).
- Covariate/IOV arms: model, gradient, and variance on one combination × occasion pass.

### Bug fixes

- PK/PD catalogue: two-compartment infusion closed forms inline hybrid rates
  (no free `alpha`/`beta`); first-order SS micro-constants use `ka`, `k`,
  `k12`, `k21`, `V` (not `Cl`/`Q`/`V1`/`V2`); 1-cpt Michaelis–Menten infusion
  is the vignette ODE (`MichaelisMenten1InfusionSingleDose_VmKmV`; former
  2-cpt name kept as alias); sigmoid turnover keys list `gamma` and use `^`.
  2-cpt MM infusion C1 includes peripheral transfer; quadratic PD uses `Aquad`;
  1-cpt MM bolus vignette name `..._VmKmV` is an alias of `..._VmKm`.
  Two-compartment bolus closed forms wrap hybrid rates (`B/beta` is not
  `B/0.5 * (...)`). Exponential baseline `increase` rises to `S0` and
  `decrease` decays from `S0`. 1-cpt infusion SS Cl uses `/Cl` (not
  `/((Cl/V)*V)`). 2-cpt Michaelis–Menten ODEs use central volume `V1`
  (catalogue suffix `V1V2`), not `V`.
  Unknown `modelFromLibrary` names error instead of evaluating `NULL` equations.
  Constructors reject `NA` / non-finite `tau`, `outcome`, `gamma`, `omega`,
  residual sigmas / `cError`, and arm/design size instead of crashing later
  in `if (NA)`. Cache scope, option names, analytic time names, and report
  HTML treat `NA` like empty (`nzchar(NA)` is TRUE).
  `Arm(initialCondition=)` is an official alias of `initialConditions`.
- `prop<-` on `Optimization` now writes nested project fields (`name`,
  `fimType`, `designs`, …) and returns the updated object (S7 value semantics).
- ODE bolus \code{y} follows \code{Deriv_*} declaration order, not the dosed
  compartment; first-order absorption with the depot declared second now runs.
  Event-based bolus still zeroes the dosed compartment IC at \code{t = 0}
  (mass comes from the event table). Library 2-compartment ODE names keep the
  unmapped catalogue state (no \code{Deriv_NA}).
- `Administration(outcome=)` accepts a compartment, a `dose_*` target, or an
  `outputs` alias that maps to a state (library `"RespPK"` included).
- Unknown `initialConditions` names are rejected. Categorical covariate
  effects and IOV sequences must name declared categories; proportions must
  be non-negative.
- ODE initial-condition expressions evaluate in a sealed environment and reject
  unknown symbols (no session leakage).
- `Report()` restores knitr chunk options and does not close the user's
  graphics devices. `saveCovariateTest()` returns the output path invisibly.
- `src/init.c` uses `(void)` prototypes for zero-argument routines.

- Population FIM with a fixed typical value and unfixed IIV
  (`fixedMu = TRUE`, `\omega > 0`): the finite-difference gradient of that
  `\mu` is now computed, so `V = F \Omega F^\top + R` keeps the `\omega^2`
  term. The previous zero gradient dropped that IIV from `V` (about 72%
  error versus the reference with `\mu_V` fixed and `\omega_V^2 = 0.1`).
- Population FIM: `\omega^2` stays estimable when `\mu` is fixed
  (`fixedMu` and `fixedOmega` are independent). Previously, fixing `\mu`
  also dropped the `\omega^2` row.
- Individual / Bayesian: a rank-deficient subject FIM (e.g. one sample for
  two `\mu`) now reports Inf SE on the unidentifiable directions. The old
  `\varepsilon\cdot p\cdot\lambda_{\max}` cutoff treated a `10^{-13}`
  roundoff eigenvalue as information and returned huge finite SE.

- Residual-error FIM columns are labelled $\sigma$ (SD scale), not $\sigma^2$.
  Values and SE/RSE were already $\partial V/\partial\sigma$; only the display
  prefix was wrong (console, plots, HTML).
- `Dcriterion()` on an `Optimization` after `run()` no longer errors: the nested
  project `fim` is filled from the optimal evaluation, `fisherMatrix` is
  validated as a square matrix, and `.fimDcriterionFromMatrix()` returns 0 on
  empty / non-matrix input.
- `facet_wrap(space = "free_x")` is omitted on ggplot2 3.5.x (the argument is
  ggplot2 >= 4.0.0). The declared floor stays `ggplot2 (>= 3.5.0)`.
- Armadillo: `-DARMA_WARN_LEVEL=1`, `solve_opts::no_approx`, and explicit
  symmetrization before Cholesky (stderr warnings no longer bypass R conditions).
  Non-finite matrices are rejected before `chol()` (NA fails Armadillo's
  symmetry check because `NA != NA`).
- PK/PD covariate tests accept ASCII / Unicode / `<U+03BC>` parameter labels
  under `LC_CTYPE=C`.
- HTML model-error table no longer emits a stray `NA` column after
  `varianceForm` was added to `getModelErrorData()`.
- Condition numbers of rank-deficient FIM blocks are `Inf` (near-zero
  eigenvalues are no longer dropped to produce a finite pseudo-κ).
- Finite-difference `|mu| < 1e-4` warns once per distinct tiny-μ signature
  that the absolute step floor can dominate truncation error; scale parameters
  to O(1). A second badly scaled design in the same session still warns.
- S7 devel (#409): `new_object()` parents are instances (`S7_object()` /
  `ModelError(...)`), not class objects, so constructors work on S7 0.2.2.9000.
- ggplot2 version for `facet_wrap(space=)` is cached in `.pfimSession`.
- G-CSF / PopED Design 1 RSE is checked on every `R CMD check`
  (`tests/testthat/test-gcsf-poped.R`, tolerance `1e-3`).
- Fast Rd examples (constructors, residual variance, covariate links) run
  during `R CMD check`; heavy FIM / optimiser examples stay `\donttest`.

- Joint discrete plots: Multiplicative / Fedorov-Wynn mixture bars are
  protocol-level when `optimalArms` is expanded; `plotWeights()` /
  `plotFrequencies()` error clearly on the wrong optimizer.
- Joint Multiplicative weights: shrink to the winning protocol when the
  constraint cell is joint (multi-outcome). Empty or zero-arm cells are a no-op.
- Population joint Multiplicative: study `N` is split across arms of the
  winning protocol via `.allocProportionalSubjects` (aligned with Fedorov-Wynn).
- ODE finite-difference gradients: no longer replace ε^(1/3) by
  `odeTol^(1/3)` (that inflated FD steps and biased the FIM); steps are
  only raised when they fall below `10 * atol/rtol`.
- Population covariate/IOV FIM: C++ kernels (`computePopFimCombo_Rcpp`, `computeMFVar_Rcpp`) with R cross-tests; `chol_inv` / `safe_solve` pure Armadillo (match R `chol2inv`/`chol` and `solve`+jitter); `src/init.c` auto-generated from Rcpp exports; IOV block accumulation fix.
- Fedorov-Wynn C++: stray text removed from `FedorovWynnAlgorithm.cpp` (Windows build).
- Fedorov-Wynn: `fisher` buffers sized to FIM dimension (`ndimFim`).
- Individual/Bayesian covariate×occasion FIM: named occasion helpers and
  one BLAS `tcrossprod` for FO IOV inflation of V (same algebra as rank-1
  updates). Population mixture is `reduce` of combination blocks.

### Packaging / CRAN

- Public API freeze: pipeline generics (`evaluateDesign`, `optimizeDesign`,
  `generateFimsFromConstraints`, `setEvaluationFim`, \ldots) are internal.
  User verbs remain `run()`, `Report()`, `get*` / `plot*` (including
  `plotSensitivityIndices()`), and constructors (`FedorovWynnAlgorithm`,
  `MultiplicativeAlgorithm`, PSO / PGBO / Simplex). Tests and Rd examples of
  unexported internals call `PFIM:::`. `pfim_set_option()` rejects unknown names.
- User-facing errors go through `.pfimStop()` (`PFIM: \ldots`); no internal
  `.fn:` prefixes leak from subject-level FIM mixing.
- Frozen gold numbers in `inst/fixtures/gold-fim.R` (G-CSF Design 1 RSE and
  cas1–cas6 D-criteria). Do not edit the table to silence a failing test.

- Analytic Model* evaluation: shared time×outcome grid, dose tables, and
  infusion superposition (`R/model-analytic-eval.R`). Infusion steady-state
  no longer uses `assign()` into the evaluator environment.

- S7 constructors pass the parent to `new_object()` positionally (CRAN S7
  names that argument `.parent`; S7 devel uses `_parent`).
- Optional GitHub Actions job `s7-devel` installs `pak::pak("RConsortium/S7")`
  and is allowed to fail.
- S7 validators on remaining structural classes (`SamplingTimes`,
  `AdministrationConstraints`, `Distribution`, `Covariate`, `Model`,
  `PFIMProject`).
- Split `Fim.R` (labels / discrete arms / HTML render) and `pfim-fim-cache.R`
  (continuous-optimizer helpers).
- `stats` added to `Imports` (aligned with `NAMESPACE`).
- `stats::setNames()` everywhere (R CMD check NOTE).
- `tost` documented as alias of `covariateTest()`.
- Version floors: `ggplot2 (>= 3.5.0)` (`facet_wrap(space=)` gated for 4.0.0),
  `S7 (>= 0.1.1)`, `testthat (>= 3.2.3)` in `DESCRIPTION`.
- `@examples` on `run()`, `Design`, `Arm` (see `inst/examples/evaluation-minimal.R`).
- FD grid includes estimable `\mu`, any `\omega > 0`, and any `\gamma > 0`;
  population FIM keeps `\omega^2` when `\mu` is fixed (`fixedOmega` is the
  flag that drops IIV).
- `mfvar_mixed`: rank-1 omega/gamma columns use `(v' V^{-1} w)^2` instead of
  dense `p x p` outer products; wired in `computePopFimCombo_Rcpp` and the
  simple population path. R reference `.computePopFimCombo_R` aligned.
- `constraints.maxTasks` documented (`?pfim_set_option`); subsampling on large
  constraint grids is deterministic (equispaced per dose), not RNG-based.
- Default ODE tolerances when `odeSolverParameters = list()`.
- `optimizerParameters` validated at `Optimization()` construction; dead config
  properties removed from `*Algorithm` S7 classes (settings live on `Optimization` only).
- Optional `fim.cache.maxEntries` LRU cap on FIM caches.
- Vignette HTML / figures default to `tempdir()` (`PFIM_VIGNETTE_RESULTS` to
  keep a local `results/` copy). `results/` is not shipped in the tarball.
- Population FIM variance block uses `adjustGradient()` per distribution
  (Normal: no μ scaling; LogNormal: chain rule on η). Residual variance
  `(σ_inter + σ_slope·f^c)^2` with `cError`; `Constant` / `Proportional`
  constructors reject mixed σ parameters.

### Maintenance

- FIM cache scopes: `pfim#N` on `cacheScope` (`pfim-fim-cache.R`).
- C++ pop FIM / linalg kernels; `src/init.c` routine registration.
- Maintainer (`cre`): France Mentré. Author (`aut`): Romain Leroux, who
  coded PFIM 7.1 in full and stops development after this release.

## History

Romain Leroux coded PFIM 7.1 in full and stops development for subsequent
versions. PFIM 7.1 is the S7-based release documented in the package
vignettes and `?PFIM-package`.
