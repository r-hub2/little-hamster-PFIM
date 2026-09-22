# Example03: evaluation, optimization, plots, and reports.
utils = system.file("vignette-scripts", "pfim-vignette-utils.R", package = "PFIM")
if (!nzchar(utils)) stop("pfim-vignette-utils.R not found.", call. = FALSE)
source(utils, local = TRUE)
library(PFIM)

paths = pfimVignetteSetupPaths()
plotOptions = list(unitTime = c("hour"), unitOutcomes = c("mcg/mL"))
reportOptions = list(unitTime = "", unitOutcomes = "")

showOutputMult = showOutputSimplex = NULL
plotMult_SE = plotMult_RSE = plotSimplex_SE = plotSimplex_RSE = NULL

modelEquations = list(
  "RespPK" = "dose_RespPK/V * ka/(ka - Cl/V) * (exp(-Cl/V * t) - exp(-ka * t))"
)
modelError = list(Constant(output = "RespPK", sigmaInter = 0.1))
modelParameters = list(
  ModelParameter(name = "ka", distribution = LogNormal(mu = 1,   omega = sqrt(0.09)), gamma = sqrt(0.0225)),
  ModelParameter(name = "V",  distribution = LogNormal(mu = 3.5, omega = sqrt(0.09)), gamma = sqrt(0.0225)),
  ModelParameter(name = "Cl", distribution = LogNormal(mu = 2,   omega = sqrt(0.09)), gamma = sqrt(0.0225))
)
sex = Covariate(
  name = "Sex",
  categories = c("M", "F"),
  categoriesProportions = c(0.5, 0.5),
  effects = list("F" = c("V" = log(1.2)))
)
treatment = Covariate(
  name = "Treatment",
  categories = c("R", "T"),
  sequences = list(c("R", "T"), c("T", "R")),
  sequencesProportions = c(0.5, 0.5),
  effects = list("T" = c("Cl" = log(1.1)))
)
administrationRespPK = Administration(outcome = "RespPK", timeDose = c(0), dose = c(30))
samplingTimesRespPK  = SamplingTimes(outcome = "RespPK", samplings = c(0.5, 2, 4, 6, 8))

# --- Evaluation ---
arm1    = Arm(name = "arm1", size = 40,
               administrations = list(administrationRespPK),
               samplingTimes   = list(samplingTimesRespPK))
design1 = Design(name = "design1", arms = list(arm1))

evaluation = Evaluation(
  name            = "Vignette03_Eval",
  modelEquations  = modelEquations,
  modelParameters = modelParameters,
  modelCovariates = list(sex, treatment),
  modelCovariatesEquation = "exponential",
  modelError      = modelError,
  designs         = list(design1),
  fimType         = "population",
  outputs         = list("RespPK" = "RespPK")
)
evaluation = run(evaluation)

showOutputEvaluation = pfimCapture({
  show(evaluation)
  getFisherMatrix(evaluation)
  getCorrelationMatrix(evaluation)
  getSE(evaluation)
  getRSE(evaluation)
  getDeterminant(evaluation)
  getDcriterion(evaluation)
})
writeLines(showOutputEvaluation,
           file.path(paths$outputs, "vignette3_evaluation_PopFIM_show.txt"))
invisible(pfimSafeReport(evaluation, paths$reports,
                         "vignette3_evaluation_popFim_report.html", reportOptions))

# --- Covariate test on the initial design ---
resultsTests = covariateTest(evaluation)
showOutputTests = pfimCapture(show(resultsTests))
writeLines(showOutputTests,
           file.path(paths$outputs, "vignette3_covariateTest_initialDesign_show.txt"))

plotsEval3_eval = plotEvaluation(evaluation, plotOptions)
plotsEval3_si   = plotSensitivityIndices(evaluation, plotOptions)
plotOutcomesEvaluationRespPK    = plotsEval3_eval[["design1"]][["arm1"]][["RespPK"]]
plotSensitivityIndice_RespPK_V  = plotsEval3_si[["design1"]][["arm1"]][["RespPK"]][["V"]]
plotSensitivityIndice_RespPK_Cl = plotsEval3_si[["design1"]][["arm1"]][["RespPK"]][["Cl"]]
plotEval_SE  = PFIM::plotSE(evaluation)
plotEval_RSE = PFIM::plotRSE(evaluation)

# --- Optimizations ---
finalizeOpt = function(opt, reportName, showName) {
  pfimFinalizeOpt(opt, paths, reportName, showName, reportOptions)
}

administrationConstraintsRespPK = AdministrationConstraints(
  outcome = "RespPK", doses = list(30))
samplingConstraintsRespPK = SamplingTimeConstraints(
  outcome = "RespPK",
  initialSamplings = c(0.5, 2, 4, 6, 8),
  numberOfsamplingsOptimisable = 3)
armMult = Arm(name = "armOpt", size = 40,
               administrations            = list(administrationRespPK),
               samplingTimes              = list(samplingTimesRespPK),
               administrationsConstraints = list(administrationConstraintsRespPK),
               samplingTimesConstraints   = list(samplingConstraintsRespPK))
designMult = Design(name = "design1", arms = list(armMult))

optimizationMult = pfimRunOrLoad(
  file.path(paths$data, "vignette3_optimization_Mult_populationFIM.RDS"),
  function() {
    opt = Optimization(
      name                    = "Multiplicative",
      modelEquations          = modelEquations,
      modelParameters         = modelParameters,
      modelCovariates         = list(treatment, sex),
      modelCovariatesEquation = "exponential",
      numberOfOccasions       = 2,
      modelError              = modelError,
      optimizer               = "MultiplicativeAlgorithm",
      optimizerParameters     = list(lambda             = 0.99,
                                     numberOfIterations = 1000,
                                     weightThreshold    = 0.01,
                                     delta              = 1e-04,
                                     showProcess        = FALSE),
      designs                 = list(designMult),
      fimType                 = "population",
      outputs                 = list("RespPK" = "RespPK"))
    run(opt)
  })
optimizationMult = finalizeOpt(
  optimizationMult,
  "vignette3_optimization_Mult_populationFIM_report.html",
  "vignette3_optimization_Mult_populationFIM_show.txt")
showOutputMult = attr(optimizationMult, "showOutput")
plotMult_SE  = pfimSafePlot(PFIM::plotSE(optimizationMult))
plotMult_RSE = pfimSafePlot(PFIM::plotRSE(optimizationMult))

samplingTimesRespPK_simplex = SamplingTimes(outcome = "RespPK", samplings = c(0.5, 4, 8))
samplingConstraintsRespPK_simplex = SamplingTimeConstraints(
  outcome = "RespPK",
  initialSamplings       = c(0.5, 4, 8),
  samplingsWindows       = list(c(0, 8)),
  numberOfTimesByWindows = c(3),
  minSampling            = c(0.5))
armSimplex = Arm(name = "armOpt", size = 40,
                  administrations          = list(administrationRespPK),
                  samplingTimes            = list(samplingTimesRespPK_simplex),
                  samplingTimesConstraints = list(samplingConstraintsRespPK_simplex))
designSimplex = Design(name = "design1", arms = list(armSimplex))

optimizationSimplex = pfimRunOrLoad(
  file.path(paths$data, "vignette3_optimization_Simplex_populationFIM.RDS"),
  function() {
    opt = Optimization(
      name                    = "Simplex",
      modelEquations          = modelEquations,
      modelParameters         = modelParameters,
      modelCovariates         = list(treatment, sex),
      modelCovariatesEquation = "exponential",
      modelError              = modelError,
      optimizer               = "SimplexAlgorithm",
      optimizerParameters     = list(pctInitialSimplexBuilding = 20,
                                     maxIteration              = 200,
                                     tolerance                 = 1e-6,
                                     showProcess               = FALSE),
      designs                 = list(designSimplex),
      fimType                 = "population",
      outputs                 = list("RespPK" = "RespPK"))
    run(opt)
  })
optimizationSimplex = finalizeOpt(
  optimizationSimplex,
  "vignette3_optimization_Simplex_populationFIM_report.html",
  "vignette3_optimization_Simplex_populationFIM_show.txt")
showOutputSimplex = attr(optimizationSimplex, "showOutput")
plotSimplex_SE  = pfimSafePlot(PFIM::plotSE(optimizationSimplex))
plotSimplex_RSE = pfimSafePlot(PFIM::plotRSE(optimizationSimplex))

# --- Covariate test on the optimal design (Simplex) ---
optimalTestsShowFile = file.path(paths$outputs,
                                  "vignette3_covariateTest_optimalDesign_show.txt")
showOutputOptimalTests = tryCatch({
  optimisationDesign      = prop(optimizationSimplex, "optimisationDesign")
  evaluationOptimalDesign = optimisationDesign$evaluationOptimalDesign
  optimalTests            = covariateTest(evaluationOptimalDesign)
  pfimCapture(show(optimalTests))
}, error = function(e) {
  if (file.exists(optimalTestsShowFile)) {
    paste(readLines(optimalTestsShowFile, warn = FALSE), collapse = "\n")
  } else {
    paste0(
      "[Covariate test on the optimal design unavailable for this RDS object.\n",
      "Run run() and saveRDS() with the current PFIM version to regenerate.\n",
      conditionMessage(e)
    )
  }
})
writeLines(showOutputOptimalTests, optimalTestsShowFile)

invisible(NULL)