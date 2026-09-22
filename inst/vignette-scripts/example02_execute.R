# Example02: evaluation, optimization, plots, and reports.
utils = system.file("vignette-scripts", "pfim-vignette-utils.R", package = "PFIM")
if (!nzchar(utils)) stop("pfim-vignette-utils.R not found.", call. = FALSE)
source(utils, local = TRUE)
library(PFIM)

paths = pfimVignetteSetupPaths()
plotOptions = list(unitTime = c("hour"), unitOutcomes = c("mcg/mL"))

showOutputPSO = showOutputPGBO = showOutputSimplex = NULL
plotPSO_SE = plotPSO_RSE = plotPGBO_SE = plotPGBO_RSE = NULL
plotSimplex_SE = plotSimplex_RSE = NULL

modelFromLibrary = list("PKModel" = "Linear1InfusionSingleDose_ClV")
modelParameters = list(
  ModelParameter(name = "V",  distribution = LogNormal(mu = 50, omega = sqrt(0.26))),
  ModelParameter(name = "Cl", distribution = LogNormal(mu = 5,  omega = sqrt(0.34)))
)
modelError = list(Combined1(output = "RespPK", sigmaInter = 0.5, sigmaSlope = sqrt(0.15)))
administrationRespPK = Administration(
  outcome = "RespPK", Tinf = rep(1, 5),
  timeDose = seq(0, 96, 24), dose = c(400, rep(200, 4))
)
samplingTimesRespPK = SamplingTimes(outcome = "RespPK", samplings = c(1, 12, 24, 44, 72, 120))
arm1 = Arm(name = "arm1", size = 150,
            administrations = list(administrationRespPK),
            samplingTimes = list(samplingTimesRespPK))
design1 = Design(name = "design1", arms = list(arm1))
evalCommon = list(
  modelFromLibrary = modelFromLibrary, modelParameters = modelParameters,
  modelError = modelError, outputs = list("RespPK"), designs = list(design1),
  odeSolverParameters = list(atol = 1e-8, rtol = 1e-8)
)

evaluationPopResults = do.call(Evaluation, c(list(name = "evaluation", fimType = "population"), evalCommon))
evaluationPopResults = run(evaluationPopResults)
evaluationFIMPop = evaluationPopResults
showOutputFIMPop = pfimCapture({
  show(evaluationPopResults)
  getFisherMatrix(evaluationPopResults)
  getCorrelationMatrix(evaluationPopResults)
  getSE(evaluationPopResults)
  getRSE(evaluationPopResults)
  getShrinkage(evaluationPopResults)
  getDeterminant(evaluationPopResults)
  getDcriterion(evaluationPopResults)
})
writeLines(showOutputFIMPop,
           file.path(paths$outputs, "vignette2_evaluation_PopFIM_show.txt"))
invisible(pfimSafeReport(evaluationPopResults, paths$reports, "vignette2_evaluation_popFim_report.html", plotOptions))

evaluationIndResults = do.call(Evaluation, c(list(name = "evaluation", fimType = "individual"), evalCommon))
evaluationIndResults = run(evaluationIndResults)
evaluationFIMInd = evaluationIndResults
showOutputFIMInd = pfimCapture(show(evaluationIndResults))
writeLines(showOutputFIMInd,
           file.path(paths$outputs, "vignette2_evaluation_IndFIM_show.txt"))
invisible(pfimSafeReport(evaluationIndResults, paths$reports, "vignette2_evaluation_indFim_report.html", plotOptions))

evaluationBayResults = do.call(Evaluation, c(list(name = "evaluation", fimType = "Bayesian"), evalCommon))
evaluationBayResults = run(evaluationBayResults)
evaluationFIMBay = evaluationBayResults
showOutputFIMBay = pfimCapture(show(evaluationBayResults))
writeLines(showOutputFIMBay,
           file.path(paths$outputs, "vignette2_evaluation_BayFIM_show.txt"))
invisible(pfimSafeReport(evaluationBayResults, paths$reports, "vignette2_evaluation_BayFim_report.html", plotOptions))

plotsEval2_eval = plotEvaluation(evaluationPopResults, plotOptions)
plotsEval2_si   = plotSensitivityIndices(evaluationPopResults, plotOptions)
plotEval_SE = PFIM::plotSE(evaluationPopResults)
plotEval_RSE = PFIM::plotRSE(evaluationPopResults)
plotOutcomesEvaluationRespPK = plotsEval2_eval[["design1"]][["arm1"]][["RespPK"]]
plotSensitivityIndice_RespPK_V = plotsEval2_si[["design1"]][["arm1"]][["RespPK"]][["V"]]
plotSensitivityIndice_RespPK_Cl = plotsEval2_si[["design1"]][["arm1"]][["RespPK"]][["Cl"]]

# --- Optimisations ---
samplingTimesRespPK_opt = SamplingTimes(outcome = "RespPK", samplings = c(1, 48, 72, 120))
samplingConstraintsRespPK = SamplingTimeConstraints(
  outcome = "RespPK", initialSamplings = c(1, 48, 72, 120),
  numberOfTimesByWindows = c(2, 2),
  samplingsWindows = list(c(1, 48), c(72, 120)), minSampling = 5)
arm2 = Arm(name = "arm2", size = 150,
            administrations = list(administrationRespPK),
            samplingTimes = list(samplingTimesRespPK_opt),
            samplingTimesConstraints = list(samplingConstraintsRespPK))
design2 = Design(name = "design2", arms = list(arm2), numberOfArms = 150)

optCommon = list(
  name = "", modelFromLibrary = modelFromLibrary,
  modelParameters = modelParameters, modelError = modelError,
  designs = list(design2), fimType = "population", outputs = list("RespPK")
)

finalizeOpt = function(opt, reportName, showName) {
  pfimFinalizeOpt(opt, paths, reportName, showName, plotOptions)
}

optimizationPSO = pfimRunOrLoad(
  file.path(paths$data, "vignette2_optimization_PSO_populationFIM.RDS"),
  function() {
    opt = do.call(Optimization, c(optCommon, list(
      optimizer = "PSOAlgorithm",
      optimizerParameters = list(
        maxIteration = 100, populationSize = 50,
        personalLearningCoefficient = 2.05, globalLearningCoefficient = 2.05,
        seed = 42, showProcess = FALSE))))
    run(opt)
  })
optimizationPSO = finalizeOpt(
  optimizationPSO,
  "vignette2_optimization_PSO_populationFIM_report.html",
  "vignette2_optimization_PSO_populationFIM_show.txt")
showOutputPSO = attr(optimizationPSO, "showOutput")
plotPSO_SE = pfimSafePlot(PFIM::plotSE(optimizationPSO))
plotPSO_RSE = pfimSafePlot(PFIM::plotRSE(optimizationPSO))

optimizationPGBO = pfimRunOrLoad(
  file.path(paths$data, "vignette2_optimization_PGBO_populationFIM.RDS"),
  function() {
    opt = do.call(Optimization, c(optCommon, list(
      optimizer = "PGBOAlgorithm",
      optimizerParameters = list(
        N = 30, muteEffect = 0.65, maxIteration = 1000,
        purgeIteration = 200, seed = 42, showProcess = FALSE))))
    run(opt)
  })
optimizationPGBO = finalizeOpt(
  optimizationPGBO,
  "vignette2_optimization_PGBO_populationFIM_report.html",
  "vignette2_optimization_PGBO_populationFIM_show.txt")
showOutputPGBO = attr(optimizationPGBO, "showOutput")
plotPGBO_SE = pfimSafePlot(PFIM::plotSE(optimizationPGBO))
plotPGBO_RSE = pfimSafePlot(PFIM::plotRSE(optimizationPGBO))

optimizationSimplex = pfimRunOrLoad(
  file.path(paths$data, "vignette2_optimization_Simplex_populationFIM.RDS"),
  function() {
    opt = do.call(Optimization, c(optCommon, list(
      optimizer = "SimplexAlgorithm",
      optimizerParameters = list(
        pctInitialSimplexBuilding = 10, maxIteration = 1000,
        tolerance = 1e-10, showProcess = FALSE))))
    run(opt)
  })
optimizationSimplex = finalizeOpt(
  optimizationSimplex,
  "vignette2_optimization_Simplex_populationFIM_report.html",
  "vignette2_optimization_Simplex_populationFIM_show.txt")
showOutputSimplex = attr(optimizationSimplex, "showOutput")
plotSimplex_SE = pfimSafePlot(PFIM::plotSE(optimizationSimplex))
plotSimplex_RSE = pfimSafePlot(PFIM::plotRSE(optimizationSimplex))

invisible(NULL)
