# Example01: evaluation, optimization, plots, and reports.
utils = system.file("vignette-scripts", "pfim-vignette-utils.R", package = "PFIM")
if (!nzchar(utils)) stop("pfim-vignette-utils.R not found.", call. = FALSE)
source(utils, local = TRUE)
library(PFIM)

paths = pfimVignetteSetupPaths()
plotOptions = list(unitTime = c("hour"), unitOutcomes = c("mcg/mL", "DI%"))

showOutputFW = showOutputMult = NULL
plotFWFrequencies = plotFW_SE = plotFW_RSE = NULL
plotMultWeights = plotMult_SE = plotMult_RSE = NULL

modelEquations = list(
  "Deriv_Cc" = "dose_RespPK/V*ka*exp(-ka*t) - Cl/V*Cc",
  "Deriv_E"  = "Rin*(1-Imax*(Cc**gamma)/(Cc**gamma + IC50**gamma))-kout*E"
)
modelParameters = list(
  ModelParameter(name = "V",     distribution = LogNormal(mu = 0.74,  omega = 0.316)),
  ModelParameter(name = "Cl",    distribution = LogNormal(mu = 0.28,  omega = 0.456)),
  ModelParameter(name = "ka",    distribution = LogNormal(mu = 10,    omega = 0), fixedMu = TRUE),
  ModelParameter(name = "kout",  distribution = LogNormal(mu = 6.14,  omega = 0.947)),
  ModelParameter(name = "Rin",   distribution = LogNormal(mu = 614,   omega = 0), fixedMu = TRUE),
  ModelParameter(name = "Imax",  distribution = LogNormal(mu = 0.76,  omega = 0.439)),
  ModelParameter(name = "IC50",  distribution = LogNormal(mu = 9.22,  omega = 0.452)),
  ModelParameter(name = "gamma", distribution = LogNormal(mu = 2.77,  omega = 1.761))
)
modelError = list(
  Combined1(output = "RespPK", sigmaInter = 0,   sigmaSlope = 0.21),
  Constant(output = "RespPD", sigmaInter = 9.6)
)
samplingTimesRespPK = SamplingTimes(outcome = "RespPK",
  samplings = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4))
samplingTimesRespPD = SamplingTimes(outcome = "RespPD",
  samplings = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4))
mkArm = function(name, dose) {
  Arm(name = name, size = 6,
      administrations = list(Administration(outcome = "RespPK", timeDose = 0, dose = dose)),
      samplingTimes = list(samplingTimesRespPK, samplingTimesRespPD),
      initialCondition = list("Cc" = 0, "E" = 100))
}
design1 = Design(name = "design1", arms = list(
  mkArm("0.2mg Arm", 0.2), mkArm("0.64mg Arm", 0.64), mkArm("2mg Arm", 2),
  mkArm("6.32mg Arm", 6.32), mkArm("11.24mg Arm", 11.24), mkArm("20mg Arm", 20)
))
evalCommon = list(
  modelEquations = modelEquations, modelParameters = modelParameters,
  modelError = modelError,
  outputs = list("RespPK" = "Cc", "RespPD" = "E"),
  designs = list(design1),
  odeSolverParameters = list(atol = 1e-8, rtol = 1e-8)
)

evaluationPopResults = do.call(Evaluation, c(list(name = "evaluation", fimType = "population"), evalCommon))
evaluationPopResults = run(evaluationPopResults)
evaluationPop = evaluationPopResults
showOutputEvaluationPop = pfimCapture({
  show(evaluationPopResults)
  getFisherMatrix(evaluationPopResults)
  getCorrelationMatrix(evaluationPopResults)
  getSE(evaluationPopResults)
  getRSE(evaluationPopResults)
  getShrinkage(evaluationPopResults)
  getDeterminant(evaluationPopResults)
  getDcriterion(evaluationPopResults)
})
writeLines(showOutputEvaluationPop,
           file.path(paths$outputs, "vignette1_evaluation_populationFIM_show.txt"))
invisible(pfimSafeReport(evaluationPopResults, paths$reports, "vignette1_evaluation_popFIM_report.html", plotOptions))

evaluationBayesianResults = do.call(Evaluation, c(list(name = "evaluation", fimType = "Bayesian"), evalCommon))
evaluationBayesianResults = run(evaluationBayesianResults)
evaluationBay = evaluationBayesianResults
showOutputEvaluationBay = pfimCapture({
  show(evaluationBayesianResults)
  getFisherMatrix(evaluationBayesianResults)
  getCorrelationMatrix(evaluationBayesianResults)
  getSE(evaluationBayesianResults)
  getRSE(evaluationBayesianResults)
  getShrinkage(evaluationBayesianResults)
  getDeterminant(evaluationBayesianResults)
  getDcriterion(evaluationBayesianResults)
})
writeLines(showOutputEvaluationBay,
           file.path(paths$outputs, "vignette1_evaluation_BayesianFIM_show.txt"))

plotsEval1_eval = plotEvaluation(evaluationPopResults, plotOptions)
plotsEval1_si   = plotSensitivityIndices(evaluationPopResults, plotOptions)
plotEval_SE = PFIM::plotSE(evaluationPopResults)
plotEval_RSE = PFIM::plotRSE(evaluationPopResults)
plotOutcomesEvaluationRespPK = plotsEval1_eval[["design1"]][["20mg Arm"]][["RespPK"]]
plotOutcomesEvaluationRespPD = plotsEval1_eval[["design1"]][["20mg Arm"]][["RespPD"]]
plotSensitivityIndice_RespPK_Cl = plotsEval1_si[["design1"]][["20mg Arm"]][["RespPK"]][["Cl"]]
plotSensitivityIndice_RespPK_V = plotsEval1_si[["design1"]][["20mg Arm"]][["RespPK"]][["V"]]

# --- Optimizations ---
pfim_set_option(constraints.maxTasks = 500)

administrationRespPK = Administration(outcome = "RespPK", timeDose = c(0), dose = c(6.32))
samplingTimesRespPK_opt = SamplingTimes(
  outcome = "RespPK", samplings = c(0.25, 0.75, 1, 1.5, 2, 4, 6))
samplingTimesRespPD_opt = SamplingTimes(
  outcome = "RespPD", samplings = c(0.25, 0.75, 1.5, 2, 3, 6, 8, 12))
samplingConstraintsRespPK = SamplingTimeConstraints(
  outcome = "RespPK",
  initialSamplings = c(0.25, 0.75, 1, 1.5, 2, 4, 6),
  fixedTimes = c(0.25, 4), numberOfsamplingsOptimisable = 4)
samplingConstraintsRespPD = SamplingTimeConstraints(
  outcome = "RespPD",
  initialSamplings = c(0.25, 0.75, 1.5, 2, 3, 6, 8, 12),
  fixedTimes = c(2, 6), numberOfsamplingsOptimisable = 4)
initialElementaryProtocols = list(
  list(
    c(0.25, 0.75, 1, 4),
    c(1.5, 2, 6, 12)
  )
)
administrationConstraintsRespPK = AdministrationConstraints(
  outcome = "RespPK", doses = list(0.2, 0.64, 2, 6.32, 11.24, 20))
armConstraint = Arm(
  name = "armConstraint", size = 30,
  administrations = list(administrationRespPK),
  samplingTimes = list(samplingTimesRespPK_opt, samplingTimesRespPD_opt),
  administrationsConstraints = list(administrationConstraintsRespPK),
  samplingTimesConstraints = list(samplingConstraintsRespPK, samplingConstraintsRespPD),
  initialCondition = list("Cc" = 0, "E" = "Rin/kout"))
designConstraint = Design(name = "designConstraint", arms = list(armConstraint), numberOfArms = 30)
numberOfSubjects = c(30)
proportionsOfSubjects = c(30) / 30

optCommon = list(
  name = "PKPD_ODE_multi_doses_populationFIM",
  modelEquations = modelEquations, modelParameters = modelParameters,
  modelError = modelError, designs = list(designConstraint),
  fimType = "population",
  outputs = list("RespPK" = "Cc", "RespPD" = "E"),
  odeSolverParameters = list(atol = 1e-8, rtol = 1e-8)
)

finalizeOpt = function(opt, reportName, showName) {
  pfimFinalizeOpt(opt, paths, reportName, showName, plotOptions)
}

filePathFW = file.path(paths$data, "vignette_1_optimization_FedorovWynn_populationFIM.RDS")
optimizationFWPopFIM = pfimRunOrLoad(filePathFW, function() {
  opt = do.call(Optimization, c(optCommon, list(
    optimizer = "FedorovWynnAlgorithm",
    optimizerParameters = list(
      elementaryProtocols = initialElementaryProtocols,
      numberOfSubjects = numberOfSubjects,
      proportionsOfSubjects = proportionsOfSubjects,
      showProcess = FALSE))))
  run(opt)
})
optimizationFWPopFIM = finalizeOpt(
  optimizationFWPopFIM,
  "vignette1_optimization_FedorovWynn_populationFIM_report.html",
  "vignette1_optimization_FedorovWynn_populationFIM_show.txt")
showOutputFW = attr(optimizationFWPopFIM, "showOutput")
plotFWFrequencies = pfimSafePlot(PFIM::plotFrequencies(optimizationFWPopFIM))
plotFW_SE = pfimSafePlot(PFIM::plotSE(optimizationFWPopFIM))
plotFW_RSE = pfimSafePlot(PFIM::plotRSE(optimizationFWPopFIM))

filePathMult = file.path(paths$data, "vignette_1_optimization_multiplicativeAlgorithm_populationFIM.RDS")
optimizationMultPopFIM = pfimRunOrLoad(filePathMult, function() {
  opt = do.call(Optimization, c(optCommon, list(
    optimizer = "MultiplicativeAlgorithm",
    optimizerParameters = list(
      lambda = 0.99, numberOfIterations = 1000,
      weightThreshold = 0.01, delta = 1e-04, showProcess = FALSE))))
  run(opt)
})
optimizationMultPopFIM = finalizeOpt(
  optimizationMultPopFIM,
  "vignette1_optimization_MultiplicativeAlgorithm_populationFIM_report.html",
  "vignette1_optimization_MultiplicativeAlgorithm_populationFIM_show.txt")
showOutputMult = attr(optimizationMultPopFIM, "showOutput")
plotMultWeights = pfimSafePlot(PFIM::plotWeights(optimizationMultPopFIM))
plotMult_SE = pfimSafePlot(PFIM::plotSE(optimizationMultPopFIM))
plotMult_RSE = pfimSafePlot(PFIM::plotRSE(optimizationMultPopFIM))

invisible(NULL)
