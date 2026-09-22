# @keywords internal
#' @name Additive
#' @examples
#' Additive(beta = 1, combinedEffect = 0.2)
NULL

#' @name Administration
#' @examples
#' Administration(outcome = "RespPK", timeDose = 0, dose = 100)
NULL

#' @name AdministrationConstraints
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(AdministrationConstraints)
#' }
NULL

#' @name armAdministration
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(armAdministration)
#' }
NULL

#' @name CategoricalCovariate
#' @examples
#' CategoricalCovariate(
#'   name = "Sex", categories = c("M", "F"),
#'   categoriesProportions = c(0.5, 0.5)
#' )
NULL

#' @name CategoricalCovariateWithIOV
#' @examples
#' CategoricalCovariateWithIOV(
#'   name = "Trt",
#'   categories = c("A", "B"),
#'   sequences = list(c("A", "B"), c("B", "A")),
#'   sequencesProportions = c(0.5, 0.5)
#' )
NULL

#' @name checkSamplingTimeConstraintsForMetaheuristic
#' @examples
#' \dontrun{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(checkSamplingTimeConstraintsForMetaheuristic)
#' }
NULL

#' @name checkValiditySamplingConstraint
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(checkValiditySamplingConstraint)
#' }
NULL

#' @name Combined1
#' @examples
#' Combined1(output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.1)
NULL

#' @name Combined2
#' @examples
#' Combined2(output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.1)
NULL

#' @name computeCovariateValue
#' @examples
#' PFIM:::computeCovariateValue(Additive(), beta = 2, combinedEffect = 0.1)
NULL

#' @name Constant
#' @examples
#' Constant(output = "RespPK", sigmaInter = 0.1)
NULL

#' @name constraintsTableForReport
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(constraintsTableForReport)
#' }
NULL

#' @name convertPKModelAnalyticToPKModelODE
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(convertPKModelAnalyticToPKModelODE)
#' }
NULL

#' @name Covariate
#' @examples
#' Covariate(name = "Sex", categories = c("M", "F"),
#' categoriesProportions = c(0.5, 0.5), effects = list(F = c(V = log(1.2))))
NULL

#' @name CovariateModelEquation
#' @examples
#' CovariateModelEquation(beta = 1, combinedEffect = 0.2)
NULL

#' @name createEffectVector
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(createEffectVector)
#' }
NULL

#' @name defineCovariatesData
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(defineCovariatesData)
#' }
NULL

#' @name defineFim
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' defineFim(ev)
#' }
NULL

#' @name defineModelAdministration
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(defineModelAdministration)
#' }
NULL

#' @name defineModelType
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(defineModelType)
#' }
NULL

#' @name defineModelWrapper
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(defineModelWrapper)
#' }
NULL

#' @name defineOptimizationAlgorithm
#' @examples
#' \donttest{
#' # help(defineOptimizationAlgorithm)
#' }
NULL

#' @name Distribution
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(Distribution)
#' }
NULL

#' @name evaluateArm
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(evaluateArm)
#' }
NULL

#' @name evaluateCovariatesEffects
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(evaluateCovariatesEffects)
#' }
NULL

#' @name evaluateErrorModelDerivatives
#' @examples
#' err = Combined2(output = "RespPK", sigmaInter = 0.2, sigmaSlope = 0.1)
#' PFIM:::evaluateErrorModelDerivatives(err, c(1, 2))
NULL

#' @name evaluateInitialConditions
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(evaluateInitialConditions)
#' }
NULL

#' @name evaluateModel
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(evaluateModel)
#' }
NULL

#' @name evaluateOmegaMatrixFromCovariates
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(evaluateOmegaMatrixFromCovariates)
#' }
NULL

#' @name Exponential
#' @examples
#' Exponential(beta = 1, combinedEffect = 0.2)
NULL

#' @name Fim
#' @examples
#' PopulationFim()
NULL

#' @name finiteDifferenceHessian
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(finiteDifferenceHessian)
#' }
NULL

#' @name generateCovariatesCombination
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(generateCovariatesCombination)
#' }
NULL

#' @name generateDosesCombination
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(generateDosesCombination)
#' }
NULL

#' @name generateFimsFromConstraints
#' @examples
#' \dontrun{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(generateFimsFromConstraints)
#' }
NULL

#' @name generateReportEvaluation
#' @examples
#' \dontrun{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # Report(ev, outputPath = tempdir())
#' }
NULL

#' @name generateReportOptimization
#' @examples
#' \dontrun{
#' \dontrun{vignette("Example01")}
#' }
NULL

#' @name generateSamplingsFromSamplingConstraints
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(generateSamplingsFromSamplingConstraints)
#' }
NULL

#' @name generateSamplingTimesCombination
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(generateSamplingTimesCombination)
#' }
NULL

#' @name getArmConstraints
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getArmConstraints)
#' }
NULL

#' @name getArmData
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getArmData)
#' }
NULL

#' @name getCategoryOfReference
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getCategoryOfReference)
#' }
NULL

#' @name getCorrelationMatrix
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getCorrelationMatrix)
#' }
NULL

#' @name getCovariateEffects
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getCovariateEffects)
#' }
NULL

#' @name getCovariateTestTables
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getCovariateTestTables)
#' }
NULL

#' @name getDcriterion
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getDcriterion(ev)
#' }
NULL

#' @name getDeterminant
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getDeterminant(ev)
#' }
NULL

#' @name getEvaluationDesign
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getEvaluationDesign)
#' }
NULL

#' @name getFim
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getFim(ev)
#' }
NULL

#' @name getFisherMatrix
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getFisherMatrix(ev)
#' }
NULL

#' @name getModelErrorData
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getModelErrorData)
#' }
NULL

#' @name getModelParametersData
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getModelParametersData)
#' }
NULL

#' @name getRSE
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getRSE(ev)
#' }
NULL

#' @name getSamplingData
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getSamplingData)
#' }
NULL

#' @name getSE
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getSE(ev)
#' }
NULL

#' @name getShrinkage
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(getShrinkage)
#' }
NULL

#' @name hasCovariates
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(hasCovariates)
#' }
NULL

#' @name Linear2BolusSingleDose_ClQV1V2
#' @examples
#' nzchar(Linear2BolusSingleDose_ClQV1V2())
NULL

#' @name Linear2BolusSingleDose_kk12k21V
#' @examples
#' nzchar(Linear2BolusSingleDose_kk12k21V())
NULL

#' @name Linear2BolusSteadyState_ClQV1V2tau
#' @examples
#' nzchar(Linear2BolusSteadyState_ClQV1V2tau())
NULL

#' @name Linear2BolusSteadyState_kk12k21Vtau
#' @examples
#' nzchar(Linear2BolusSteadyState_kk12k21Vtau())
NULL

#' @name Linear2FirstOrderSingleDose_kaClQV1V2
#' @examples
#' nzchar(Linear2FirstOrderSingleDose_kaClQV1V2())
NULL

#' @name Linear2FirstOrderSingleDose_kakk12k21V
#' @examples
#' nzchar(Linear2FirstOrderSingleDose_kakk12k21V())
NULL

#' @name Linear2FirstOrderSteadyState_kaClQV1V2tau
#' @examples
#' nzchar(Linear2FirstOrderSteadyState_kaClQV1V2tau())
NULL

#' @name Linear2FirstOrderSteadyState_kakk12k21Vtau
#' @examples
#' nzchar(Linear2FirstOrderSteadyState_kakk12k21Vtau())
NULL

#' @name Linear2InfusionSingleDose_ClQV1V2
#' @examples
#' nzchar(Linear2InfusionSingleDose_ClQV1V2()$duringInfusion$RespPK)
NULL

#' @name Linear2InfusionSingleDose_kk12k21V
#' @examples
#' nzchar(Linear2InfusionSingleDose_kk12k21V()$duringInfusion$RespPK)
NULL

#' @name Linear2InfusionSteadyState_ClQV1V2tau
#' @examples
#' nzchar(Linear2InfusionSteadyState_ClQV1V2tau()$duringInfusion$RespPK)
NULL

#' @name Linear2InfusionSteadyState_kk12k21Vtau
#' @examples
#' nzchar(Linear2InfusionSteadyState_kk12k21Vtau()$duringInfusion$RespPK)
NULL

#' @name LogNormal
#' @examples
#' LogNormal(mu = 1, omega = 0.3)
NULL

#' @name Model
#' @examples
#' Model()
NULL

#' @name ModelAnalytic
#' @examples
#' ModelAnalytic()
NULL

#' @name ModelAnalyticInfusion
#' @examples
#' ModelAnalyticInfusion()
NULL

#' @name ModelAnalyticInfusionSteadyState
#' @examples
#' ModelAnalyticInfusionSteadyState()
NULL

#' @name ModelAnalyticSteadyState
#' @examples
#' ModelAnalyticSteadyState()
NULL

#' @name ModelError
#' @examples
#' ModelError()
NULL

#' @name ModelInfusion
#' @examples
#' ModelInfusion()
NULL

#' @name ModelODE
#' @examples
#' ModelODE()
NULL

#' @name ModelODEBolus
#' @examples
#' ModelODEBolus()
NULL

#' @name ModelODEDoseInEquations
#' @examples
#' ModelODEDoseInEquations()
NULL

#' @name ModelODEDoseNotInEquations
#' @examples
#' ModelODEDoseNotInEquations()
NULL

#' @name ModelODEInfusion
#' @examples
#' ModelODEInfusion()
NULL

#' @name ModelODEInfusionDoseInEquation
#' @examples
#' ModelODEInfusionDoseInEquation()
NULL

#' @name ModelParameter
#' @examples
#' ModelParameter(name = "Cl", distribution = LogNormal(mu = 1, omega = 0.3))
NULL

#' @name modelParametersWithCovariates
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(modelParametersWithCovariates)
#' }
NULL

#' @name Normal
#' @examples
#' Normal(mu = 0, omega = 1)
NULL

#' @name optimizeDesign
#' @examples
#' \dontrun{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(optimizeDesign)
#' }
NULL

#' @name pfim_cache_stats
#' @examples
#' PFIM:::pfim_cache_stats()
NULL

#' @name pfim_get_option
#' @examples
#' pfim_get_option("fim.cache.maxEntries")
NULL

#' @name pfim_list_extensions
#' @examples
#' pfim_list_extensions()
NULL

#' @name pfim_register_fim_type
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(pfim_register_fim_type)
#' }
NULL

#' @name pfim_register_model_class
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(pfim_register_model_class)
#' }
NULL

#' @name pfim_register_optimizer
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(pfim_register_optimizer)
#' }
NULL

#' @name pfim_reset_session
#' @examples
#' pfim_reset_session()
NULL

#' @name pfim_resolve_model_class
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(pfim_resolve_model_class)
#' }
NULL

#' @name pfim_set_option
#' @examples
#' pfim_set_option(fim.cache.maxEntries = 2048L)
NULL

#' @name PFIMProject
#' @examples
#' PFIMProject(name = "demo")
NULL

#' @name plotFrequencies
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(plotFrequencies)
#' }
NULL

#' @name plotShrinkage
#' @examples
#' \dontrun{
#' source(system.file("examples", "covariate-test-minimal.R", package = "PFIM"))
#' # plotShrinkage(prop(ev, "fim"), ev)  # BayesianFim
#' }
NULL

#' @name plotWeights
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(plotWeights)
#' }
NULL

#' @name plotWeightsMultiplicativeAlgorithm
#' @examples
#' \dontrun{
#' # help(plotWeightsMultiplicativeAlgorithm)
#' }
NULL

#' @name projectOf
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' projectOf(ev)
#' }
NULL

#' @name Proportional
#' @examples
#' Proportional(output = "RespPK", sigmaSlope = 0.15)
NULL

#' @name replaceVariablesLibraryOfModels
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(replaceVariablesLibraryOfModels)
#' }
NULL

#' @name SamplingTimeConstraints
#' @examples
#' SamplingTimeConstraints(
#'   outcome = "RespPK",
#'   initialSamplings = c(0.5, 2, 8),
#'   numberOfsamplingsOptimisable = 2
#' )
NULL

#' @name SamplingTimes
#' @examples
#' SamplingTimes(outcome = "RespPK", samplings = c(1, 4, 8))
NULL

#' @name saveCovariateTest
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(saveCovariateTest)
#' }
NULL

#' @name setOptimalArms
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(setOptimalArms)
#' }
NULL

#' @name setSamplingConstraintForOptimization
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(setSamplingConstraintForOptimization)
#' }
NULL

#' @name tablesForReport
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(tablesForReport)
#' }
NULL

#' @name updateSamplingTimes
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' # help(updateSamplingTimes)
#' }
NULL

