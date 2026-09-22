#' @title ModelAnalyticInfusionSteadyState
#' @description Analytic infusion model at steady state.
#' @inheritParams ModelInfusion
#' @param wrapperModelAnalyticInfusion                   Wrapper for the analytic solver.
#' @param functionArgumentsModelAnalyticInfusion         A list with the function arguments.
#' @param functionArgumentsSymbolModelAnalyticInfusion   A list with the function argument symbols.
#' @param solverInputs                                   A list with the solver inputs.
#' @inheritParams Model
#' @include ModelInfusion.R
#' @include ModelAnalytic.R
#' @return An S7 object of class \code{ModelAnalyticInfusionSteadyState}.
#' @export

ModelAnalyticInfusionSteadyState = new_class( "ModelAnalyticInfusionSteadyState",
                                              package = "PFIM",
                                              parent  = ModelInfusion,
                                              properties = list(
                                                wrapperModelAnalyticInfusion                 = new_property(class_list, default = list()),
                                                functionArgumentsModelAnalyticInfusion       = new_property(class_list, default = list()),
                                                functionArgumentsSymbolModelAnalyticInfusion = new_property(class_list, default = list()),
                                                solverInputs                                 = new_property(class_list, default = list())
                                              ))


#' Compile analytic wrappers for during/after infusion (with and without admin).
#'
#' Splits library equations into administered vs non-administered outcomes, then
#' builds four R functions via \code{.buildAnalyticWrapper}. Steady state adds
#' \code{tau} to the shared formal argument list.
#' @return Updated model with compiled infusion steady-state wrappers.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelAnalyticInfusionSteadyState ) = function( model, evaluation ) {

  # Outcomes that receive dosing vs library names used in equation keys.
  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )
  libraryOutcomesWithAdmin   = .administeredLibraryOutcomeNames( evaluation )

  # Formal names injected into each analytic wrapper: dose_*, Tinf_*, t_*, params, tau.
  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNames      = paste0( "t_",    outcomesWithAdministration )
  TinfNames      = paste0( "Tinf_", outcomesWithAdministration )

  # Split during/after equations into administered vs passive outcomes.
  equations                  = prop( evaluation, "modelEquations" )
  equationsDuringWithAdmin   = equations$duringInfusion[  names( equations$duringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithAdmin    = equations$afterInfusion[   names( equations$afterInfusion  ) %in% libraryOutcomesWithAdmin ]
  equationsDuringWithNoAdmin = equations$duringInfusion[ !names( equations$duringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithNoAdmin  = equations$afterInfusion[  !names( equations$afterInfusion  ) %in% libraryOutcomesWithAdmin ]

  outputDuringAdmin   = names( equationsDuringWithAdmin )
  outputDuringNoAdmin = names( equationsDuringWithNoAdmin )
  outputAfterAdmin    = names( equationsAfterWithAdmin )
  outputAfterNoAdmin  = names( equationsAfterWithNoAdmin )
  timeNamesNoAdmin = unique( c(
    .libraryEquationTimeNames( evaluation, outputDuringNoAdmin ),
    .libraryEquationTimeNames( evaluation, outputAfterNoAdmin )
  ) )

  # Steady-state dosing interval "tau" is required by the closed-form expressions.
  functionArgumentsAdmin   = unique( c( doseNames, TinfNames, parameterNames, timeNames, "tau" ) )
  functionArgumentsNoAdmin = unique( c( outcomesWithAdministration, parameterNames, timeNamesNoAdmin, "tau" ) )

  # Four evaluators: (during|after) x (administered|passive) outcomes.
  prop( model, "wrapperModelAnalyticInfusion" ) = list(
    functionDefinitionDuringInfusionWithAdmin   = .buildAnalyticWrapper(
      equationsDuringWithAdmin, functionArgumentsAdmin,
      .libraryEquationTimeNames( evaluation, names( equationsDuringWithAdmin ) ), outputDuringAdmin ),
    functionDefinitionDuringInfusionWithNoAdmin = .buildAnalyticWrapper(
      equationsDuringWithNoAdmin, functionArgumentsNoAdmin,
      .libraryEquationTimeNames( evaluation, names( equationsDuringWithNoAdmin ) ), outputDuringNoAdmin ),
    functionDefinitionAfterInfusionWithAdmin    = .buildAnalyticWrapper(
      equationsAfterWithAdmin, functionArgumentsAdmin,
      .libraryEquationTimeNames( evaluation, names( equationsAfterWithAdmin ) ), outputAfterAdmin ),
    functionDefinitionAfterInfusionWithNoAdmin  = .buildAnalyticWrapper(
      equationsAfterWithNoAdmin, functionArgumentsNoAdmin,
      .libraryEquationTimeNames( evaluation, names( equationsAfterWithNoAdmin ) ), outputAfterNoAdmin )
  )
  prop( model, "functionArgumentsModelAnalyticInfusion" ) = list(
    functionArguments = functionArgumentsAdmin,
    functionArgumentsNoAdmin = functionArgumentsNoAdmin
  )
  prop( model, "functionArgumentsSymbolModelAnalyticInfusion" ) = list(
    functionArgumentsSymbol = map( functionArgumentsAdmin, as.symbol )
  )
  prop( model, "outputNames" ) = unlist( names( equations$duringInfusion ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


#' Build per-outcome dose / infusion / sampling tables for the analytic solver.
#'
#' For \code{tau != 0}, expands a single steady-state dose into a regular grid
#' \code{0, tau, 2*tau, ...} up to the last sampling time. Labels each sampling
#' as \code{duringInfusion} or \code{afterInfusion}, and records which dose index
#' is active at that time (for superposition of past infusions).
#' @return Updated model with \code{samplings} and \code{solverInputs}.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelAnalyticInfusionSteadyState ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  samplings                  = .analyticSamplingGrid( arm )

  solverInputs = map( administrations, function( adm ) {
    dosing = .alignAdministrationDosing( adm )
    .analyticInfusionWindowTable(
      samplings, dosing$timeDose, dosing$dose, dosing$Tinf, prop( adm, "tau" )
    )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}


#' Evaluate analytic infusion steady-state concentrations at all sampling times.
#'
#' For each observation time and administered outcome:
#' \enumerate{
#'   \item If inside an infusion window: evaluate the during-infusion formula for
#'     the current dose; if earlier doses exist, add their after-infusion remnants.
#'   \item If after an infusion: evaluate the after-infusion formula for the
#'     current dose, again superposing remnants of previous doses.
#' }
#' Passive (non-admin) equations are evaluated once the administered outcome value
#' is assigned into the shared evaluation environment.
#' @param model A \code{ModelAnalyticInfusionSteadyState} object.
#' @param arm   An \code{Arm} object.
#' @return Named list of output data frames at requested sampling times.
#' @keywords internal
evaluateAnalyticInfusionSteadyStateCore = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = map_chr( administrations, ~ prop( .x, "outcome" ) )
  outputNames                = prop( model, "outputNames" ) |> unlist()
  solverInputs               = prop( model, "solverInputs" )
  samplings                  = prop( model, "samplings" )

  raw = prop( model, "wrapperModelAnalyticInfusion" )
  wrappers = .analyticBindUserEnv(
    raw$functionDefinitionDuringInfusionWithAdmin,
    raw$functionDefinitionDuringInfusionWithNoAdmin,
    raw$functionDefinitionAfterInfusionWithAdmin,
    raw$functionDefinitionAfterInfusionWithNoAdmin
  )
  fnDuringAdmin  = wrappers[[ 1L ]]
  fnAfterAdmin   = wrappers[[ 3L ]]
  fnAfterNoAdmin = wrappers[[ 4L ]]
  mu = .extractMu( prop( model, "modelParameters" ) )

  tmp = .analyticEvalGrid(
    samplings, outcomesWithAdministration, as.list( mu ),
    function( iterTime, outcome, args ) {
      data = solverInputs[[ outcome ]]$data
      args[[ "tau" ]] = solverInputs[[ outcome ]]$tau
      admin = .analyticEvalInfusionAdmin(
        data$duringAndAfter[ iterTime ],
        data$indicesDoses[ iterTime ],
        .analyticInfusionSampleTimes( data, iterTime ),
        solverInputs[[ outcome ]]$dose,
        solverInputs[[ outcome ]]$Tinf,
        args,
        paste0( "t_", outcome ),
        paste0( "dose_", outcome ),
        paste0( "Tinf_", outcome ),
        fnDuringAdmin,
        fnAfterAdmin
      )
      args[[ outcome ]] = admin
      argsNoAdmin = .pfimFillMissingTimeFormals( fnAfterNoAdmin, args, samplings[[ iterTime ]] )
      list( admin = admin, noAdmin = do.call( fnAfterNoAdmin, argsNoAdmin )[[ 1L ]], args = args )
    }
  )
  .analyticFinishEvaluation( tmp, outputNames, arm )
}


#' Dispatch: covariate/occasion structure -> specialised path; else core evaluator.
#' @return Named list of output data frames at requested sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelAnalyticInfusionSteadyState ) = function( model, arm ) {
  .analyticDispatchEvaluate( model, arm, evaluateAnalyticInfusionSteadyStateCore )
}


#' Return PK equations stored on the model (library / user analytic infusion SS).
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List of PK equations from \code{pkModel}.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelAnalyticInfusionSteadyState, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}
