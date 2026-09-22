#' ModelAnalyticSteadyState
#' @title ModelAnalyticSteadyState
#' @description Analytic model at steady state (tau in equations).
#' @inheritParams ModelAnalytic
#' @param wrapperModelAnalytic                   Wrapper for the analytic solver.
#' @param functionArgumentsModelAnalytic         A list with the function arguments.
#' @param functionArgumentsSymbolModelAnalytic   A list with the function argument symbols.
#' @param solverInputs                           A list with the solver inputs.
#' @inheritParams Model
#' @include ModelAnalytic.R
#' @return An S7 object of class \code{ModelAnalyticSteadyState}.
#' @export

ModelAnalyticSteadyState = new_class( "ModelAnalyticSteadyState",
                                      package = "PFIM",
                                      parent  = ModelAnalytic,
                                      properties = list(
                                        wrapperModelAnalytic                 = new_property(class_list, default = list()),
                                        functionArgumentsModelAnalytic       = new_property(class_list, default = list()),
                                        functionArgumentsSymbolModelAnalytic = new_property(class_list, default = list()),
                                        solverInputs                         = new_property(class_list, default = list())
                                      ))


#' Compile bolus / oral analytic wrappers at steady state (\code{tau} in formals).
#' @return Updated model with compiled analytic wrappers.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelAnalyticSteadyState ) = function( model, evaluation ) {

  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )

  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNamesAdmin = paste0( "t_",    outcomesWithAdministration )

  # Administered vs passive equations (PD typically has no dose_*).
  equations            = prop( evaluation, "modelEquations" )
  equationsWithAdmin   = equations[  names( equations ) %in% outcomesWithAdministration ]
  equationsWithNoAdmin = equations[ !names( equations ) %in% outcomesWithAdministration ]

  outputAdmin   = names( equationsWithAdmin )
  outputNoAdmin = names( equationsWithNoAdmin )
  timeNamesNoAdmin = .libraryEquationTimeNames( evaluation, outputNoAdmin )

  # Steady-state closed forms need the dosing interval tau.
  functionArgumentsWithAdmin   = unique( c( doseNames, parameterNames, timeNamesAdmin, "tau" ) )
  functionArgumentsWithNoAdmin = unique( c( outcomesWithAdministration, parameterNames, timeNamesNoAdmin, "tau" ) )

  prop( model, "wrapperModelAnalytic" ) = list(
    functionDefinitionWithAdmin   = .buildAnalyticWrapper( equationsWithAdmin,   functionArgumentsWithAdmin,   timeNamesAdmin,   outputAdmin   ),
    functionDefinitionWithNoAdmin = .buildAnalyticWrapper( equationsWithNoAdmin, functionArgumentsWithNoAdmin, timeNamesNoAdmin, outputNoAdmin )
  )
  prop( model, "functionArgumentsModelAnalytic" ) = list(
    functionArgumentsWithAdmin   = functionArgumentsWithAdmin,
    functionArgumentsWithNoAdmin = functionArgumentsWithNoAdmin
  )
  prop( model, "outputNames" ) = unlist( names( equations ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


#' Expand steady-state doses and relative sampling times for each administration.
#' @return Updated model with solver input structures.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelAnalyticSteadyState ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  samplings                  = .analyticSamplingGrid( arm )

  solverInputs = map( administrations, function( adm ) {
    dosing = .alignAdministrationDosing( adm )
    .analyticRelativeDoseTable(
      samplings, dosing$timeDose, dosing$dose, prop( adm, "tau" )
    )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}


#' Evaluate analytic steady-state predictions (sum over contributing doses).
#'
#' Same superposition pattern as \code{evaluateAnalyticCore}, but each wrapper
#' call also receives \code{tau} (dosing interval) required by the closed forms.
#' @param model A \code{ModelAnalyticSteadyState} object.
#' @param arm   An \code{Arm} object.
#' @return Named list of output data frames at requested sampling times.
#' @keywords internal
evaluateAnalyticSteadyStateCore = function( model, arm ) {

  parameters                 = prop( model, "modelParameters" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  outputNames                = prop( model, "outputNames" )
  samplings                  = prop( model, "samplings" )
  solverInputs               = prop( model, "solverInputs" )

  wrappers = .analyticBindUserEnv(
    prop( model, "wrapperModelAnalytic" )$functionDefinitionWithAdmin,
    prop( model, "wrapperModelAnalytic" )$functionDefinitionWithNoAdmin
  )
  fnAdmin   = wrappers[[ 1L ]]
  fnNoAdmin = wrappers[[ 2L ]]
  mu        = .extractMu( parameters )

  tmp = .analyticEvalGrid(
    samplings, outcomesWithAdministration, as.list( mu ),
    function( iterTime, outcome, args ) {
      data         = solverInputs[[ outcome ]]$data
      dose         = solverInputs[[ outcome ]]$dose
      tau          = solverInputs[[ outcome ]]$tau
      indicesDoses = data$indicesDoses[ iterTime ]
      argsAdmin    = c(
        args,
        set_names(
          list(
            as.numeric( data[ iterTime, seq_len( indicesDoses ) ] ),
            as.numeric( dose[ seq_len( indicesDoses ) ] )
          ),
          c( paste0( "t_", outcome ), paste0( "dose_", outcome ) )
        )
      )
      argsAdmin[[ "tau" ]] = tau
      admin = sum( do.call( fnAdmin, argsAdmin )[[ 1L ]] )
      args[[ outcome ]] = admin
      args[[ "tau" ]]   = tau
      argsNoAdmin = .pfimFillMissingTimeFormals( fnNoAdmin, args, samplings[[ iterTime ]] )
      list( admin = admin, noAdmin = do.call( fnNoAdmin, argsNoAdmin )[[ 1L ]], args = args )
    }
  )
  .analyticFinishEvaluation( tmp, outputNames, arm )
}


#' Dispatch: covariate/occasion structure -> specialised path; else SS core evaluator.
#' @return Named list of output data frames at requested sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelAnalyticSteadyState ) = function( model, arm ) {
  .analyticDispatchEvaluate( model, arm, evaluateAnalyticSteadyStateCore )
}


#' Return PK equations stored on the steady-state analytic model.
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List of PK equations from \code{pkModel}.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelAnalyticSteadyState, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}


#' Concatenate steady-state analytic PK with analytic PD equation lists.
#' @param pkModel First argument of generic.
#' @param pdModel Analytic PD model.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return Concatenated analytic PK/PD equation list.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelAnalyticSteadyState, ModelAnalytic, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {
    c( prop( pkModel, "modelEquations" ), prop( pdModel, "modelEquations" ) )
  }


#' Convert SS analytic PK to ODE, append ODE PD, and remap library compartments.
#' @param pkModel First argument of generic.
#' @param pdModel ODE PD model.
#' @param pfimproject \code{PFIMProject} used for remapping.
#' @return Named list of remapped ODE equations.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelAnalyticSteadyState, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    # Analytic PK -> ODE derivative, then join with PD and remap RespPK/E tokens.
    equations = c(
      convertPKModelAnalyticToPKModelODE( pkModel ),
      prop( pdModel, "modelEquations" )
    )
    eq = remapPkpdLibraryEquations( equations, pfimproject )
    set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ), names( eq ) ) )
  }
