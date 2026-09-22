#' @title ModelAnalyticInfusion
#' @description
#' Closed-form model with infusion (dose in equations).
#' Equations come as \code{duringInfusion} / \code{afterInfusion} pairs; evaluation
#' superposes past doses (after formula) with the active infusion (during formula)
#' using half-open windows \eqn{[t_{\mathrm{dose}},\, t_{\mathrm{dose}}+T_{\mathrm{inf}})}.
#' @inheritParams ModelInfusion
#' @param wrapperModelAnalyticInfusion                   Wrapper for the analytic solver.
#' @param functionArgumentsModelAnalyticInfusion         A list with the function arguments.
#' @param functionArgumentsSymbolModelAnalyticInfusion   A list with the function argument symbols.
#' @param solverInputs                                   A list with the solver inputs.
#' @inheritParams Model
#' @include ModelInfusion.R
#' @include ModelAnalytic.R
#' @return The requested PFIM result.
#' @export

ModelAnalyticInfusion = new_class(
  "ModelAnalyticInfusion",
  package = "PFIM",
  parent  = ModelInfusion,
  properties = list(
    wrapperModelAnalyticInfusion                 = new_property(class_list, default = list()),
    functionArgumentsModelAnalyticInfusion       = new_property(class_list, default = list()),
    functionArgumentsSymbolModelAnalyticInfusion = new_property(class_list, default = list()),
    solverInputs                                 = new_property(class_list, default = list())
  ))


#' Compile during/after infusion analytic wrappers (no steady-state tau formal).
#' @return Updated \code{ModelAnalyticInfusion} object with compiled wrappers.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelAnalyticInfusion ) = function( model, evaluation ) {

  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )
  libraryOutcomesWithAdmin   = .administeredLibraryOutcomeNames( evaluation )

  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNames      = paste0( "t_",    outcomesWithAdministration )
  TinfNames      = paste0( "Tinf_", outcomesWithAdministration )

  # Split equations: administered PK vs passive outcomes (e.g. PD).
  equations               = prop( evaluation, "modelEquations" )
  equationsDuringInfusion = equations$duringInfusion
  equationsAfterInfusion  = equations$afterInfusion

  equationsDuringWithAdmin   = equationsDuringInfusion[  names( equationsDuringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithAdmin    = equationsAfterInfusion[   names( equationsAfterInfusion  ) %in% libraryOutcomesWithAdmin ]
  equationsDuringWithNoAdmin = equationsDuringInfusion[ !names( equationsDuringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithNoAdmin  = equationsAfterInfusion[  !names( equationsAfterInfusion  ) %in% libraryOutcomesWithAdmin ]

  outputDuringAdmin   = names( equationsDuringWithAdmin )
  outputDuringNoAdmin = names( equationsDuringWithNoAdmin )
  outputAfterAdmin    = names( equationsAfterWithAdmin )
  outputAfterNoAdmin  = names( equationsAfterWithNoAdmin )
  timeNamesNoAdmin = unique( c(
    .libraryEquationTimeNames( evaluation, outputDuringNoAdmin ),
    .libraryEquationTimeNames( evaluation, outputAfterNoAdmin )
  ) )

  functionArgumentsAdmin   = unique( c( doseNames, TinfNames, parameterNames, timeNames ) )
  functionArgumentsNoAdmin = unique( c( outcomesWithAdministration, parameterNames, timeNamesNoAdmin ) )

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
  prop( model, "outputNames" ) = unlist( names( equationsDuringInfusion ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


#' Build dose / infusion / sampling tables (optional tau expands a regular grid).
#'
#' Labels each sampling as \code{duringInfusion} or \code{afterInfusion} on
#' half-open windows \eqn{[t_{\mathrm{dose}},\, t_{\mathrm{dose}}+T_{\mathrm{inf}})},
#' and records the active dose index plus relative times since each past dose
#' (for superposition in \code{evaluateAnalyticInfusionCore}).
#' @return Updated \code{ModelAnalyticInfusion} object with solver inputs.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelAnalyticInfusion ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  samplings                  = .analyticSamplingGrid( arm )

  solverInputs = map( administrations, function( adm ) {
    dosing = .alignAdministrationDosing( adm )
    tbl    = .analyticInfusionWindowTable(
      samplings, dosing$timeDose, dosing$dose, dosing$Tinf, prop( adm, "tau" )
    )
    list( data = tbl$data, dose = tbl$dose, Tinf = tbl$Tinf )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}


#' Evaluate analytic infusion predictions (superposition of during/after doses).
#'
#' Unlike the steady-state infusion core, arguments are passed as named lists
#' (vectorised over past doses when several infusions contribute).
#' @param model A \code{ModelAnalyticInfusion} object.
#' @param arm   An \code{Arm} object.
#' @return Named list of output data frames at sampling times.
#' @keywords internal
evaluateAnalyticInfusionCore = function( model, arm ) {

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
  fnDuringAdmin = wrappers[[ 1L ]]
  fnAfterAdmin  = wrappers[[ 3L ]]
  fnAfterNoAdmin = wrappers[[ 4L ]]
  mu = .extractMu( prop( model, "modelParameters" ) )

  tmp = .analyticEvalGrid(
    samplings, outcomesWithAdministration, as.list( mu ),
    function( iterTime, outcome, args ) {
      data = solverInputs[[ outcome ]]$data
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


#' Default method for \code{ModelAnalyticInfusion}.
#' @return Named list of output data frames at sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelAnalyticInfusion ) = function( model, arm ) {
  .analyticDispatchEvaluate( model, arm, evaluateAnalyticInfusionCore )
}


#' Default method for \code{ModelAnalyticInfusion}.
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List of PK equations from \code{pkModel}.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelAnalyticInfusion, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}


#' Default method for \code{ModelAnalyticInfusion}.
#' @param pkModel First argument of generic.
#' @param pdModel Analytic PD model.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List with during/after infusion equation sets.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelAnalyticInfusion, ModelAnalytic, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {
    pkEq = prop( pkModel, "modelEquations" )
    pdEq = prop( pdModel, "modelEquations" )
    list(
      duringInfusion = c( pkEq$duringInfusion, pdEq ),
      afterInfusion  = c( pkEq$afterInfusion,  pdEq )
    )
  }


#' Default method for \code{ModelAnalyticInfusion}.
#' @param pkModel First argument of generic.
#' @param pdModel ODE PD model.
#' @param pfimproject \code{PFIMProject} used for equation remapping.
#' @return List with during/after infusion combined ODE equations.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelAnalyticInfusion, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    pkEqODE = remapPkpdLibraryEquations(
      convertPKModelAnalyticToPKModelODE( pkModel ),
      pfimproject
    )
    pdEq = remapPkpdLibraryEquations(
      prop( pdModel, "modelEquations" ),
      pfimproject
    )
    derivNames = .derivativeNamesFromCompartments(
      pfimproject, 2L,
      c( names( pkEqODE$duringInfusion ), names( pdEq ) )
    )

    list(
      duringInfusion = .combineInfusionPkPdEquations(
        pkEqODE$duringInfusion, pdEq, derivNames
      ),
      afterInfusion = .combineInfusionPkPdEquations(
        pkEqODE$afterInfusion, pdEq, derivNames
      )
    )
  }
