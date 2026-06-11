#' @title ModelODEInfusionDoseInEquation
#' @description ODE infusion model with dose terms in the equations.
#' @inheritParams ModelODEInfusion
#' @param modelODE              An object \code{modelODE}.
#' @param wrapperModelInfusion  Wrapper for solver.
#' @param solverInputs          A list with the solver inputs.
#' @include ModelODEInfusion.R
#' @export

ModelODEInfusionDoseInEquation = new_class( "ModelODEInfusionDoseInEquation",
                                            package = "PFIM",
                                            parent  = ModelODEInfusion,
                                            properties = list(
                                              modelODE             = new_property(class_function, default = NULL),
                                              wrapperModelInfusion = new_property(class_list,     default = list()),
                                              solverInputs         = new_property(class_list,     default = list())
                                            ))



method( defineModelWrapper, ModelODEInfusionDoseInEquation ) = function( model, evaluation ) {

  outcomesWithAdministration = prop( evaluation, "designs" ) |>
    map( \(d) prop( d, "arms" ) ) |>
    list_flatten() |>
    map( \(arm) prop( arm, "administrations" ) ) |>
    list_flatten() |>
    map_chr( \(adm) prop( adm, "outcome" ) ) |>
    unique()

  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration
  prop( model, "wrapperModelInfusion" ) = prop( evaluation, "modelEquations" )

  outputs = prop( evaluation, "outputs" )
  prop( model, "outputNames" ) = names( outputs )
  .setModelOutputFormulas( model, outputs )

  return( model )
}


method( defineModelAdministration, ModelODEInfusionDoseInEquation ) = function( model, arm ) {

  wrapperModelInfusion       = prop( model, "wrapperModelInfusion" )
  wrapperModelDuringInfusion = wrapperModelInfusion$duringInfusion
  wrapperModelAfterInfusion  = wrapperModelInfusion$afterInfusion
  variableDerivativeNames    = names( wrapperModelDuringInfusion )

  parameters               = prop( model, "modelParameters" )
  parameterNames           = map_chr( parameters, "name" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )

  samplingTimes = prop( arm, "samplingTimes" )
  samplings     = map( samplingTimes, ~ prop( .x, "samplings" ) ) |>
    unlist() |> sort() |> unique() |> (\(s) unique( c( 0, s ) ))()

  maxSampling   = max( unlist( map( samplingTimes, ~ prop( .x, "samplings" ) ) ) )
  outputFormula = .getOutputFormulaParsed( model )

  administrations = prop( arm, "administrations" )
  solverInputs    = map( administrations, function( adm ) {
    outcome  = prop( adm, "outcome" )
    tau      = prop( adm, "tau" )
    dosing   = .alignAdministrationDosing( adm )
    timeDose = dosing$timeDose
    dose     = dosing$dose
    Tinf     = dosing$Tinf

    if ( tau != 0 ) {
      timeDose = seq( 0, maxSampling, tau )
      dose     = rep( dose, length( timeDose ) )
      Tinf     = rep( Tinf, length( timeDose ) )
    }
    administrationTime = cbind( timeDose, timeDose + Tinf ) |> unname()
    set_names( list( list( administrationTime = administrationTime, dose = dose, Tinf = Tinf ) ), outcome )
  }) |> list_flatten()

  initialConditions = evaluateInitialConditions( model, arm )

  mu = .extractMu( parameters )
  list2env( as.list( mu ), envir = environment() )

  variableNames = names( initialConditions )
  doseNames     = paste0( "dose_", outcomesWithAdministration )
  tinfNames     = paste0( "Tinf_", outcomesWithAdministration )
  timeNames     = paste0( "t_",    outcomesWithAdministration )

  wrapperModelDuringInfusion = .odeSubstituteBareT(
    wrapperModelDuringInfusion, outcomesWithAdministration
  )
  wrapperModelAfterInfusion = .odeSubstituteBareT(
    wrapperModelAfterInfusion, outcomesWithAdministration
  )

  functionArguments = unique( c( doseNames, tinfNames, timeNames, parameterNames, variableNames ) )
  solverInputs$functionArguments        = functionArguments
  solverInputs$functionArgumentsSymbols = map( functionArguments, as.symbol )

  # Pre-compile wrappers once (not at every ODE step).
  bodyDuring = c(
    map_chr( names( wrapperModelDuringInfusion ), ~ sprintf( "%s = %s", .x, wrapperModelDuringInfusion[[.x]] ) ),
    sprintf( "return(list(c(%s)))", paste( variableDerivativeNames, collapse = ", " ) )
  ) |> paste( collapse = "\n" )
  argsDuring = unique( c( doseNames, tinfNames, timeNames, parameterNames, variableNames ) )
  wrapperDuring = eval( parse( text = sprintf( "function(%s) { %s }",
                                               paste( argsDuring, collapse = ", " ), bodyDuring ) ) )

  bodyAfter = c(
    map_chr( names( wrapperModelAfterInfusion ), ~ sprintf( "%s = %s", .x, wrapperModelAfterInfusion[[.x]] ) ),
    sprintf( "return(list(c(%s)))", paste( variableDerivativeNames, collapse = ", " ) )
  ) |> paste( collapse = "\n" )
  # After-infusion PK equations still reference dose_<outcome> and Tinf_<outcome>.
  argsAfter = argsDuring
  wrapperAfter = eval( parse( text = sprintf( "function(%s) { %s }",
                                              paste( argsAfter, collapse = ", " ), bodyAfter ) ) )

  .infusionCallArgs = function( argNames, state, doseTimeVars ) {
    args = c( as.list( mu ), state, doseTimeVars )
    stats::setNames( lapply( argNames, \( nm) args[[ nm ]] ), argNames )
  }

  modelODEInfusion = function( t, y, parms ) {

    state = set_names( as.list( y ), variableNames )

    doseTimeVars = map( outcomesWithAdministration, \( outcome ) {
      aTime = parms[[ outcome ]]$administrationTime
      idxInf = which( t >= aTime[, 1L] & t < aTime[, 2L] )
      idxDose = which( t >= aTime[, 1L] )
      idx = if ( length( idxInf ) > 0L ) idxInf[1L] else if ( length( idxDose ) > 0L ) idxDose[ length( idxDose ) ] else 1L
      set_names(
        list(
          parms[[ outcome ]]$dose[ idx ],
          parms[[ outcome ]]$Tinf[ idx ],
          t - aTime[ idx, 1L ]
        ),
        c(
          paste0( "dose_", outcome ),
          paste0( "Tinf_", outcome ),
          paste0( "t_", outcome )
        )
      )
    }) |> list_flatten()

    inInfusion = any( map_lgl( outcomesWithAdministration, \( outcome ) {
      aTime = parms[[ outcome ]]$administrationTime
      any( t >= aTime[, 1L] & t < aTime[, 2L] )
    } ) )

    argNames = if ( inInfusion ) argsDuring else argsAfter
    callArgs = .infusionCallArgs( argNames, state, doseTimeVars )

    evaluationModel = do.call(
      if ( inInfusion ) wrapperDuring else wrapperAfter,
      callArgs
    )

    env = list2env( c( as.list( mu ), state, doseTimeVars ), parent = baseenv() )
    evaluationOutputs = map( outputFormula, ~ eval( .x, envir = env ) )
    c( evaluationModel, evaluationOutputs )
  }

  prop( model, "initialConditions" ) = initialConditions
  prop( model, "samplings" ) = samplings
  prop( model, "modelODE" ) = modelODEInfusion
  prop( model, "solverInputs" ) = solverInputs
  model
}


method( evaluateModel, ModelODEInfusionDoseInEquation ) = function( model, arm ) {

  odeSolverParameters = prop( model, "odeSolverParameters" )
  outputNames         = prop( model, "outputNames" )
  samplingTimes       = prop( arm,   "samplingTimes" )

  evaluationModelTmp = ode(
    prop( model, "initialConditions" ),
    prop( model, "samplings" ),
    prop( model, "modelODE" ),
    prop( model, "solverInputs" ),
    hmax = 0.0,
    atol = odeSolverParameters$atol,
    rtol = odeSolverParameters$rtol
  ) |> as.data.frame()

  .odeExtractOutputAtSamplingTimes(
    evaluationModelTmp, samplingTimes, outputNames, .getModelOutputFormulas( model )
  )
}


method( definePKModel, list( ModelODEInfusionDoseInEquation, PFIMProject ) ) = function( pkModel, pfimproject ) {

  pkModelEquations = prop( pkModel, "modelEquations" )
  nPk = length( pkModelEquations$duringInfusion )
  derivNames = .derivativeNamesFromCompartments( pfimproject, nPk )

  pkModelEquations$duringInfusion = remapOdePkLibraryEquations(
    pkModelEquations$duringInfusion, pfimproject
  ) |> set_names( derivNames )

  pkModelEquations$afterInfusion = remapOdePkLibraryEquations(
    pkModelEquations$afterInfusion, pfimproject
  ) |> set_names( derivNames )

  pkModelEquations
}


method( definePKPDModel, list( ModelODEInfusionDoseInEquation, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    pkModelEquations = remapOdePkLibraryEquations(
      prop( pkModel, "modelEquations" ),
      pfimproject
    )
    pdEq = remapPkpdLibraryEquations(
      prop( pdModel, "modelEquations" ),
      pfimproject
    )
    derivNames = .derivativeNamesFromCompartments( pfimproject, 2L )

    list(
      duringInfusion = .combineInfusionPkPdEquations(
        pkModelEquations$duringInfusion, pdEq, derivNames
      ),
      afterInfusion = .combineInfusionPkPdEquations(
        pkModelEquations$afterInfusion, pdEq, derivNames
      )
    )
  }
