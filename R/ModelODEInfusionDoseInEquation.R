#' @title ModelODEInfusionDoseInEquation
#' @description ODE infusion model with dose terms in the equations.
#' @inheritParams ModelODEInfusion
#' @param modelODE              An object \code{modelODE}.
#' @param wrapperModelInfusion  Wrapper for solver.
#' @param solverInputs          A list with the solver inputs.
#' @include ModelODEInfusion.R
#' @return An S7 object of class \code{ModelODEInfusionDoseInEquation}.
#' @export

ModelODEInfusionDoseInEquation = new_class( "ModelODEInfusionDoseInEquation",
                                            package = "PFIM",
                                            parent  = ModelODEInfusion,
                                            properties = list(
                                              modelODE             = new_property(class_function, default = NULL),
                                              wrapperModelInfusion = new_property(class_list,     default = list()),
                                              solverInputs         = new_property(class_list,     default = list())
                                            ))



#' Store during/after infusion equations and output mapping on the model.
#'
#' Does not compile wrappers yet; compilation happens in
#' \code{defineModelAdministration} once arm dosing is known.
#' @return Updated model with wrapper metadata.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelODEInfusionDoseInEquation ) = function( model, evaluation ) {

  # Outcomes linked to administrations across all designs/arms.
  outcomesWithAdministration = prop( evaluation, "designs" ) |>
    map( \(d) prop( d, "arms" ) ) |>
    list_flatten() |>
    map( \(arm) prop( arm, "administrations" ) ) |>
    list_flatten() |>
    map_chr( \(adm) prop( adm, "outcome" ) ) |>
    unique()

  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration
  eqs = prop( evaluation, "modelEquations" )
  during = eqs$duringInfusion
  if ( is.null( during ) ) during = eqs
  prop( model, "variableNames" ) = .pfimOdeCompartmentNames( during )
  prop( model, "wrapperModelInfusion" ) = eqs

  outputs = prop( evaluation, "outputs" )
  prop( model, "outputNames" ) = names( outputs )
  .setModelOutputFormulas( model, outputs )

  return( model )
}


#' Build infusion ODE RHS closure from precompiled wrappers and current mu.
#'
#' At each \code{deSolve} step, selects during- vs after-infusion wrapper based
#' on whether \code{t} falls inside any administration window
#' (\eqn{[t_{\mathrm{dose}},\, t_{\mathrm{dose}}+T_{\mathrm{inf}})} - same
#' half-open convention as analytic infusion), and exposes \code{dose_*},
#' \code{Tinf_*}, and relative \code{t_*} for each outcome.
#' @return Function suitable as \code{deSolve::ode} RHS.
#' @noRd
#' @keywords internal
.odeInfusionModelOde = function(
    mu, wrapperDuring, wrapperAfter, argsDuring, argsAfter,
    variableNames, variableDerivativeNames, outcomesWithAdministration,
    outputFormula ) {
  # Align named args with wrapper formals from mu + state + dose/time vars.
  .infusionCallArgs = function( argNames, state, doseTimeVars ) {
    args = c( as.list( mu ), state, doseTimeVars )
    stats::setNames( lapply( argNames, \( nm) args[[ nm ]] ), argNames )
  }

  function( t, y, parms ) {
    state = set_names( as.list( y ), variableNames )

    # Active dose index: prefer current infusion window, else last started dose.
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

    # Any administered outcome currently inside [t_dose, t_dose + Tinf).
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

    # Evaluate output formulas in an env with mu, states, and dose/time vars.
    env = list2env( c( as.list( mu ), state, doseTimeVars ), parent = baseenv() )
    evaluationOutputs = map( outputFormula, ~ eval( .x, envir = env ) )
    c( evaluationModel, evaluationOutputs )
  }
}

#' Build infusion dose windows passed to \code{deSolve::ode}.
#'
#' For \code{tau != 0}, expands doses onto a regular grid up to the last
#' sampling time. Each outcome stores \code{administrationTime} as
#' \code{[t_dose, t_dose + Tinf)} plus dose/Tinf vectors.
#' @param model \code{ModelODEInfusionDoseInEquation} object.
#' @param arm \code{Arm} object with administrations and sampling times.
#' @return Named list of per-outcome dosing windows plus function argument metadata.
#' @noRd
#' @keywords internal
.odeInfusionSolverInputsFromArm = function( model, arm ) {
  samplingTimes = prop( arm, "samplingTimes" )
  maxSampling   = max( unlist( map( samplingTimes, ~ prop( .x, "samplings" ) ) ) )
  parameterNames = map_chr( prop( model, "modelParameters" ), "name" )
  initialConditions = evaluateInitialConditions( model, arm )
  variableNames = .pfimOdeStateNames( model )
  if ( !length( variableNames ) )
    variableNames = names( initialConditions )
  initialConditions = .pfimAlignOdeStates( initialConditions, variableNames )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  doseNames = paste0( "dose_", outcomesWithAdministration )
  tinfNames = paste0( "Tinf_", outcomesWithAdministration )
  timeNames = paste0( "t_", outcomesWithAdministration )

  solverInputs = map( prop( arm, "administrations" ), function( adm ) {
    outcome  = prop( adm, "outcome" )
    tau      = prop( adm, "tau" )
    dosing   = .alignAdministrationDosing( adm )
    timeDose = dosing$timeDose
    dose     = dosing$dose
    Tinf     = dosing$Tinf

    # Steady-state / multi-dose: replicate the unit infusion every tau.
    if ( tau != 0 ) {
      timeDose = seq( 0, maxSampling, tau )
      dose     = rep( dose, length( timeDose ) )
      Tinf     = rep( Tinf, length( timeDose ) )
    }
    administrationTime = cbind( timeDose, timeDose + Tinf ) |> unname()
    set_names( list( list( administrationTime = administrationTime, dose = dose, Tinf = Tinf ) ), outcome )
  }) |> list_flatten()

  functionArguments = unique( c( doseNames, tinfNames, timeNames, parameterNames, variableNames ) )
  solverInputs$functionArguments        = functionArguments
  solverInputs$functionArgumentsSymbols = map( functionArguments, as.symbol )
  solverInputs
}

#' Compile during/after infusion ODE wrappers and administration inputs for one arm.
#'
#' Substitutes bare \code{t} with \code{t_<outcome>}, builds executable during/after
#' wrappers, and constructs the deSolve RHS closure with current mu.
#' @param model \code{ModelODEInfusionDoseInEquation} object.
#' @param arm \code{Arm} object.
#' @return Named list of administration parts (IC, samplings, wrappers, RHS, ...).
#' @noRd
#' @keywords internal
.odeInfusionAdminParts = function( model, arm ) {
  wrapperModelInfusion       = prop( model, "wrapperModelInfusion" )
  wrapperModelDuringInfusion = wrapperModelInfusion$duringInfusion
  wrapperModelAfterInfusion  = wrapperModelInfusion$afterInfusion
  variableDerivativeNames    = names( wrapperModelDuringInfusion )

  parameters                 = prop( model, "modelParameters" )
  parameterNames             = map_chr( parameters, "name" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )

  samplingTimes = prop( arm, "samplingTimes" )
  samplings     = map( samplingTimes, ~ prop( .x, "samplings" ) ) |>
    unlist() |> sort() |> unique() |> (\(s) unique( c( 0, s ) ))()

  outputFormula = .getOutputFormulaParsed( model )
  solverInputs  = .odeInfusionSolverInputsFromArm( model, arm )
  initialConditions = evaluateInitialConditions( model, arm )
  mu = .extractMu( parameters )

  variableNames = .pfimOdeStateNames( model )
  if ( !length( variableNames ) )
    variableNames = str_remove( variableDerivativeNames, "Deriv_" )
  initialConditions = .pfimAlignOdeStates( initialConditions, variableNames )
  prop( model, "variableNames" ) = variableNames
  doseNames     = paste0( "dose_", outcomesWithAdministration )
  tinfNames     = paste0( "Tinf_", outcomesWithAdministration )
  timeNames     = paste0( "t_",    outcomesWithAdministration )

  wrapperModelDuringInfusion = .odeSubstituteBareT(
    wrapperModelDuringInfusion, outcomesWithAdministration
  )
  wrapperModelAfterInfusion = .odeSubstituteBareT(
    wrapperModelAfterInfusion, outcomesWithAdministration
  )

  argsDuring = unique( c( doseNames, tinfNames, timeNames, parameterNames, variableNames ) )
  bodyDuring = c(
    map_chr( names( wrapperModelDuringInfusion ), ~ sprintf( "%s = %s", .x, wrapperModelDuringInfusion[[.x]] ) ),
    sprintf( "return(list(c(%s)))", paste( variableDerivativeNames, collapse = ", " ) )
  ) |> paste( collapse = "\n" )
  wrapperDuring = eval( parse( text = sprintf( "function(%s) { %s }",
                                               paste( argsDuring, collapse = ", " ), bodyDuring ) ) )

  bodyAfter = c(
    map_chr( names( wrapperModelAfterInfusion ), ~ sprintf( "%s = %s", .x, wrapperModelAfterInfusion[[.x]] ) ),
    sprintf( "return(list(c(%s)))", paste( variableDerivativeNames, collapse = ", " ) )
  ) |> paste( collapse = "\n" )
  argsAfter = argsDuring
  wrapperAfter = eval( parse( text = sprintf( "function(%s) { %s }",
                                              paste( argsAfter, collapse = ", " ), bodyAfter ) ) )

  modelODE = .odeInfusionModelOde(
    mu, wrapperDuring, wrapperAfter, argsDuring, argsAfter,
    variableNames, variableDerivativeNames, outcomesWithAdministration,
    outputFormula
  )

  list(
    initialConditions = initialConditions,
    samplings = samplings,
    solverInputs = solverInputs,
    modelODE = modelODE,
    wrapperDuring = wrapperDuring,
    wrapperAfter = wrapperAfter,
    argsDuring = argsDuring,
    argsAfter = argsAfter,
    variableNames = variableNames,
    variableDerivativeNames = variableDerivativeNames,
    outcomesWithAdministration = outcomesWithAdministration,
    outputFormula = outputFormula
  )
}

#' Apply compiled infusion administration slots to a model.
#' @noRd
#' @keywords internal
.odeInfusionApplyParts = function( model, parts ) {
  prop( model, "initialConditions" ) = parts$initialConditions
  prop( model, "samplings" ) = parts$samplings
  prop( model, "modelODE" ) = parts$modelODE
  prop( model, "solverInputs" ) = parts$solverInputs
  model
}

#' Cacheable infusion administration entry (compiled wrappers + arm layout).
#' @noRd
#' @keywords internal
.pfimInfusionAdminEntryFromParts = function( parts ) {
  list(
    type = "infusionDoseInEq",
    samplings = parts$samplings,
    wrapperDuring = parts$wrapperDuring,
    wrapperAfter = parts$wrapperAfter,
    argsDuring = parts$argsDuring,
    argsAfter = parts$argsAfter,
    variableDerivativeNames = parts$variableDerivativeNames,
    outcomesWithAdministration = parts$outcomesWithAdministration
  )
}

#' Apply cached infusion administration parts with refreshed mu.
#' @noRd
#' @keywords internal
.odeInfusionApplyAdminEntry = function( model, arm, entry ) {
  mu = .extractMu( prop( model, "modelParameters" ) )
  variableNames = .pfimOdeStateNames( model )
  if ( !length( variableNames ) )
    variableNames = str_remove( entry$variableDerivativeNames, "Deriv_" )
  initialConditions = .pfimAlignOdeStates(
    evaluateInitialConditions( model, arm ), variableNames
  )
  prop( model, "variableNames" ) = variableNames
  prop( model, "initialConditions" ) = initialConditions
  prop( model, "samplings" ) = entry$samplings
  prop( model, "solverInputs" ) = .odeInfusionSolverInputsFromArm( model, arm )
  prop( model, "modelODE" ) = .odeInfusionModelOde(
    mu,
    entry$wrapperDuring, entry$wrapperAfter,
    entry$argsDuring, entry$argsAfter,
    variableNames, entry$variableDerivativeNames,
    entry$outcomesWithAdministration,
    .getOutputFormulaParsed( model )
  )
  model
}


#' Compile infusion wrappers and bind arm dosing to the model.
#' @return Updated model with ODE solver function and inputs.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelODEInfusionDoseInEquation ) = function( model, arm ) {
  .odeInfusionApplyParts( model, .odeInfusionAdminParts( model, arm ) )
}


#' Integrate ODE infusion trajectories and extract outputs at sampling times.
#' @param model \code{ModelODEInfusionDoseInEquation} object.
#' @param arm \code{Arm} object used for evaluation.
#' @return Named list of output data frames at requested sampling times.
#' @noRd
#' @keywords internal
.odeInfusionEvaluateModelCore = function( model, arm ) {

  odeSolverParameters = prop( model, "odeSolverParameters" )
  outputNames         = prop( model, "outputNames" )
  samplingTimes       = prop( arm,   "samplingTimes" )
  tol                 = .pfimDeSolveTolerances( odeSolverParameters )

  # RHS switches during/after infusion using solverInputs administration windows.
  evaluationModelTmp = ode(
    prop( model, "initialConditions" ),
    prop( model, "samplings" ),
    prop( model, "modelODE" ),
    prop( model, "solverInputs" ),
    hmax = 0.0,
    atol = tol$atol,
    rtol = tol$rtol
  ) |> as.data.frame()

  .odeExtractOutputAtSamplingTimes(
    evaluationModelTmp, samplingTimes, outputNames, .getModelOutputFormulas( model )
  )
}


#' Dispatch: covariate/occasion structure -> specialised path; else core evaluator.
#' @return Named list of output data frames at requested sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelODEInfusionDoseInEquation ) = function( model, arm ) {
  if ( usesCovariateOccasionStructure( model ) )
    evaluateModelWithCovariates( model, arm, .odeInfusionEvaluateModelCore )
  else
    .odeInfusionEvaluateModelCore( model, arm )
}


#' Remap library ODE infusion PK equations onto project compartments.
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} used for compartment remapping.
#' @return List of remapped PK equations for infusion ODE models.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelODEInfusionDoseInEquation, PFIMProject ) ) = function( pkModel, pfimproject ) {

  pkModelEquations = prop( pkModel, "modelEquations" )
  nPk = length( pkModelEquations$duringInfusion )
  derivNames = .derivativeNamesFromCompartments(
    pfimproject, nPk, names( pkModelEquations$duringInfusion )
  )

  # Remap C1/C2 tokens and rename to Deriv_<compartment> for during and after.
  pkModelEquations$duringInfusion = remapOdePkLibraryEquations(
    pkModelEquations$duringInfusion, pfimproject
  ) |> set_names( derivNames )

  pkModelEquations$afterInfusion = remapOdePkLibraryEquations(
    pkModelEquations$afterInfusion, pfimproject
  ) |> set_names( derivNames )

  pkModelEquations
}


#' Combine remapped infusion PK with ODE PD for during and after phases.
#' @param pkModel First argument of generic.
#' @param pdModel ODE PD model.
#' @param pfimproject \code{PFIMProject} used for equation remapping.
#' @return List with combined PK/PD equations for during and after infusion.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelODEInfusionDoseInEquation, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    pkOrig = prop( pkModel, "modelEquations" )
    pdOrig = prop( pdModel, "modelEquations" )
    catalog = c( names( pkOrig$duringInfusion ), names( pdOrig ) )
    if ( is.null( catalog ) || !length( catalog ) )
      catalog = names( pkOrig )
    pkModelEquations = remapOdePkLibraryEquations(
      pkOrig,
      pfimproject
    )
    pdEq = remapPkpdLibraryEquations(
      pdOrig,
      pfimproject
    )
    derivNames = .derivativeNamesFromCompartments( pfimproject, 2L, catalog )

    list(
      duringInfusion = .combineInfusionPkPdEquations(
        pkModelEquations$duringInfusion, pdEq, derivNames
      ),
      afterInfusion = .combineInfusionPkPdEquations(
        pkModelEquations$afterInfusion, pdEq, derivNames
      )
    )
  }
