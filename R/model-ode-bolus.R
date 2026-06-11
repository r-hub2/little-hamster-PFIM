#' @title ModelODEBolus
#' @description ODE bolus model with doses in initial conditions.
#' @inheritParams ModelODE
#' @param modelODE \code{deSolve} right-hand side built at administration time.
#' @param doseEvent Bolus dose event table passed to \code{deSolve}.
#' @param solverInputs Reserved; not used for bolus initial-condition dosing.
#' @include ModelODE.R
#' @export
ModelODEBolus = new_class(
  "ModelODEBolus", package = "PFIM", parent = ModelODE,
  properties = list(
    modelODE     = new_property( class_function, default = NULL ),
    doseEvent    = new_property( class_list,     default = list() ),
    solverInputs = new_property( class_list,     default = list() )
  )
)

#' @title ModelODEDoseInEquations
#' @description ODE model with bolus dose terms inside the differential equations.
#' @inheritParams ModelODE
#' @param modelODEDoseInEquations \code{deSolve} right-hand side with \code{dose_} terms.
#' @param solverInputs Per-outcome administration windows and dose levels.
#' @include Model.R
#' @export
ModelODEDoseInEquations = new_class(
  "ModelODEDoseInEquations", package = "PFIM", parent = ModelODE,
  properties = list(
    modelODEDoseInEquations = new_property( class_function, default = NULL ),
    solverInputs            = new_property( class_list,     default = list() )
  )
)

#' @title ModelODEDoseNotInEquations
#' @description ODE model with bolus doses added as compartment events.
#' @inheritParams ModelODE
#' @param modelODE \code{deSolve} right-hand side built at administration time.
#' @param doseEvent Bolus dose event table passed to \code{deSolve}.
#' @param solverInputs Reserved; not used for compartment-event bolus dosing.
#' @include Model.R
#' @export
ModelODEDoseNotInEquations = new_class(
  "ModelODEDoseNotInEquations", package = "PFIM", parent = ModelODE,
  properties = list(
    modelODE     = new_property( class_function, default = NULL ),
    doseEvent    = new_property( class_list,     default = list() ),
    solverInputs = new_property( class_list,     default = list() )
  )
)

# Bolus ODE core: three administration modes share one solver path (deSolve).
#   bolusIc    — dose added to compartment state at t = 0 (initial conditions)
#   doseInEq   — dose_* terms inside the RHS, per-outcome administration windows
#   doseEvent  — compartment events via doseEvent table

.odeBolusMode = function( model ) {
  if ( S7::S7_inherits( model, ModelODEBolus ) ) return( "bolusIc" )
  if ( S7::S7_inherits( model, ModelODEDoseInEquations ) ) return( "doseInEq" )
  if ( S7::S7_inherits( model, ModelODEDoseNotInEquations ) ) return( "doseEvent" )
  stop( "not a bolus ODE model", call. = FALSE )
}

.odeSubstituteBareT = function( equations, outcomes ) {
  map( equations, function( eq ) {
    reduce( outcomes, function( text, outcome ) {
      if ( str_detect( text, paste0( "dose_", outcome ) ) )
        str_replace_all( text, "\\bt\\b", paste0( "t_", outcome ) )
      else
        text
    }, .init = eq )
  } )
}

.odeDefineWrapper = function( model, evaluation, mode ) {
  equations                  = prop( evaluation, "modelEquations" )
  variableNames              = str_remove( names( equations ), "Deriv_" )
  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )
  parameterNames             = map_chr( prop( evaluation, "modelParameters" ), "name" )

  if ( mode == "doseInEq" ) {
    doseNames = paste0( "dose_", outcomesWithAdministration )
    timeNames = paste0( "t_", outcomesWithAdministration )
    functionArguments = unique( c( doseNames, parameterNames, variableNames, timeNames ) )
    equations = .odeSubstituteBareT( equations, outcomesWithAdministration )
  } else {
    functionArguments = unique( c( parameterNames, variableNames, "t" ) )
  }

  outputs = prop( evaluation, "outputs" )
  prop( model, "outputNames" ) = names( outputs )
  if ( mode %in% c( "bolusIc", "doseInEq" ) )
    prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration
  prop( model, "wrapper" ) = .buildODEWrapper( equations, functionArguments )
  prop( model, "functionArguments" ) = functionArguments
  prop( model, "functionArgumentsSymbol" ) = map( functionArguments, as.symbol )
  .setModelOutputFormulas( model, outputs )
  model
}

.parseInitialConditionExprs = function( initialConditionsTmp ) {
  lapply( initialConditionsTmp, function( ic ) {
    if ( is.character( ic ) ) parse( text = ic )[[ 1L ]] else NULL
  } )
}

.applyInitialConditionsToEvent = function( doseEvent, parsedExprs, env ) {
  for ( iter in seq_len( nrow( doseEvent ) ) ) {
    expr = parsedExprs[[ doseEvent$var[ iter ] ]]
    if ( is.null( expr ) ) next
    assign( paste0( "dose_", doseEvent$var[ iter ] ), doseEvent$value[ iter ], envir = env )
    doseEvent$value[ iter ] = eval( expr, envir = env )
  }
  doseEvent
}

.odeDoseEventFromArm = function( arm, samplings ) {
  map( prop( arm, "administrations" ), function( adm ) {
    outcome  = prop( adm, "outcome" )
    tau      = prop( adm, "tau" )
    dosing   = .alignAdministrationDosing( adm )
    timeDose = dosing$timeDose
    dose     = dosing$dose
    if ( tau != 0 ) {
      timeDose = seq( 0, max( samplings ), tau )
      dose     = rep( dose, length( timeDose ) )
    }
    data.frame(
      var    = rep( outcome, length( timeDose ) ),
      time   = timeDose,
      value  = dose,
      method = ifelse( timeDose == 0, "replace", "add" )
    )
  } ) |> list_rbind() |> ( function( df ) df[ order( df$time ), ] )()
}

.odeSolverInputsFromArm = function( arm, samplings ) {
  map( prop( arm, "administrations" ), function( adm ) {
    outcome  = prop( adm, "outcome" )
    tau      = prop( adm, "tau" )
    dosing   = .alignAdministrationDosing( adm )
    timeDose = dosing$timeDose
    dose     = dosing$dose
    administrationTime = if ( tau != 0 ) {
      tSeq = seq( 0, max( samplings ), tau )
      dose = rep( dose, length( tSeq ) )
      cbind( tSeq[ -length( tSeq ) ], tSeq[ -1 ] )
    } else if ( length( timeDose ) == 1L ) {
      matrix( rep( timeDose, 2 ), nrow = 1 )
    } else {
      cbind( timeDose, c( timeDose[ -1 ], max( samplings ) ) )
    }
    set_names( list( list( administrationTime = administrationTime, dose = dose ) ), outcome )
  } ) |> list_flatten()
}

# Standard RHS: wrapper returns d(state)/dt; outputFormula evaluated in the same env.
.odeRhsStandard = function( wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu ) {
  function( t, y, parms ) {
    with( as.list( c( t = t, y, mu ) ), {
      evaluationModel   = do.call( wrapper, set_names( functionArgumentsSymbols, functionArguments ) )
      evaluationOutputs = map( outputFormula, ~ eval( .x ) )
      c( evaluationModel, evaluationOutputs )
    } )
  }
}

# doseInEq RHS: for each outcome, locate the active infusion window and expose
# dose_<outcome> and t_<outcome> (time since window start) to the wrapper.
.odeRhsDoseInEq = function(
    wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu,
    outcomesWithAdministration ) {
  function( t, y, parms ) {
    doseTimeVars = list()
    for ( outcome in outcomesWithAdministration ) {
      adminTime = parms[[ outcome ]]$administrationTime
      idx       = which( t > adminTime[ , 1 ] & t <= adminTime[ , 2 ] )
      doseTimeVars[[ paste0( "dose_", outcome ) ]] =
        if ( length( idx ) > 0 ) parms[[ outcome ]]$dose[ idx ] else parms[[ outcome ]]$dose[ 1 ]
      doseTimeVars[[ paste0( "t_", outcome ) ]] =
        if ( length( idx ) > 0 ) t - adminTime[ idx, 1 ] else t
    }
    with( as.list( c( list( t = t ), as.list( y ), as.list( mu ), doseTimeVars ) ), {
      evaluationModel   = do.call( wrapper, set_names( functionArgumentsSymbols, functionArguments ) )
      evaluationOutputs = map( outputFormula, ~ eval( .x ) )
      c( evaluationModel, evaluationOutputs )
    } )
  }
}

.odeFinalizeAdministration = function(
    model, mode, mu, initialConditions, samplings,
    wrapper, functionArguments, functionArgumentsSymbols, outputFormula,
    doseEvent = NULL, solverInputs = NULL, outcomesWithAdministration = NULL ) {

  rhs = if ( mode == "doseInEq" ) {
    .odeRhsDoseInEq(
      wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu,
      outcomesWithAdministration
    )
  } else {
    .odeRhsStandard( wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu )
  }

  prop( model, "initialConditions" ) = initialConditions
  prop( model, "samplings" ) = samplings
  if ( mode == "doseInEq" ) {
    prop( model, "modelODEDoseInEquations" ) = rhs
    prop( model, "solverInputs" ) = solverInputs
  } else {
    prop( model, "modelODE" ) = rhs
    prop( model, "doseEvent" ) = doseEvent
  }
  model
}

# bolusIc: doses at t = 0 enter the IC env as dose_<compartment>; IC exprs may
# reference them before deSolve integrates (remaining doses use doseEvent).
.odeInitialConditionsBolus = function( model, arm, doseEvent ) {
  outcomes = prop( model, "outcomesWithAdministration" )
  dosesAtT0 = doseEvent[ doseEvent$time == 0, , drop = FALSE ]
  doses = set_names(
    as.list( dosesAtT0$value[ match( outcomes, dosesAtT0$var ) ] ),
    paste0( "dose_", outcomes )
  )
  mu  = .extractMu( prop( model, "modelParameters" ) )
  env = list2env( c( as.list( mu ), doses ), parent = environment() )
  map( prop( arm, "initialConditions" ), ~ {
    if ( is.numeric( .x ) ) .x else eval( parse( text = .x ), envir = env )
  } ) |> unlist()
}

# doseEvent: zero compartment states at t = 0; all mass via deSolve events table.
.odeInitialConditionsDoseEvent = function( model, arm, doseEvent ) {
  ic = evaluateInitialConditions( model, arm )
  admin = set_names(
    rep( 0, length( unique( doseEvent$var[ doseEvent$time == 0 ] ) ) ),
    unique( doseEvent$var[ doseEvent$time == 0 ] )
  )
  ic = c( admin, ic )
  ic[ !duplicated( names( ic ) ) ]
}

# Bind arm dosing to model slots; branches on bolusIc / doseInEq / doseEvent.
.odeDefineAdministration = function( model, arm, mode ) {
  samplingTimes            = prop( arm, "samplingTimes" )
  samplings                = .buildSamplings( samplingTimes )
  wrapper                  = prop( model, "wrapper" )
  mu                       = .extractMu( prop( model, "modelParameters" ) )
  functionArguments        = prop( model, "functionArguments" )
  functionArgumentsSymbols = prop( model, "functionArgumentsSymbol" )
  outputFormula            = .getOutputFormulaParsed( model )

  if ( mode == "doseInEq" ) {
    .odeFinalizeAdministration(
    model, mode, mu,
    evaluateInitialConditions( model, arm ), samplings,
    wrapper, functionArguments, functionArgumentsSymbols, outputFormula,
    solverInputs = .odeSolverInputsFromArm( arm, samplings ),
    outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
    )
  } else {
    doseEvent = .odeDoseEventFromArm( arm, samplings )
    initialConditions = if ( mode == "bolusIc" ) {
      ic = .odeInitialConditionsBolus( model, arm, doseEvent )
      env = list2env( as.list( mu ), parent = environment() )
      doseEvent = .applyInitialConditionsToEvent(
        doseEvent, .parseInitialConditionExprs( prop( arm, "initialConditions" ) ), env
      )
      ic
    } else {
      .odeInitialConditionsDoseEvent( model, arm, doseEvent )
    }
    .odeFinalizeAdministration(
      model, mode, mu, initialConditions, samplings,
      wrapper, functionArguments, functionArgumentsSymbols, outputFormula,
      doseEvent = doseEvent
    )
  }
}

.odeExtractOutputAtSamplingTimes = function(
    evaluationModelTmp, samplingTimes, outputNames, outputFormula = list() ) {
  samplings_list = set_names(
    map( samplingTimes, ~ prop( .x, "samplings" ) ),
    outputNames
  )
  set_names(
    map( outputNames, function( out_name ) {
      col_name  = .odeColumnForOutput( out_name, evaluationModelTmp, outputFormula )
      req_times = samplings_list[[ out_name ]]
      idx = vapply( req_times, function( rt ) {
        which.min( abs( evaluationModelTmp$time - rt ) )
      }, integer( 1L ), USE.NAMES = FALSE )
      df = data.frame( time = req_times, val = evaluationModelTmp[ idx, col_name, drop = FALSE ] )
      names( df )[ 2L ] = out_name
      df
    } ),
    outputNames
  )
}

# deSolve call shared by the three bolus modes; events only for doseEvent/bolusIc.
.odeSimulateBolus = function( model, arm, mode ) {
  odeSolverParameters = prop( model, "odeSolverParameters" )
  raw_samplings       = prop( model, "samplings" )
  events_df           = if ( mode != "doseInEq" ) prop( model, "doseEvent" ) else NULL
  sim_key             = if ( mode == "doseInEq" ) {
    paste( "doseInEq", paste( raw_samplings, collapse = "," ), sep = "::" )
  } else {
    paste(
      mode, paste( raw_samplings, collapse = "," ),
      paste( events_df$time, events_df$value, sep = ":", collapse = "," ),
      sep = "::"
    )
  }
  sim_times = .pfimOdeSimTimesCached( sim_key, raw_samplings, events_df )
  ode_fn    = if ( mode == "doseInEq" )
    prop( model, "modelODEDoseInEquations" ) else prop( model, "modelODE" )
  parms = if ( mode == "doseInEq" ) prop( model, "solverInputs" ) else NULL

  ode(
    prop( model, "initialConditions" ),
    sim_times,
    ode_fn,
    parms,
    events = if ( !is.null( events_df ) ) list( data = events_df ) else NULL,
    atol   = odeSolverParameters$atol,
    rtol   = odeSolverParameters$rtol
  ) |> as.data.frame()
}

.odeEvaluateBolus = function( model, arm ) {
  mode = .odeBolusMode( model )
  core = function( model, arm ) {
    evaluationModelTmp = .odeSimulateBolus( model, arm, mode )
    .odeExtractOutputAtSamplingTimes(
      evaluationModelTmp,
      prop( arm, "samplingTimes" ),
      prop( model, "outputNames" ),
      .getModelOutputFormulas( model )
    )
  }
  if ( usesCovariateOccasionStructure( model ) )
    evaluateModelWithCovariates( model, arm, core )
  else
    core( model, arm )
}

# Fast path for gradient cache: reapply stored administration template to new mu.
.odeApplyAdminEntry = function( model, arm, entry ) {
  mu = .extractMu( prop( model, "modelParameters" ) )
  mode = entry$type
  if ( mode == "doseInEq" ) {
    .odeFinalizeAdministration(
      model, mode, mu,
      evaluateInitialConditions( model, arm ), entry$samplings,
      entry$wrapper, entry$functionArguments, entry$functionArgumentsSymbols, entry$outputFormula,
      solverInputs = entry$solverInputs,
      outcomesWithAdministration = entry$outcomesWithAdministration
    )
  } else {
    doseEvent = entry$doseEventTemplate
    initialConditions = if ( mode == "bolusIc" ) {
      ic = .odeInitialConditionsBolus( model, arm, doseEvent )
      env = list2env( as.list( mu ), parent = environment() )
      doseEvent = .applyInitialConditionsToEvent(
        doseEvent, entry$initialConditionsParsed, env
      )
      ic
    } else {
      .odeInitialConditionsDoseEvent( model, arm, doseEvent )
    }
    .odeFinalizeAdministration(
      model, mode, mu, initialConditions, entry$samplings,
      entry$wrapper, entry$functionArguments, entry$functionArgumentsSymbols, entry$outputFormula,
      doseEvent = doseEvent
    )
  }
}

method( defineModelWrapper, ModelODEBolus ) = function( model, evaluation ) {
  .odeDefineWrapper( model, evaluation, "bolusIc" )
}
method( defineModelWrapper, ModelODEDoseInEquations ) = function( model, evaluation ) {
  .odeDefineWrapper( model, evaluation, "doseInEq" )
}
method( defineModelWrapper, ModelODEDoseNotInEquations ) = function( model, evaluation ) {
  .odeDefineWrapper( model, evaluation, "doseEvent" )
}

method( evaluateInitialConditions, ModelODEBolus ) = function( model, arm, doseEvent ) {
  .odeInitialConditionsBolus( model, arm, doseEvent )
}

method( defineModelAdministration, ModelODEBolus ) = function( model, arm ) {
  .odeDefineAdministration( model, arm, "bolusIc" )
}
method( defineModelAdministration, ModelODEDoseInEquations ) = function( model, arm ) {
  .odeDefineAdministration( model, arm, "doseInEq" )
}
method( defineModelAdministration, ModelODEDoseNotInEquations ) = function( model, arm ) {
  .odeDefineAdministration( model, arm, "doseEvent" )
}

method( evaluateModel, ModelODEBolus ) = function( model, arm ) {
  .odeEvaluateBolus( model, arm )
}
method( evaluateModel, ModelODEDoseInEquations ) = function( model, arm ) {
  .odeEvaluateBolus( model, arm )
}
method( evaluateModel, ModelODEDoseNotInEquations ) = function( model, arm ) {
  .odeEvaluateBolus( model, arm )
}

method( definePKModel, list( ModelODEBolus, PFIMProject ) ) = function( pkModel, pfimproject ) {
  eq = remapOdePkLibraryEquations( prop( pkModel, "modelEquations" ), pfimproject )
  set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ) ) )
}

method( definePKModel, list( ModelODEDoseInEquations, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}

method( definePKModel, list( ModelODEDoseNotInEquations, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}
