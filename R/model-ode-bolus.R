#' @title ModelODEBolus
#' @description ODE bolus model with doses in initial conditions.
#' @inheritParams ModelODE
#' @param modelODE \code{deSolve} right-hand side built at administration time.
#' @param doseEvent Bolus dose event table passed to \code{deSolve}.
#' @param solverInputs Reserved; not used for bolus initial-condition dosing.
#' @include ModelODE.R
#' @return An S7 object of class \code{ModelODEBolus}.
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
#' @return An S7 object of class \code{ModelODEDoseInEquations}.
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
#' @return An S7 object of class \code{ModelODEDoseNotInEquations}.
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
#   bolusIc    - dose added to compartment state at t = 0 (initial conditions);
#                later doses still go through doseEvent
#   doseInEq   - dose_* / t_* terms inside the RHS; windows are (t_start, t_end]
#                (note: open on the left - differs from infusion [t, t+Tinf))
#   doseEvent  - all mass via deSolve events; ICs zeroed for dosed compartments
# Event method: replace at t==0, add for later times (avoids double-counting IC).

#' Resolve administration mode for bolus ODE model classes.
#' @param model \code{ModelODE} subclass instance.
#' @return Character scalar mode identifier.
#' @noRd
#' @keywords internal
.odeBolusMode = function( model ) {
  if ( S7::S7_inherits( model, ModelODEBolus ) ) return( "bolusIc" )
  if ( S7::S7_inherits( model, ModelODEDoseInEquations ) ) return( "doseInEq" )
  if ( S7::S7_inherits( model, ModelODEDoseNotInEquations ) ) return( "doseEvent" )
  stop( "not a bolus ODE model", call. = FALSE )
}

#' Replace bare time variable in dose-in-equation formulas.
#' @param equations Named character vector of ODE equations.
#' @param outcomes Character vector of administered outcome names.
#' @return Character vector of transformed equations.
#' @noRd
#' @keywords internal
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

#' Build ODE wrapper and formula metadata for bolus modes.
#' @param model \code{ModelODE} object to update.
#' @param evaluation \code{PFIMProject} or \code{Evaluation} source object.
#' @param mode Character mode among bolus administration strategies.
#' @return Updated \code{ModelODE} object.
#' @noRd
#' @keywords internal
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
  prop( model, "variableNames" ) = variableNames
  if ( mode %in% c( "bolusIc", "doseInEq" ) )
    prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration
  prop( model, "wrapper" ) = .buildODEWrapper( equations, functionArguments )
  prop( model, "functionArguments" ) = functionArguments
  prop( model, "functionArgumentsSymbol" ) = map( functionArguments, as.symbol )
  .setModelOutputFormulas( model, outputs )
  model
}

#' Parse textual initial-condition expressions.
#' @param initialConditionsTmp List of initial-condition values or expressions.
#' @return List of parsed expressions keyed by state variable.
#' @noRd
#' @keywords internal
.parseInitialConditionExprs = function( initialConditionsTmp ) {
  lapply( initialConditionsTmp, function( ic ) {
    if ( is.character( ic ) ) parse( text = ic )[[ 1L ]] else NULL
  } )
}

#' Apply initial-condition expressions to dose events.
#' @param doseEvent Data frame of event doses.
#' @param parsedExprs Named list of parsed initial-condition expressions.
#' @param env Environment containing current parameter and dose values.
#' @return Updated \code{doseEvent} data frame.
#' @noRd
#' @keywords internal
.applyInitialConditionsToEvent = function( doseEvent, parsedExprs, env ) {
  # Sequential: later IC exprs may read prior dose_* assigned into env.
  reduce(
    seq_len( nrow( doseEvent ) ),
    function( de, iter ) {
      expr = parsedExprs[[ de$var[ iter ] ]]
      if ( is.null( expr ) ) return( de )
      assign( paste0( "dose_", de$var[ iter ] ), de$value[ iter ], envir = env )
      used = unique( all.names( expr ) )
      allowed = c( .pfimIcSafeNames, ls( envir = env, all.names = TRUE ) )
      bad = setdiff( used, allowed )
      if ( length( bad ) )
        .pfimStop(
          "initial condition '", de$var[ iter ], "' uses unknown name(s): ",
          paste( bad, collapse = ", " ), "."
        )
      de$value[ iter ] = eval( expr, envir = env )
      de
    },
    .init = doseEvent
  )
}

#' Build deSolve event table from arm administrations.
#' @param arm \code{Arm} object containing administrations.
#' @param samplings Numeric vector of simulation times.
#' @return Ordered data frame of event rows for \code{deSolve}.
#' @noRd
#' @keywords internal
.odeDoseEventFromArm = function( arm, samplings, model = NULL ) {
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
  } ) |> list_rbind() |> ( function( df ) {
    df = df[ order( df$time ), ]
    if ( !is.null( model ) && nrow( df ) ) {
      states = .pfimOdeStateNames( model )
      aliases = .pfimOutputStateAliases( model )
      df$var = vapply(
        as.character( df$var ),
        function( v ) .pfimResolveOutcomeToState( v, states, aliases ),
        character( 1L )
      )
    }
    df
  } )()
}

#' Build dose-window inputs for dose-in-equation mode.
#'
#' Window convention is \eqn{(t_{\mathrm{start}},\, t_{\mathrm{end}}]} (see RHS
#' \code{t > a[,1] & t <= a[,2]}). For \code{tau != 0}, consecutive grid points
#' form those windows; a single bolus uses a degenerate 1-row matrix.
#' @param arm \code{Arm} object containing administrations.
#' @param samplings Numeric vector of simulation times.
#' @return Named list of administration windows and doses by outcome.
#' @noRd
#' @keywords internal
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
#' Build standard ODE RHS closure for deSolve integration.
#' @param wrapper Compiled model wrapper function.
#' @param functionArguments Character vector of wrapper argument names.
#' @param functionArgumentsSymbols Symbol list aligned with wrapper arguments.
#' @param outputFormula Named list of parsed output expressions.
#' @param mu Named numeric vector of parameter means.
#' @return Function suitable as \code{deSolve::ode} RHS.
#' @noRd
#' @keywords internal
.odeRhsStandard = function(
    wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu,
    variableNames ) {
  function( t, y, parms ) {
    if ( length( variableNames ) )
      y = set_names( as.numeric( y ), variableNames )
    with( as.list( c( t = t, y, mu ) ), {
      evaluationModel   = do.call( wrapper, set_names( functionArgumentsSymbols, functionArguments ) )
      evaluationOutputs = map( outputFormula, ~ eval( .x ) )
      c( evaluationModel, evaluationOutputs )
    } )
  }
}

# doseInEq RHS: for each outcome, locate the active infusion window and expose
# dose_<outcome> and t_<outcome> (time since window start) to the wrapper.
#' Build dose-aware RHS closure for dose-in-equation models.
#' @param wrapper Compiled model wrapper function.
#' @param functionArguments Character vector of wrapper argument names.
#' @param functionArgumentsSymbols Symbol list aligned with wrapper arguments.
#' @param outputFormula Named list of parsed output expressions.
#' @param mu Named numeric vector of parameter means.
#' @param outcomesWithAdministration Character outcomes receiving doses.
#' @return Function suitable as \code{deSolve::ode} RHS.
#' @noRd
#' @keywords internal
.odeRhsDoseInEq = function(
    wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu,
    outcomesWithAdministration, variableNames ) {
  # Precompute names once; RHS itself stays allocation-light every ODE step.
  dose_names = paste0( "dose_", outcomesWithAdministration )
  t_names    = paste0( "t_", outcomesWithAdministration )
  function( t, y, parms ) {
    doseTimeVars = flatten( map( seq_along( outcomesWithAdministration ), function( i ) {
      outcome   = outcomesWithAdministration[[ i ]]
      adminTime = parms[[ outcome ]]$administrationTime
      idx       = which( t > adminTime[ , 1 ] & t <= adminTime[ , 2 ] )
      set_names(
        list(
          if ( length( idx ) > 0 ) parms[[ outcome ]]$dose[ idx ] else parms[[ outcome ]]$dose[ 1 ],
          if ( length( idx ) > 0 ) t - adminTime[ idx, 1 ] else t
        ),
        c( dose_names[[ i ]], t_names[[ i ]] )
      )
    } ) )
    if ( length( variableNames ) )
      y = set_names( as.numeric( y ), variableNames )
    with( as.list( c( list( t = t ), as.list( y ), as.list( mu ), doseTimeVars ) ), {
      evaluationModel   = do.call( wrapper, set_names( functionArgumentsSymbols, functionArguments ) )
      evaluationOutputs = map( outputFormula, ~ eval( .x ) )
      c( evaluationModel, evaluationOutputs )
    } )
  }
}

#' Finalize administration slots for ODE model evaluation.
#' @param model \code{ModelODE} object to mutate.
#' @param mode Character mode among bolus administration strategies.
#' @param mu Named numeric parameter vector.
#' @param initialConditions Named numeric initial conditions.
#' @param samplings Numeric simulation time grid.
#' @param wrapper Compiled model wrapper function.
#' @param functionArguments Character wrapper argument names.
#' @param functionArgumentsSymbols Symbol list aligned to arguments.
#' @param outputFormula Named list of parsed output expressions.
#' @param doseEvent Optional event table for event-based modes.
#' @param solverInputs Optional dose-window inputs for dose-in-equation mode.
#' @param outcomesWithAdministration Optional administered outcomes.
#' @return Updated \code{ModelODE} object.
#' @noRd
#' @keywords internal
.odeFinalizeAdministration = function(
    model, mode, mu, initialConditions, samplings,
    wrapper, functionArguments, functionArgumentsSymbols, outputFormula,
    doseEvent = NULL, solverInputs = NULL, outcomesWithAdministration = NULL ) {

  variableNames = .pfimOdeStateNames( model )
  initialConditions = .pfimAlignOdeStates( initialConditions, variableNames )
  rhs = if ( mode == "doseInEq" ) {
    .odeRhsDoseInEq(
      wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu,
      outcomesWithAdministration, variableNames
    )
  } else {
    .odeRhsStandard(
      wrapper, functionArguments, functionArgumentsSymbols, outputFormula, mu,
      variableNames
    )
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
#' Evaluate initial conditions for bolus-at-time-zero mode.
#' @param model \code{ModelODEBolus} object.
#' @param arm \code{Arm} object with initial condition definitions.
#' @param doseEvent Dose-event data frame used to expose \code{dose_} variables.
#' @return Named numeric vector of initial compartment values.
#' @noRd
#' @keywords internal
.odeInitialConditionsBolus = function( model, arm, doseEvent ) {
  outcomes = prop( model, "outcomesWithAdministration" )
  dosesAtT0 = doseEvent[ doseEvent$time == 0, , drop = FALSE ]
  doses = set_names(
    as.list( dosesAtT0$value[ match( outcomes, dosesAtT0$var ) ] ),
    paste0( "dose_", outcomes )
  )
  mu  = .extractMu( prop( model, "modelParameters" ) )
  env = .odeIcEvalEnv( mu, doses )
  imap( prop( arm, "initialConditions" ), function( ic, compartment ) {
    if ( is.numeric( ic ) ) ic else .pfimEvalIcExpr( ic, env, compartment )
  } ) |> unlist()
}

# doseEvent: zero compartment states at t = 0; all mass via deSolve events table.
#' Evaluate initial conditions for event-based bolus mode.
#' @param model \code{ModelODEDoseNotInEquations} object.
#' @param arm \code{Arm} object with initial condition definitions.
#' @param doseEvent Dose-event data frame.
#' @return Named numeric vector of initial compartment values.
#' @noRd
#' @keywords internal
.odeInitialConditionsDoseEvent = function( model, arm, doseEvent ) {
  states = .pfimOdeStateNames( model )
  ic = .pfimAlignOdeStates( evaluateInitialConditions( model, arm ), states )
  dosed = intersect(
    states,
    unique( as.character( doseEvent$var[ doseEvent$time == 0 ] ) )
  )
  if ( length( dosed ) )
    ic[ dosed ] = 0
  ic
}

# Bind arm dosing to model slots; branches on bolusIc / doseInEq / doseEvent.
#' Bind arm dosing and sampling to an ODE model instance.
#' @param model \code{ModelODE} object to configure.
#' @param arm \code{Arm} object providing administration schedules.
#' @param mode Character mode among bolus administration strategies.
#' @return Updated \code{ModelODE} object.
#' @noRd
#' @keywords internal
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
    doseEvent = .odeDoseEventFromArm( arm, samplings, model )
    initialConditions = if ( mode == "bolusIc" ) {
      ic = .odeInitialConditionsBolus( model, arm, doseEvent )
      env = .odeIcEvalEnv( as.list( mu ) )
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

#' Map each output name to its sampling-time vector from the arm.
#'
#' Resolves via output-formula targets when the sampling outcome differs from
#' the reported output name.
#' @param samplingTimes List of \code{SamplingTimes} objects.
#' @param outputNames Character vector of model outputs.
#' @param outputFormula Optional output-formula mapping.
#' @return Named list of numeric sampling vectors.
#' @noRd
#' @keywords internal
.odeSamplingsByOutput = function( samplingTimes, outputNames, outputFormula = list() ) {
  samp_outcomes = map_chr( samplingTimes, ~ prop( .x, "outcome" ) )
  set_names(
    map( outputNames, function( out_name ) {
      target = out_name
      if ( length( outputFormula ) && out_name %in% names( outputFormula ) ) {
        x = outputFormula[[ out_name ]]
        if ( is.character( x ) && length( x ) == 1L ) target = x
        else if ( is.symbol( x ) || is.name( x ) ) target = as.character( x )
        else target = tryCatch( as.character( x ), error = function( e ) out_name )
      }
      idx = match( target, samp_outcomes )
      if ( is.na( idx ) ) idx = match( out_name, samp_outcomes )
      if ( is.na( idx ) ) idx = 1L
      prop( samplingTimes[[ idx ]], "samplings" )
    } ),
    outputNames
  )
}

#' Extract outputs at requested sampling times.
#' @param evaluationModelTmp Data frame returned by ODE integration.
#' @param samplingTimes List of \code{SamplingTimes} objects.
#' @param outputNames Character vector of requested outputs.
#' @param outputFormula Optional parsed output formulas.
#' @return Named list of output data frames by outcome.
#' @noRd
#' @keywords internal
.odeExtractOutputAtSamplingTimes = function(
    evaluationModelTmp, samplingTimes, outputNames, outputFormula = list() ) {
  samplings_list = .odeSamplingsByOutput( samplingTimes, outputNames, outputFormula )
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

# deSolve call (three bolus modes); events for doseEvent/bolusIc only.
#' Simulate ODE trajectories for bolus administration modes.
#' @param model \code{ModelODE} object with prepared administration fields.
#' @param arm \code{Arm} object used for sampling metadata.
#' @param mode Character mode among bolus administration strategies.
#' @return Data frame of integrated states and outputs.
#' @noRd
#' @keywords internal
.odeSimulateBolus = function( model, arm, mode ) {
  odeSolverParameters = prop( model, "odeSolverParameters" )
  raw_samplings       = prop( model, "samplings" )
  events_df           = if ( mode != "doseInEq" ) prop( model, "doseEvent" ) else NULL
  sim_key             = .pfimOdeSimCacheKey( mode, raw_samplings, events_df )
  sim_times = .pfimOdeSimTimesCached( sim_key, raw_samplings, events_df )
  ode_fn    = if ( mode == "doseInEq" )
    prop( model, "modelODEDoseInEquations" ) else prop( model, "modelODE" )
  parms = if ( mode == "doseInEq" ) prop( model, "solverInputs" ) else NULL
  tol   = .pfimDeSolveTolerances( odeSolverParameters )

  ode(
    prop( model, "initialConditions" ),
    sim_times,
    ode_fn,
    parms,
    events = if ( !is.null( events_df ) ) list( data = events_df ) else NULL,
    atol   = tol$atol,
    rtol   = tol$rtol
  ) |> as.data.frame()
}

#' Evaluate bolus ODE model outputs at arm sampling times.
#' @param model \code{ModelODE} object with prepared administration state.
#' @param arm \code{Arm} object providing requested sampling times.
#' @return Named list of output data frames.
#' @noRd
#' @keywords internal
.odeEvaluateModelCore = function( model, arm ) {
  mode = .odeBolusMode( model )
  evaluationModelTmp = .odeSimulateBolus( model, arm, mode )
  .odeExtractOutputAtSamplingTimes(
    evaluationModelTmp,
    prop( arm, "samplingTimes" ),
    prop( model, "outputNames" ),
    .getModelOutputFormulas( model )
  )
}

#' Evaluate bolus ODE model with optional covariate expansion.
#' @param model \code{ModelODE} object.
#' @param arm \code{Arm} object used for evaluation.
#' @return Output structure from direct or covariate-expanded evaluation.
#' @noRd
#' @keywords internal
.odeEvaluateBolus = function( model, arm ) {
  if ( usesCovariateOccasionStructure( model ) )
    evaluateModelWithCovariates( model, arm, .odeEvaluateModelCore )
  else
    .odeEvaluateModelCore( model, arm )
}

# Fast path for gradient cache: reapply stored administration template to new mu.
#' Reapply cached ODE administration entry to current parameters.
#' @param model \code{ModelODE} object to update.
#' @param arm \code{Arm} object used for initial-condition context.
#' @param entry Cached administration entry list.
#' @return Updated \code{ModelODE} object.
#' @noRd
#' @keywords internal
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
      env = .odeIcEvalEnv( as.list( mu ) )
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

#' Compile ODE wrapper for bolus-in-initial-condition mode.
#' @return Updated \code{ModelODEBolus} wrapper configuration.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelODEBolus ) = function( model, evaluation ) {
  .odeDefineWrapper( model, evaluation, "bolusIc" )
}
#' Compile ODE wrapper for dose-in-equation bolus mode.
#' @return Updated \code{ModelODEDoseInEquations} wrapper configuration.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelODEDoseInEquations ) = function( model, evaluation ) {
  .odeDefineWrapper( model, evaluation, "doseInEq" )
}
#' Compile ODE wrapper for event-based bolus mode.
#' @return Updated \code{ModelODEDoseNotInEquations} wrapper configuration.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelODEDoseNotInEquations ) = function( model, evaluation ) {
  .odeDefineWrapper( model, evaluation, "doseEvent" )
}

#' Evaluate bolus IC expressions with \code{dose_*} exposed at t = 0.
#' @return Named numeric vector of evaluated initial conditions.
#' @name evaluateInitialConditions
#' @keywords internal
method( evaluateInitialConditions, ModelODEBolus ) = function( model, arm, doseEvent ) {
  .odeInitialConditionsBolus( model, arm, doseEvent )
}

#' Bind arm dosing for bolus-in-initial-condition mode.
#' @return Updated model with bolus-in-initial-condition administration setup.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelODEBolus ) = function( model, arm ) {
  .odeDefineAdministration( model, arm, "bolusIc" )
}
#' Bind arm dosing for dose-in-equation bolus mode.
#' @return Updated model with dose-in-equation administration setup.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelODEDoseInEquations ) = function( model, arm ) {
  .odeDefineAdministration( model, arm, "doseInEq" )
}
#' Bind arm dosing for event-based bolus mode.
#' @return Updated model with event-based administration setup.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelODEDoseNotInEquations ) = function( model, arm ) {
  .odeDefineAdministration( model, arm, "doseEvent" )
}

#' Dispatch: covariate/occasion structure -> specialised path; else core evaluator.
#' @return Model outputs at arm sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelODEBolus ) = function( model, arm ) {
  .odeEvaluateBolus( model, arm )
}
#' Dispatch: covariate/occasion structure -> specialised path; else core evaluator.
#' @return Model outputs at arm sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelODEDoseInEquations ) = function( model, arm ) {
  .odeEvaluateBolus( model, arm )
}
#' Dispatch: covariate/occasion structure -> specialised path; else core evaluator.
#' @return Model outputs at arm sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelODEDoseNotInEquations ) = function( model, arm ) {
  .odeEvaluateBolus( model, arm )
}

#' Remap library ODE PK equations onto project compartments (bolus IC).
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} used for compartment remapping.
#' @return Named list of remapped PK equations.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelODEBolus, PFIMProject ) ) = function( pkModel, pfimproject ) {
  eq = remapOdePkLibraryEquations( prop( pkModel, "modelEquations" ), pfimproject )
  set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ), names( eq ) ) )
}

#' Return PK equations stored on the dose-in-equation model.
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List of PK equations from \code{pkModel}.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelODEDoseInEquations, PFIMProject ) ) = function( pkModel, pfimproject ) {
  eq = remapOdePkLibraryEquations( prop( pkModel, "modelEquations" ), pfimproject )
  set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ), names( eq ) ) )
}

#' Return PK equations stored on the dose-event model.
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List of PK equations from \code{pkModel}.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelODEDoseNotInEquations, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}
