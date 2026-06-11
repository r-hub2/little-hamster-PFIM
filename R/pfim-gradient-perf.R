# Performance layer for gradient / ODE evaluation.
#
# Three independent caches (all optional via pfim_set_option):
#   perf.fdCache        — finite-difference step grid when mu is unchanged
#   perf.adminCache     — ODE administration objects across FD perturbations
#   perf.odeTimesCache  — deSolve time grid (samplings + dose event times)
#
# Invalidation: .pfimClearGradientPerfCaches() clears all env caches (tests).

.pfimFdSchemeCache   = new.env( parent = emptyenv() )
.pfimGradAdminCache  = new.env( parent = emptyenv() )
.pfimOdeSimGrid      = new.env( parent = emptyenv() )

.pfimPerfOption = function( name, default = FALSE ) {
  isTRUE( pfim_get_option( name, default ) )
}

# Cache key for FD scheme: parameter names + formatted mu (12-digit sci notation).
.pfimMuSignature = function( parameters ) {
  paste(
    map_chr( parameters, ~ prop( .x, "name" ) ),
    format(
      map_dbl( parameters, ~ prop( prop( .x, "distribution" ), "mu" ) ),
      digits = 12, scientific = TRUE
    ),
    sep = "=",
    collapse = "|"
  )
}

# Arm signature: doses, sampling times, initial-condition expressions.
.pfimArmAdminSignature = function( arm ) {
  adms = prop( arm, "administrations" )
  sts  = prop( arm, "samplingTimes" )
  ic   = prop( arm, "initialConditions" )

  dose_sig = vapply( adms, function( adm ) {
    dosing = .alignAdministrationDosing( adm )
    paste(
      prop( adm, "outcome" ),
      paste( dosing$timeDose, collapse = "," ),
      paste( dosing$dose, collapse = "," ),
      prop( adm, "tau" ),
      sep = ":"
    )
  }, character( 1L ) )

  samp_sig = vapply( sts, function( st ) {
    paste( prop( st, "outcome" ), paste( prop( st, "samplings" ), collapse = "," ), sep = ":" )
  }, character( 1L ) )

  ic_sig = if ( length( ic ) )
    paste( names( ic ), vapply( ic, as.character, character( 1L ) ), sep = "=", collapse = ";" )
  else ""

  paste( paste( dose_sig, collapse = ";" ), paste( samp_sig, collapse = ";" ), ic_sig, sep = "||" )
}

.pfimGradAdminCacheKey = function( model, arm ) {
  cls = class( model )[[ 1L ]]
  param_names = paste( map_chr( prop( model, "modelParameters" ), ~ prop( .x, "name" ) ), collapse = "," )
  paste( cls, param_names, .pfimArmAdminSignature( arm ), sep = "::" )
}

# Apply one column of the FD Hessian perturbation to mu on each ModelParameter.
.pfimShiftModelParameters = function( parameters, shiftedCol ) {
  map2( parameters, shiftedCol, function( param, newMu ) {
    distr = prop( param, "distribution" )
    prop( distr, "mu" ) = newMu
    prop( param, "distribution" ) = distr
    param
  } )
}

# Reuse parametersForComputingGradient when mu signature matches (FD loop hot path).
.pfimFiniteDifferenceHessianCached = function( model ) {
  if ( !.pfimPerfOption( "PFIM.perf.fdCache", TRUE ) )
    return( finiteDifferenceHessian( model ) )

  key = .pfimMuSignature( prop( model, "modelParameters" ) )
  cached = .pfimFdSchemeCache[[ key ]]
  if ( !is.null( cached ) ) {
    prop( model, "parametersForComputingGradient" ) = cached
    return( model )
  }

  model = finiteDifferenceHessian( model )
  .pfimFdSchemeCache[[ key ]] = prop( model, "parametersForComputingGradient" )
  model
}

.pfimIsOdeModelClass = function( model ) {
  S7::S7_inherits( model, ModelODE )
}

# Sorted unique simulation times; drops consecutive duplicates within 1e-8.
.pfimOdeSimTimesCached = function( cacheKey, raw_samplings, events_df ) {
  if ( .pfimPerfOption( "PFIM.perf.odeTimesCache", TRUE ) &&
       !is.null( .pfimOdeSimGrid[[ cacheKey ]] ) )
    return( .pfimOdeSimGrid[[ cacheKey ]] )

  sim_times = sort( unique( c( raw_samplings, if ( !is.null( events_df ) ) events_df$time else numeric( 0 ) ) ) )
  sim_times = sim_times[ c( TRUE, diff( sim_times ) > 1e-8 ) ]

  if ( .pfimPerfOption( "PFIM.perf.odeTimesCache", TRUE ) )
    .pfimOdeSimGrid[[ cacheKey ]] = sim_times

  sim_times
}

# First FD call builds full administration; later calls refresh mu/dose only.
.pfimDefineModelAdministrationCached = function( model, arm ) {
  if ( !.pfimPerfOption( "PFIM.perf.adminCache", TRUE ) || !.pfimIsOdeModelClass( model ) )
    return( defineModelAdministration( model, arm ) )

  key   = .pfimGradAdminCacheKey( model, arm )
  entry = .pfimGradAdminCache[[ key ]]

  if ( is.null( entry ) ) {
    model = defineModelAdministration( model, arm )
    .pfimGradAdminCache[[ key ]] = .pfimCaptureOdeAdminEntry( model, arm )
    return( model )
  }

  .pfimRefreshOdeAdministration( model, arm, entry )
}

# Snapshot of parsed wrapper/RHS/dose tables — enough to rebuild after mu shift.
.pfimCaptureOdeAdminEntry = function( model, arm ) {
  if ( !S7::S7_inherits( model, ModelODE ) ) return( list( type = "full" ) )
  mode = tryCatch( .odeBolusMode( model ), error = function( e ) "full" )
  if ( identical( mode, "full" ) ) return( list( type = "full" ) )

  entry = list(
    type                       = mode,
    samplings                  = prop( model, "samplings" ),
    wrapper                    = prop( model, "wrapper" ),
    functionArguments          = prop( model, "functionArguments" ),
    functionArgumentsSymbols   = prop( model, "functionArgumentsSymbol" ),
    outputFormula              = .getOutputFormulaParsed( model )
  )
  if ( mode == "doseInEq" ) {
    entry$solverInputs               = prop( model, "solverInputs" )
    entry$outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  } else {
    entry$doseEventTemplate = .odeDoseEventFromArm( arm, prop( model, "samplings" ) )
    if ( mode == "bolusIc" )
      entry$initialConditionsParsed = .parseInitialConditionExprs( prop( arm, "initialConditions" ) )
  }
  entry
}

.pfimRefreshOdeAdministration = function( model, arm, entry ) {
  if ( entry$type == "full" )
    defineModelAdministration( model, arm )
  else
    .odeApplyAdminEntry( model, arm, entry )
}

.pfimPrepareModelForEvaluation = function( model, arm ) {
  .pfimDefineModelAdministrationCached( model, arm )
}

.pfimClearGradientPerfCaches = function() {
  rm( list = ls( .pfimFdSchemeCache, all.names = TRUE ), envir = .pfimFdSchemeCache )
  rm( list = ls( .pfimGradAdminCache, all.names = TRUE ), envir = .pfimGradAdminCache )
  rm( list = ls( .pfimOdeSimGrid, all.names = TRUE ), envir = .pfimOdeSimGrid )
  rm( list = ls( .pfimOutputFormulaCache, all.names = TRUE ),
      envir = .pfimOutputFormulaCache )
  rm( list = ls( .pfimParsedOutputFormulaCache, all.names = TRUE ),
      envir = .pfimParsedOutputFormulaCache )
  invisible( NULL )
}

.pfimGradientPerfStatus = function() {
  list(
    adminCache        = .pfimPerfOption( "PFIM.perf.adminCache", TRUE ),
    fdCache           = .pfimPerfOption( "PFIM.perf.fdCache", TRUE ),
    odeTimesCache     = .pfimPerfOption( "PFIM.perf.odeTimesCache", TRUE ),
    fdCacheSize       = length( ls( .pfimFdSchemeCache, all.names = TRUE ) ),
    adminCacheSize    = length( ls( .pfimGradAdminCache, all.names = TRUE ) ),
    odeTimesCacheSize = length( ls( .pfimOdeSimGrid, all.names = TRUE ) )
  )
}
