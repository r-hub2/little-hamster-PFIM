# FIM design cache: keyed by project scope + arm doses/sampling signature (optimization grid).

.pfimFimDesignCache = new.env( parent = emptyenv() )

#' Copy FIM results from a finished evaluationDesign back into an Evaluation shell.
#' @keywords internal
.pfimEvaluationFromDesign = function( evaluation, design, evaluatedDesign ) {
  prop( evaluation, "designs" ) = list( design )
  prop( evaluation, "evaluationDesign" ) = list( evaluatedDesign )
  prop( evaluation, "fim" ) = prop( evaluatedDesign, "fim" )
  evaluation
}

#' Content hash for rebuildEvalModel cache (mu/omega, library, covariate count).
#' Falls back to object address when PFIM.model.cache.signature = FALSE.
#' @keywords internal
.pfimModelCacheId = function( pfimproject ) {
  if ( !isTRUE( pfim_get_option( "model.cache.signature", TRUE ) ) )
    return( rlang::obj_address( projectOf( pfimproject ) ) )

  params = projectProp( pfimproject, "modelParameters" )
  mu_sig = paste(
    vapply( params, function( p ) {
      d = prop( p, "distribution" )
      paste(
        prop( p, "name" ),
        round( prop( d, "mu" ), 8L ),
        round( prop( d, "omega" ), 8L ),
        sep = ":"
      )
    }, character( 1L ) ),
    collapse = ";"
  )
  lib = projectProp( pfimproject, "modelFromLibrary" )
  lib_key = if ( length( lib ) )
    paste( names( lib ), unlist( lib ), sep = "=", collapse = "," ) else ""
  cov_eq = projectProp( pfimproject, "modelCovariatesEquation" )
  paste(
    projectProp( pfimproject, "fimType" ),
    lib_key,
    mu_sig,
    length( projectProp( pfimproject, "modelCovariates" ) ),
    cov_eq,
    sep = "|"
  )
}

#' Doses and sampling times flattened per arm — used only for cache lookup.
#' @keywords internal
.pfimDesignSignature = function( design ) {
  arms = prop( design, "arms" )
  parts = vapply( arms, function( arm ) {
    adms = prop( arm, "administrations" )
    doses = vapply( adms, function( a ) {
      d = prop( a, "dose" )
      paste( as.numeric( d ), collapse = "," )
    }, character( 1L ) )
    sts = prop( arm, "samplingTimes" )
    samps = vapply( sts, function( s ) {
      paste( as.numeric( prop( s, "samplings" ) ), collapse = "," )
    }, character( 1L ) )
    paste( c( doses, samps ), collapse = ";" )
  }, character( 1L ) )
  paste( parts, collapse = "|" )
}

#' @keywords internal
.pfimFimCacheKey = function( evaluation ) {
  scope = pfim_get_option( "fim.cache.scope", NULL )
  if ( is.null( scope ) )
    scope = rlang::obj_address( projectOf( evaluation ) )
  design = prop( evaluation, "designs" )[[ 1L ]]
  paste0( scope, "::", .pfimDesignSignature( design ) )
}

#' Start a fresh cache namespace for one optimization run.
#' @keywords internal
.pfimFimCacheBegin = function( optimization ) {
  scope = rlang::obj_address( projectOf( optimization ) )
  pfim_set_option( fim.cache.scope = scope )
  .pfimFimCacheClear( scope )
  invisible( scope )
}

#' @keywords internal
.pfimFimCacheClear = function( scope ) {
  prefix = paste0( scope, "::" )
  keys = ls( .pfimFimDesignCache, all.names = TRUE )
  rm( list = keys[ startsWith( keys, prefix ) ], envir = .pfimFimDesignCache )
  invisible( NULL )
}

#' @keywords internal
.pfimFimCacheRegister = function( evaluation ) {
  if ( !isTRUE( pfim_get_option( "fim.cache", TRUE ) ) )
    return( invisible( NULL ) )
  key = .pfimFimCacheKey( evaluation )
  assign( key, evaluation, envir = .pfimFimDesignCache )
  invisible( key )
}

#' @keywords internal
.pfimRunEvaluationCached = function( evaluation ) {
  if ( !isTRUE( pfim_get_option( "fim.cache", TRUE ) ) )
    return( run( evaluation ) )

  key = .pfimFimCacheKey( evaluation )
  if ( exists( key, envir = .pfimFimDesignCache, inherits = FALSE ) ) {
    hits = pfim_get_option( "fim.cache.hits", 0L )
    pfim_set_option( fim.cache.hits = hits + 1L )
    return( get( key, envir = .pfimFimDesignCache ) )
  }

  result = run( evaluation )
  assign( key, result, envir = .pfimFimDesignCache )
  result
}

#' @keywords internal
.pfimFimCacheStats = function() {
  list(
    enabled = isTRUE( pfim_get_option( "fim.cache", TRUE ) ),
    scope   = pfim_get_option( "fim.cache.scope", NA ),
    size    = length( ls( .pfimFimDesignCache, all.names = TRUE ) ),
    hits    = pfim_get_option( "fim.cache.hits", 0L )
  )
}

#' Metaheuristic cost = 1 / D-criterion (kernels minimise this).
#' Singular FIM → .metaheuristicFitnessPenalty so the particle is discarded.
#' @keywords internal
.pfimMetaheuristicFitness = function( evalTemplate, design, arms, flat_pos ) {
  tempDesign = design
  prop( tempDesign, "arms" ) = .applyFlatToArms( flat_pos, arms )
  prop( evalTemplate, "designs" ) = list( tempDesign )
  fim = prop( .pfimRunEvaluationCached( evalTemplate ), "fim" )
  d   = 1 / Dcriterion( fim )
  if ( !is.finite( d ) ) .metaheuristicFitnessPenalty else d
}

#' Batch fitness for PSO (matrix of flat positions, one cost per row).
#'
#' Invalid rows (when \code{layout} is set) get the metaheuristic penalty without
#' running the FIM.
#' @keywords internal
.pfimMetaheuristicFitnessBatch = function( evalTemplate, design, arms, pos_matrix,
                                           layout = NULL ) {
  pos_matrix = as.matrix( pos_matrix )
  n = nrow( pos_matrix )
  if ( n == 0L ) return( numeric( 0 ) )
  if ( n == 1L )
    return( .pfimMetaheuristicFitness( evalTemplate, design, arms, pos_matrix[ 1L, ] ) )

  map(
    seq_len( n ),
    function( i ) {
      flat = pos_matrix[ i, ]
      if ( !is.null( layout ) && !.isFlatValid( layout, flat ) )
        return( .metaheuristicFitnessPenalty )
      .pfimMetaheuristicFitness( evalTemplate, design, arms, flat )
    }
  ) |> unlist( use.names = FALSE )
}

#' Run multiple evaluations through the FIM cache.
#' @keywords internal
.pfimRunEvaluations = function( evaluations ) {
  if ( !length( evaluations ) ) return( evaluations )
  map( seq_along( evaluations ), function( i ) .pfimRunEvaluationCached( evaluations[[ i ]] ) )
}

#' Run a final optimization design evaluation through the FIM cache.
#' @keywords internal
.pfimRunOptimizationEvaluation = function( optimization, design, name = "" ) {
  .pfimRunEvaluationCached(
    .evaluationFromOptimization( optimization, design, name = name )
  )
}
