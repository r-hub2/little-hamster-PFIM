#' @title ModelODE
#' @description Base class for models integrated with \code{deSolve}.
#' @inheritParams Model
#' @include Model.R
#' @export

ModelODE = new_class( "ModelODE", package = "PFIM", parent = Model )

.scalar_finite = function( x ) {
  length( x ) == 1L && !is.na( x ) && is.finite( x )
}

.extractMu = function( parameters ) {
  set_names(
    map_dbl( parameters, function( parameter ) {
      distribution = prop( parameter, "distribution" )
      value        = prop( parameter, "value" )

      if ( !is.null( distribution ) ) {
        omega = prop( distribution, "omega" )
        mu    = prop( distribution, "mu" )
        if ( .scalar_finite( omega ) && omega != 0 )
          return( mu )
        if ( .scalar_finite( mu ) && mu != 0 && !.scalar_finite( value ) )
          return( mu )
      }
      if ( .scalar_finite( value ) )
        return( value )
      if ( !is.null( distribution ) && .scalar_finite( prop( distribution, "mu" ) ) )
        return( prop( distribution, "mu" ) )
      NA_real_
    }),
    map_chr( parameters, ~ prop( .x, "name" ) )
  )
}

.buildSamplings = function( samplingTimes ) {
  s = map( samplingTimes, ~ prop( .x, "samplings" ) ) |> unlist() |> sort() |> unique()
  c( 0.0, s[ s != 0 ] )
}

.getOutcomesFromEvaluation = function( evaluation ) {
  prop( evaluation, "designs" ) |>
    map( \(d) prop( d, "arms" ) ) |>
    list_flatten() |>
    map( \(arm) prop( arm, "administrations" ) ) |>
    list_flatten() |>
    map_chr( \(adm) prop( adm, "outcome" ) ) |>
    unique()
}

.administeredLibraryOutcomeNames = function( evaluation ) {
  adminOutcomes = .getOutcomesFromEvaluation( evaluation )
  outputs       = prop( evaluation, "outputs" )
  if ( length( outputs ) == 0L ) return( adminOutcomes )
  outputValues = unlist( outputs, use.names = TRUE )
  outputNames  = names( outputValues )
  if ( is.null( outputNames ) ) outputNames = as.character( outputValues )
  kept = outputValues %in% adminOutcomes
  if ( any( kept ) ) outputNames[ kept ] else adminOutcomes
}

.libraryEquationTimeName = function( evaluation, libraryOutcomeName ) {
  outputs = prop( evaluation, "outputs" )
  admin   = outputs[[ libraryOutcomeName ]]
  paste0( "t_", if ( !is.null( admin ) && nzchar( admin ) ) admin else libraryOutcomeName )
}

.libraryEquationTimeNames = function( evaluation, libraryOutcomeNames ) {
  map_chr( libraryOutcomeNames, ~ .libraryEquationTimeName( evaluation, .x ) )
}

.buildODEWrapper = function( equations, functionArguments ) {
  body = paste(
    paste( map_chr( names( equations ),
                    ~ sprintf( "%s = %s", .x, equations[[.x]] ) ),
           collapse = "\n" ),
    sprintf( "return(list(c(%s)))", paste( names( equations ), collapse = ", " ) ),
    sep = "\n"
  )
  eval( parse( text = sprintf( "function(%s) { %s }",
                               paste( functionArguments, collapse = ", " ), body ) ) )
}

.pfimOutputFormulaCache = new.env( parent = emptyenv() )
.pfimParsedOutputFormulaCache = new.env( parent = emptyenv() )

.parseOutputFormulaEntry = function( x ) {
  if ( is.character( x ) ) parse( text = x ) else x
}

.odeColumnForOutput = function( out_name, evaluationModelTmp, outputFormula ) {
  cols = names( evaluationModelTmp )
  if ( out_name %in% cols ) return( out_name )
  if ( length( outputFormula ) && out_name %in% names( outputFormula ) ) {
    target = outputFormula[[ out_name ]]
    if ( is.character( target ) && length( target ) == 1L ) {
      dotted = paste0( out_name, ".", target )
      if ( dotted %in% cols ) return( dotted )
      if ( target %in% cols ) return( target )
    }
  }
  out_name
}

.pfimOutputFormulaLookupKey = function( model ) {
  paste( class( model )[[ 1L ]], paste( prop( model, "outputNames" ), collapse = "," ), sep = "::" )
}

.setModelOutputFormulas = function( model, outputs ) {
  lookupKey = .pfimOutputFormulaLookupKey( model )
  .pfimOutputFormulaCache[[ lookupKey ]] = outputs
  .pfimParsedOutputFormulaCache[[ lookupKey ]] = NULL
  invisible( model )
}

.getModelOutputFormulas = function( model ) {
  cached = .pfimOutputFormulaCache[[ .pfimOutputFormulaLookupKey( model ) ]]
  if ( !is.null( cached ) ) return( cached )
  prop( model, "outputFormula" )
}

.getOutputFormulaParsed = function( model ) {
  key = .pfimOutputFormulaLookupKey( model )
  cached = .pfimParsedOutputFormulaCache[[ key ]]
  if ( !is.null( cached ) ) return( cached )
  outputs = .getModelOutputFormulas( model )
  if ( !length( outputs ) ) return( outputs )
  parsed = map( outputs, .parseOutputFormulaEntry )
  .pfimParsedOutputFormulaCache[[ key ]] = parsed
  parsed
}

.evalInitialConditionsImpl = function( model, arm ) {
  mu = .extractMu( prop( model, "modelParameters" ) )
  list2env( as.list( mu ), envir = environment() )
  imap( prop( arm, "initialConditions" ), function( ic, compartment ) {
    if ( is.numeric( ic ) )
      return( ic )
    eval( parse( text = ic ) )
  }) |> unlist()
}

method( evaluateInitialConditions, ModelODE ) = function( model, arm ) {
  .evalInitialConditionsImpl( model, arm )
}
