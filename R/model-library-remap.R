# PK/PD library equation remapping (RespPK, dose_, Tinf_ tokens must stay intact).
#' Replace variable names in library model equations.
#'
#' Replaces word-boundary occurrences of \code{old} with \code{new} in \code{text},
#' skipping protected pharmacokinetic tokens (\code{dose_RespPK}, \code{Tinf_RespPK},
#' and names containing \code{Emax}).
#'
#' @param text Character string: model equation or expression.
#' @param old  Character string: variable name to replace.
#' @param new  Character string: replacement variable name.
#' @return Updated character string.
#' @export

replaceVariablesLibraryOfModels = function( text, old, new ) {
  protectedTerms = c( "dose_", "Tinf_", "Emax" )

  if ( any( str_detect( old, fixed( protectedTerms ) ) ) ) return( text )

  if ( old == "RespPK" ) {
    str_replace_all(
      text,
      regex( paste0( "(?<!dose_)(?<!Tinf_)", old, "(?!\\w)" ) ),
      new
    )
  } else {
    str_replace_all( text, regex( paste0( "\\b", old, "\\b" ) ), new )
  }
}

#' @keywords internal
.pkpdLibraryCompartments = function( pfimproject ) {
  outputs = prop( pfimproject, "outputs" )

  pkComp = NULL
  pdComp = NULL
  if ( length( outputs ) && "RespPK" %in% names( outputs ) ) {
    pkComp = unlist( outputs[["RespPK"]], use.names = FALSE )
    if ( !length( pkComp ) || !nzchar( pkComp ) ) pkComp = NULL
  }
  if ( length( outputs ) >= 2L && "RespPD" %in% names( outputs ) ) {
    pdComp = unlist( outputs[["RespPD"]], use.names = FALSE )
    if ( !length( pdComp ) || !nzchar( pdComp ) ) pdComp = NULL
  }
  if ( !is.null( pkComp ) && !is.null( pdComp ) )
    return( c( RespPK = pkComp, E = pdComp ) )
  if ( !is.null( pkComp ) && is.null( pdComp ) )
    return( c( RespPK = pkComp ) )

  icNames = prop( pfimproject, "designs" ) |>
    map( \(d) map( prop( d, "arms" ), \(arm) names( prop( arm, "initialConditions" ) ) ) ) |>
    unlist() |>
    unique()

  if ( length( icNames ) >= 2L )
    return( set_names( icNames[1:2], c( "RespPK", "E" ) ) )
  if ( length( icNames ) == 1L )
    return( set_names( icNames, "RespPK" ) )

  c( RespPK = "C1", E = "C2" )
}

#' @keywords internal
.remapLibraryEquations = function( equations, remapText ) {
  if ( is.list( equations ) && !is.null( equations$duringInfusion ) ) {
    return( list(
      duringInfusion = map( equations$duringInfusion, remapText ),
      afterInfusion  = map( equations$afterInfusion,  remapText )
    ) )
  }
  map( equations, remapText )
}

#' @keywords internal
.remapLibraryDoseTokens = function( text, pfimproject ) {
  outputs = prop( pfimproject, "outputs" )
  if ( !length( outputs ) || !"RespPK" %in% names( outputs ) ) return( text )

  pk_comp = unlist( outputs[["RespPK"]], use.names = FALSE )
  if ( !length( pk_comp ) || !nzchar( pk_comp ) || identical( pk_comp, "RespPK" ) )
    return( text )

  admin_outcomes = .getOutcomesFromEvaluation( pfimproject )
  if ( !pk_comp %in% admin_outcomes ) return( text )

  text = str_replace_all( text, fixed( "dose_RespPK" ), paste0( "dose_", pk_comp ) )
  str_replace_all( text, fixed( "Tinf_RespPK" ), paste0( "Tinf_", pk_comp ) )
}

#' @keywords internal
remapPkpdLibraryText = function( text, pfimproject ) {
  comp = .pkpdLibraryCompartments( pfimproject )
  text = reduce2( names( comp ), unname( comp ), replaceVariablesLibraryOfModels, .init = text )
  .remapLibraryDoseTokens( text, pfimproject )
}

#' @keywords internal
remapPkpdLibraryEquations = function( equations, pfimproject ) {
  .remapLibraryEquations( equations, \(text) remapPkpdLibraryText( text, pfimproject ) )
}

#' @keywords internal
.derivativeNamesFromCompartments = function( pfimproject, nEq ) {
  comp = unname( .pkpdLibraryCompartments( pfimproject ) )
  paste0( "Deriv_", comp[ seq_len( nEq ) ] )
}

#' @keywords internal
.pkpdOdePkCompartments = function( pfimproject ) {
  comp = .pkpdLibraryCompartments( pfimproject )
  out = c( C1 = unname( comp["RespPK"] ) )
  if ( length( comp ) >= 2L ) out = c( out, C2 = unname( comp["E"] ) )
  out
}

#' @keywords internal
remapOdePkLibraryText = function( text, pfimproject ) {
  comp = .pkpdOdePkCompartments( pfimproject )
  reduce2( names( comp ), unname( comp ), replaceVariablesLibraryOfModels, .init = text )
}

#' @keywords internal
remapOdePkLibraryEquations = function( equations, pfimproject ) {
  .remapLibraryEquations( equations, \(text) remapOdePkLibraryText( text, pfimproject ) )
}

#' @keywords internal
.combineInfusionPkPdEquations = function( pkPart, pdEq, derivNames ) {
  stats::setNames(
    as.list( c(
      unlist( pkPart, use.names = FALSE ),
      unlist( pdEq, use.names = FALSE )
    ) ),
    derivNames
  )
}
