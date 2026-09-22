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
#' @keywords internal

replaceVariablesLibraryOfModels = function( text, old, new ) {
  protectedTerms = c( "dose_", "Tinf_", "Emax" )

  # Do not remap dose_/Tinf_/Emax* tokens themselves.
  if ( any( str_detect( old, fixed( protectedTerms ) ) ) ) return( text )

  if ( old == "RespPK" ) {
    # Avoid rewriting dose_RespPK / Tinf_RespPK when remapping RespPK alone.
    str_replace_all(
      text,
      regex( paste0( "(?<!dose_)(?<!Tinf_)", old, "(?!\\w)" ) ),
      new
    )
  } else {
    str_replace_all( text, regex( paste0( "\\b", old, "\\b" ) ), new )
  }
}

#' Resolve library compartment names (RespPK / E) from project outputs or ICs.
#'
#' Prefers \code{outputs$RespPK} / \code{outputs$RespPD}; falls back to the first
#' one or two initial-condition names, then defaults \code{C1}/\code{C2}.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Named character vector (\code{RespPK}, optionally \code{E}).
#' @noRd
#' @keywords internal
.pkpdLibraryCompartments = function( pfimproject ) {
  outputs = prop( pfimproject, "outputs" )

  pkComp = NULL
  pdComp = NULL
  if ( length( outputs ) && "RespPK" %in% names( outputs ) ) {
    pkComp = unlist( outputs[["RespPK"]], use.names = FALSE )
    if ( !length( pkComp ) || any( is.na( pkComp ) | !nzchar( pkComp ) ) ) pkComp = NULL
  }
  if ( length( outputs ) >= 2L && "RespPD" %in% names( outputs ) ) {
    pdComp = unlist( outputs[["RespPD"]], use.names = FALSE )
    if ( !length( pdComp ) || any( is.na( pdComp ) | !nzchar( pdComp ) ) ) pdComp = NULL
  }
  if ( !is.null( pkComp ) && !is.null( pdComp ) )
    return( c( RespPK = pkComp, E = pdComp ) )
  if ( !is.null( pkComp ) && is.null( pdComp ) )
    return( c( RespPK = pkComp ) )

  # Fallback: first IC name(s) from designs when outputs do not map compartments.
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

#' Apply a text remapper to flat or during/after infusion equation lists.
#' @param equations Character vector or list with \code{duringInfusion}/\code{afterInfusion}.
#' @param remapText Function \code{function(text) -> text}.
#' @return Remapped equations in the same structure as \code{equations}.
#' @noRd
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

#' Remap \code{dose_RespPK} / \code{Tinf_RespPK} onto the administered PK target.
#'
#' Library ODEs hard-code \code{dose_RespPK}. When the arm administers
#' \code{C1} (or another mapped compartment), rewrite the tokens so the ODE
#' wrapper's \code{dose_<admin>} formals match the equation.
#' @param text Character equation string.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Updated character string.
#' @noRd
#' @keywords internal
.pfimLibraryPkDoseTarget = function( pfimproject, admin_outcomes ) {
  outputs = prop( pfimproject, "outputs" )
  mapped = NULL
  if ( length( outputs ) && "RespPK" %in% names( outputs ) ) {
    mapped = unlist( outputs[[ "RespPK" ]], use.names = FALSE )
    if ( !.pfimIsNonEmptyScalar( mapped ) ) mapped = NULL
  }
  if ( !is.null( mapped ) && mapped %in% admin_outcomes )
    return( mapped )
  if ( "RespPK" %in% admin_outcomes )
    return( "RespPK" )
  if ( "C1" %in% admin_outcomes )
    return( "C1" )
  NULL
}

.remapLibraryDoseTokens = function( text, pfimproject ) {
  if ( !grepl( "dose_RespPK|Tinf_RespPK", text, perl = TRUE ) )
    return( text )
  admin_outcomes = .getOutcomesFromEvaluation( pfimproject )
  target = .pfimLibraryPkDoseTarget( pfimproject, admin_outcomes )
  if ( is.null( target ) || identical( target, "RespPK" ) )
    return( text )
  text = str_replace_all( text, fixed( "dose_RespPK" ), paste0( "dose_", target ) )
  str_replace_all( text, fixed( "Tinf_RespPK" ), paste0( "Tinf_", target ) )
}

#' Remap one PK/PD library equation string onto project compartments.
#' @param text Character equation string.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Updated character string.
#' @keywords internal
remapPkpdLibraryText = function( text, pfimproject ) {
  comp = .pkpdLibraryCompartments( pfimproject )
  text = reduce2( names( comp ), unname( comp ), replaceVariablesLibraryOfModels, .init = text )
  .remapLibraryDoseTokens( text, pfimproject )
}

#' Remap a PK/PD equation list onto project compartments.
#' @param equations Character vector or during/after infusion lists.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Remapped equations.
#' @keywords internal
remapPkpdLibraryEquations = function( equations, pfimproject ) {
  .remapLibraryEquations( equations, \(text) remapPkpdLibraryText( text, pfimproject ) )
}

#' Build \code{Deriv_<compartment>} names for the first \code{nEq} compartments.
#' @param pfimproject A \code{PFIMProject} object.
#' @param nEq Integer number of derivative equations.
#' @return Character vector of derivative names.
#' @noRd
#' @keywords internal
.derivativeNamesFromCompartments = function( pfimproject, nEq, catalogNames = NULL ) {
  nEq = as.integer( nEq )[ 1L ]
  if ( !length( nEq ) || is.na( nEq ) || nEq < 1L )
    return( character( 0 ) )
  mapped = unname( .pkpdLibraryCompartments( pfimproject ) )
  mapped = mapped[ !is.na( mapped ) & nzchar( mapped ) ]
  catalogStates = character( 0 )
  if ( length( catalogNames ) ) {
    catalogStates = sub( "^Deriv_", "", as.character( catalogNames ) )
    catalogStates = catalogStates[ !is.na( catalogStates ) & nzchar( catalogStates ) ]
  }
  states = catalogStates
  if ( !length( states ) ) {
    states = mapped
  } else if ( length( mapped ) ) {
    nMap = min( length( mapped ), length( states ) )
    states[ seq_len( nMap ) ] = mapped[ seq_len( nMap ) ]
    if ( length( mapped ) > length( states ) )
      states = c( states, mapped[ seq( length( states ) + 1L, length( mapped ) ) ] )
  }
  if ( length( states ) < nEq )
    states = c( states, paste0( "C", seq( length( states ) + 1L, nEq ) ) )
  states = states[ seq_len( nEq ) ]
  if ( any( is.na( states ) | !nzchar( states ) ) )
    .pfimStop(
      "could not name ", nEq, " ODE compartments from outputs or library equations."
    )
  paste0( "Deriv_", states )
}

#' Map library ODE PK tokens \code{C1}/\code{C2} onto project PK/PD compartments.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Named character vector (\code{C1}, optionally \code{C2}).
#' @noRd
#' @keywords internal
.pkpdOdePkCompartments = function( pfimproject ) {
  comp = .pkpdLibraryCompartments( pfimproject )
  out = c( C1 = unname( comp["RespPK"] ) )
  if ( length( comp ) >= 2L ) out = c( out, C2 = unname( comp["E"] ) )
  out
}

#' Remap one ODE PK library equation string (\code{C1}/\code{C2} tokens).
#' @param text Character equation string.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Updated character string.
#' @keywords internal
remapOdePkLibraryText = function( text, pfimproject ) {
  comp = .pkpdOdePkCompartments( pfimproject )
  text = reduce2( names( comp ), unname( comp ), replaceVariablesLibraryOfModels, .init = text )
  .remapLibraryDoseTokens( text, pfimproject )
}

#' Remap an ODE PK equation list onto project compartments.
#' @param equations Character vector or during/after infusion lists.
#' @param pfimproject A \code{PFIMProject} object.
#' @return Remapped equations.
#' @keywords internal
remapOdePkLibraryEquations = function( equations, pfimproject ) {
  .remapLibraryEquations( equations, \(text) remapOdePkLibraryText( text, pfimproject ) )
}

#' Concatenate remapped PK and PD equations under \code{Deriv_*} names.
#' @param pkPart Character/list PK equations for one infusion phase.
#' @param pdEq Character/list PD equations.
#' @param derivNames Character vector of combined derivative names.
#' @return Named list of combined equations.
#' @noRd
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
