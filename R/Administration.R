#' Validate that every element of a list inherits from an S7 class.
#' @param items List to check (empty lists are OK).
#' @param class Target S7 class object.
#' @param label Property path used in error messages.
#' @param typeName Human-readable class name for the message.
#' @return \code{NULL} if valid, otherwise an error string for the S7 validator.
#' @noRd
#' @keywords internal
.validateS7List = function( items, class, label, typeName ) {
  if ( !length( items ) )
    return( NULL )
  bad = which( !vapply( items, S7::S7_inherits, logical( 1L ), class ) )
  if ( length( bad ) )
    return( sprintf( "%s[[%d]] must be a %s object.", label, bad[ 1L ], typeName ) )
  NULL
}

#' @title Administration
#' @description
#' Dosing regimen for one outcome: dose amounts, times, infusion duration, and
#' optional inter-dose interval (\code{tau}) for repeated dosing or steady state.
#'
#' When \code{timeDose} has length 1 and \code{dose} is longer, solvers treat that
#' as a shared administration time (see \code{.alignAdministrationDosing}).
#' @param outcome Character string: model output receiving the dose.
#'   Library models typically use the catalogue name (typically \code{"RespPK"}).
#'   An \code{outputs} alias is accepted when it maps to a declared state; the
#'   dose is applied to that state. User-written ODEs may use the compartment
#'   name (\code{"Cc"}) or that alias. Sampling times may use either the state
#'   or the alias.
#' @param timeDose Numeric vector: administration times.
#' @param dose Numeric vector: dose amounts (same length as \code{timeDose}, or
#'   longer when \code{timeDose} has length 1 - replicated for each dose).
#' @param Tinf Numeric vector: infusion duration (0 or omitted for bolus).
#' @param tau Numeric: dosing interval for repeated doses or steady-state models.
#' @return An S7 object of class \code{Administration}.
#' @export

Administration = new_class("Administration",
                           package = "PFIM",
                           properties = list(
                             outcome = new_property(class_character, default = character(0)),
                             timeDose = new_property(class_double, default = numeric(0)),
                             dose = new_property(class_double, default = numeric(0)),
                             Tinf = new_property(class_double, default = numeric(0)),
                             tau = new_property(class_double, default = 0.0)
                          ),
                           validator = function( self ) {
                             o = prop( self, "outcome" )
                             if ( length( o ) > 1L )
                               return( "Administration: outcome must be a single string." )
                             if ( .pfimIsBlankScalar( o ) )
                               return( "Administration: outcome must be a non-empty string." )
                             td = prop( self, "timeDose" )
                             d  = prop( self, "dose" )
                             if ( length( td ) && any( !is.finite( td ) ) )
                               return( "Administration: timeDose must be finite." )
                             if ( length( d ) && any( !is.finite( d ) ) )
                               return( "Administration: dose must be finite." )
                             # Allow length(timeDose)==1 with multiple doses (legacy scripts).
                             if ( length( td ) > 0L && length( d ) > 0L && length( td ) != length( d ) ) {
                               if ( length( td ) != 1L )
                                 return( "Administration: length(timeDose) must equal length(dose), unless timeDose has length 1." )
                             }
                             tau = prop( self, "tau" )
                             if ( length( tau ) != 1L || !is.finite( tau ) || tau < 0 )
                               return( "Administration: tau must be a non-negative finite number." )
                             tinf = prop( self, "Tinf" )
                             if ( length( tinf ) && any( !is.finite( tinf ) | tinf < 0 ) )
                               return( "Administration: Tinf must be finite and non-negative." )
                             NULL
                           })

#' Align administration dosing vectors for downstream solvers.
#'
#' Legacy PFIM scripts may specify a single \code{timeDose} with multiple
#' \code{dose} values (e.g. \code{timeDose = 0}, \code{dose = c(200, 100)}).
#' Replicates \code{timeDose} and \code{Tinf} so all three vectors share length.
#' @param administration An \code{Administration} object.
#' @return List with aligned \code{timeDose}, \code{dose}, and \code{Tinf}.
#' @noRd
#' @keywords internal
.alignAdministrationDosing = function( administration ) {
  timeDose = prop( administration, "timeDose" )
  dose     = prop( administration, "dose" )
  Tinf     = prop( administration, "Tinf" )

  if ( length( timeDose ) == 1L && length( dose ) > 1L )
    timeDose = rep( timeDose, length( dose ) )
  if ( length( Tinf ) == 1L && length( dose ) > 1L && length( timeDose ) == length( dose ) )
    Tinf = rep( Tinf, length( dose ) )

  list( timeDose = timeDose, dose = dose, Tinf = Tinf )
}
