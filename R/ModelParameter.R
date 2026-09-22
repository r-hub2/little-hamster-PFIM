#' @title ModelParameter
#' @description
#' One population model parameter: distribution (mu, omega), optional IOV (gamma), and fix flags.
#' @param name          Character string: name of the parameter.
#' @param gamma         Numeric: the SD for inter-occasion variability (IOV).
#'                      \code{gamma^2} is the IOV variance component; \code{0} means no IOV.
#' @param distribution  An object of class \code{Distribution} for this parameter.
#' @param fixedMu       Logical: TRUE if \code{mu} is fixed (not estimated).
#' @param fixedOmega    Logical: TRUE if \code{omega} is fixed (not estimated).
#' @param value         Reserved (unused); use \code{prop(distribution, "mu")} for reports and FIM.
#' @include Distribution.R
#' @include LogNormal.R
#' @return An S7 object of class \code{ModelParameter}.
#' @export

ModelParameter = new_class( "ModelParameter",
                            package = "PFIM",
                            properties = list(
                              name         = new_property(class_character, default = character(0)),
                              gamma        = new_property(class_double,    default = 0),
                              distribution = new_property(NULL | Distribution, default = NULL),
                              fixedMu      = new_property(class_logical,   default = FALSE),
                              fixedOmega   = new_property(class_logical,   default = FALSE),
                              value        = new_property(class_double,    default = numeric(0))
                            ),
                            validator = function( self ) {
                              n = prop( self, "name" )
                              if ( length( n ) > 1L )
                                return( "ModelParameter: name must be a single string." )
                              if ( .pfimIsBlankScalar( n ) )
                                return( "ModelParameter: name must be a non-empty string." )
                              g = prop( self, "gamma" )
                              if ( length( g ) && ( is.na( g ) || g < 0 ) )
                                return( "ModelParameter: gamma must be non-negative." )
                              d = prop( self, "distribution" )
                              if ( !is.null( d ) ) {
                                om = prop( d, "omega" )
                                if ( length( om ) && ( is.na( om ) || om < 0 ) )
                                  return( "ModelParameter: omega must be non-negative." )
                              }
                              NULL
                            })

#' Read one field from a parameter's distribution slot.
#' @param parameter A \code{ModelParameter} object.
#' @param field Character field name (\code{"mu"} or \code{"omega"}).
#' @return Value of the requested distribution field.
#' @noRd
#' @keywords internal
.paramDist = function( parameter, field ) {
  prop( prop( parameter, "distribution" ), field )
}

#' Whether \code{mu} is estimable (unfixed).
#'
#' \code{Normal(mu = 0)} stays estimable. \code{LogNormal(mu = 0)} is degenerate
#' (\eqn{\theta = 0}) and is treated as non-estimable with a one-shot warning.
#' @param parameter A \code{ModelParameter} object.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.paramMuEstimable = function( parameter ) {
  if ( isTRUE( prop( parameter, "fixedMu" ) ) ) return( FALSE )
  d  = prop( parameter, "distribution" )
  if ( is.null( d ) ) return( FALSE )
  mu = prop( d, "mu" )
  if ( S7::S7_inherits( d, LogNormal ) && length( mu ) && isTRUE( mu == 0 ) ) {
    .pfimWarnOnce(
      "param.mu.zero.lognormal",
      "LogNormal(mu = 0) for parameter '", prop( parameter, "name" ),
      "' is non-estimable (degenerate typical value).",
      id = prop( parameter, "name" )
    )
    return( FALSE )
  }
  TRUE
}

#' Whether \code{mu} is treated as fixed in the FIM.
#' @param parameter A \code{ModelParameter} object.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.paramMuFixed = function( parameter ) {
  !.paramMuEstimable( parameter )
}

#' Whether a parameter carries IIV (\eqn{\omega > 0}), ignoring fix flags.
#'
#' Population FIM drops \eqn{\mu} via \code{fixedMu} and \eqn{\omega^2} via
#' \code{fixedOmega} independently. Bayesian MAP still has
#' \eqn{\eta\sim N(0,\omega^2)} whenever \eqn{\omega>0}.
#' @noRd
#' @keywords internal
.paramHasIiv = function( parameter ) {
  .paramDist( parameter, "omega" ) > 0
}

#' Whether a parameter carries IOV (\eqn{\gamma > 0}), ignoring fix flags.
#' @noRd
#' @keywords internal
.paramHasIov = function( parameter ) {
  pluck( parameter, "gamma", .default = 0 ) > 0
}

#' Bayesian eta-parameter list: every parameter with \eqn{\omega > 0}.
#' @noRd
#' @keywords internal
.pfimBayesianEtaParameters = function( parameters ) {
  parameters[ map_lgl( parameters, .paramHasIiv ) ]
}

#' Whether \code{omega} is estimable (\eqn{\omega > 0} and not \code{fixedOmega}).
#'
#' Independent of \code{fixedMu}: a typical value may be fixed while its IIV
#' is still estimated.
#' @param parameter A \code{ModelParameter} object.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.paramOmegaEstimable = function( parameter ) {
  !isTRUE( prop( parameter, "fixedOmega" ) ) &&
    .paramHasIiv( parameter )
}

#' Whether \code{omega} is treated as fixed in the FIM.
#' @param parameter A \code{ModelParameter} object.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.paramOmegaFixed = function( parameter ) {
  !.paramOmegaEstimable( parameter )
}

#' Whether \code{gamma} (IOV SD) is estimable.
#'
#' Requires \code{gamma > 0} and unfixed omega (IOV shares the omega fix flag).
#' @param parameter A \code{ModelParameter} object.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.paramGammaEstimable = function( parameter ) {
  pluck( parameter, "gamma", .default = 0 ) > 0 &&
    !isTRUE( prop( parameter, "fixedOmega" ) )
}

#' Whether any parameter carries positive IOV (\code{gamma > 0}).
#' @param parameters List of \code{ModelParameter} objects.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.hasIovParameters = function( parameters ) {
  any( map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) ) > 0 )
}

#' Names of parameters matching \code{predicate}, with optional prefix.
#' @param parameters List of \code{ModelParameter} objects.
#' @param predicate Function returning logical for one parameter.
#' @param prefix Character prefix (e.g. \code{"mu_"}).
#' @return Character vector of labelled names.
#' @noRd
#' @keywords internal
.paramLabelNames = function( parameters, predicate, prefix = "" ) {
  parameters |>
    keep( predicate ) |>
    map_chr( ~ paste0( prefix, prop( .x, "name" ) ) )
}

#' Names of estimable \code{mu} parameters (\code{prefix} optional: \code{"mu_"}, greek, etc.).
#' @param parameters List of \code{ModelParameter} objects.
#' @param prefix Character prefix prepended to each name.
#' @return Character vector of estimable mu labels.
#' @noRd
#' @keywords internal
.estimableMuNames = function( parameters, prefix = "" ) {
  .paramLabelNames( parameters, .paramMuEstimable, prefix )
}

#' Typical values (\code{mu}) for estimable fixed effects.
#' @param parameters List of \code{ModelParameter} objects.
#' @return Numeric vector of mu values.
#' @noRd
#' @keywords internal
.paramMuValues = function( parameters ) {
  parameters |>
    keep( .paramMuEstimable ) |>
    map_dbl( ~ .paramDist( .x, "mu" ) )
}

#' IIV variances (\code{omega^2}) for estimable random effects.
#' @param parameters List of \code{ModelParameter} objects.
#' @return Numeric vector of omega-squared values.
#' @noRd
#' @keywords internal
.paramOmegaSqValues = function( parameters ) {
  parameters |>
    keep( .paramOmegaEstimable ) |>
    map_dbl( ~ .paramDist( .x, "omega" )^2 )
}

#' IOV variances (\code{gamma^2}) for estimable occasion effects.
#' @param parameters List of \code{ModelParameter} objects.
#' @return Numeric vector of gamma-squared values.
#' @noRd
#' @keywords internal
.paramGammaSqValues = function( parameters ) {
  parameters |>
    keep( .paramGammaEstimable ) |>
    map_dbl( ~ pluck( .x, "gamma" )^2 )
}

getModelParametersData = new_generic( "getModelParametersData", c( "modelParameter" ) )

#' Parameter table for reports (mu, omega2, gamma2, fix flags, distribution).
#' @param modelParameter A \code{ModelParameter} object.
#' @return A data frame of model parameter values.
#' @name getModelParametersData
#' @export
method( getModelParametersData, ModelParameter ) = function( modelParameter ) {

  list( modelParameter ) |>
    map( function( parameter ) {
      dist = prop( parameter, "distribution" )
      list(
        Parameters   = prop( parameter, "name"      ),
        mu           = as.character( prop( dist, "mu"    ) ),
        omega2       = as.character( prop( dist, "omega" )^2 ),
        gamma2       = as.character( prop( parameter, "gamma" )^2 ),
        Distribution = str_remove( class( dist )[1], "PFIM::" ),
        mu_fixed     = as.character( .paramMuFixed( parameter ) ),
        omega2_fixed = as.character( .paramOmegaFixed( parameter ) )
      )
    }) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
}
