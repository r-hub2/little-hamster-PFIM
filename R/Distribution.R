#' @title Distribution
#' @description
#' Base class for parameter distributions (\code{Normal}, \code{LogNormal}).
#' @param name Character string: distribution label.
#' @param mu Typical value on the natural (untransformed) scale. For
#'   \code{LogNormal}, IIV is additive on \eqn{\log\theta} with
#'   \eqn{\theta=\mu\,e^{\eta}} at \eqn{\eta=0}.
#' @param omega Standard deviation of random effects (IIV).
#' @return An S7 object of class \code{Distribution}.
#' @export
#'
Distribution = new_class(
  "Distribution",
  package = "PFIM",
  properties = list(
    name = new_property(class_character, default = character(0)),
    mu = new_property(class_double, default = 0.0),
    omega = new_property(class_double, default = 0.0)
  ),
  validator = function( self ) {
    om = prop( self, "omega" )
    if ( length( om ) && ( any( is.na( om ) ) || any( om < 0 ) ) )
      return( "Distribution: omega must be non-negative." )
    mu = prop( self, "mu" )
    if ( length( mu ) && any( !is.finite( mu ) ) )
      return( "Distribution: mu must be finite." )
    NULL
  })


#' Scale gradients by the distribution link (identity or log-normal).
#'
#' Used when mapping \eqn{\partial f / \partial \theta} to
#' \eqn{\partial f / \partial \eta} at the typical value for the population FIM.
#' @param distribution A \code{Distribution} object.
#' @param gradient Numeric gradient vector before distribution adjustment.
#' @param thetaValue Numeric parameter value on the natural scale.
#' @name adjustGradient
#' @keywords internal
adjustGradient = new_generic( "adjustGradient", c( "distribution" ),
                              function( distribution, gradient, thetaValue ) {
                                S7_dispatch()
                              })

#' Fallback for base \code{Distribution} (treat as log-normal link).
#'
#' Multiplies the gradient by \code{thetaValue} (same as \code{LogNormal}).
#' @param distribution First argument of generic.
#' @param gradient Numeric gradient vector before distribution adjustment.
#' @param thetaValue Numeric parameter value on the natural scale.
#' @return Gradient scaled by \code{thetaValue}.
#' @name adjustGradient
#' @keywords internal
method( adjustGradient, Distribution ) = function( distribution, gradient, thetaValue ) {
  gradient * thetaValue
}
