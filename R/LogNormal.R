#' @title LogNormal
#' @description Log-normal distribution for positive PK/PD parameters (default in PFIM).
#' @inheritParams Distribution
#' @include Distribution.R
#' @return An S7 object of class \code{LogNormal}.
#' @export

LogNormal = new_class( "LogNormal", package = "PFIM", parent = Distribution )

#' Log-normal link gradient adjustment.
#'
#' Scales by \code{thetaValue} so \eqn{\partial f / \partial \eta =
#' \theta \cdot \partial f / \partial \theta} at the typical value.
#' @param distribution First argument of generic.
#' @param gradient Numeric gradient vector before distribution adjustment.
#' @param thetaValue Numeric parameter value on the natural scale.
#' @return Gradient scaled by \code{thetaValue}.
#' @name adjustGradient
#' @keywords internal
method( adjustGradient, LogNormal ) = function( distribution, gradient, thetaValue ) {
  return( gradient * thetaValue )
}
