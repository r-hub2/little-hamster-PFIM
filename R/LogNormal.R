#' @title LogNormal
#' @description Log-normal distribution for positive PK/PD parameters (default in PFIM).
#' @inheritParams Distribution
#' @include Distribution.R
#' @export

LogNormal = new_class( "LogNormal", package = "PFIM", parent = Distribution )

#' Adjust gradients for fixed parameters
#' @name adjustGradient
#' @export

method( adjustGradient, LogNormal ) = function( distribution, gradient, thetaValue ) {
  return( gradient * thetaValue )
}
