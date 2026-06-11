#' @title Normal
#' @description Normal distribution for model parameters (identity link on the mean).
#' @inheritParams Distribution
#' @include Distribution.R
#' @export

Normal = new_class( "Normal", package = "PFIM", parent = Distribution )

method( adjustGradient, Normal ) = function( distribution, gradient, thetaValue ) {
  return( gradient )
}
