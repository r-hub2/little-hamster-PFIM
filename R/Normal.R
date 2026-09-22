#' @title Normal
#' @description Normal distribution for model parameters (identity link on the mean).
#' @inheritParams Distribution
#' @include Distribution.R
#' @return An S7 object of class \code{Normal}.
#' @export

Normal = new_class( "Normal", package = "PFIM", parent = Distribution )

#' Identity-link gradient adjustment (no scaling).
#'
#' For a normal random-effect model, \eqn{\partial f / \partial \eta = \partial f / \partial \theta}.
#' @param distribution First argument of generic.
#' @param gradient Numeric gradient vector before distribution adjustment.
#' @param thetaValue Numeric parameter value (unused for normal link).
#' @return Unchanged numeric gradient vector.
#' @name adjustGradient
#' @keywords internal
method( adjustGradient, Normal ) = function( distribution, gradient, thetaValue ) {
  return( gradient )
}
