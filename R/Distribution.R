#' @title Distribution
#' @description
#' Base class for parameter distributions (\code{Normal}, \code{LogNormal}).
#' @param name Character string: distribution label.
#' @param mu Mean on the natural scale (or log-scale mean for log-normal).
#' @param omega Standard deviation of random effects (IIV).
#' @export
#'
Distribution = new_class(
  "Distribution",
  package = "PFIM",
  properties = list(
    name = new_property(class_character, default = character(0)),
    mu = new_property(class_double, default = 0.0),
    omega = new_property(class_double, default = 0.0)
  ))


adjustGradient = new_generic( "adjustGradient", c( "distribution" ),
                              function( distribution, gradient, thetaValue ) {
                                S7_dispatch()
                              })
