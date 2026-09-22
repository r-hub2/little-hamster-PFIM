#' @title Exponential
#' @description
#' Exponential covariate link: theta(cov) = mu * exp(beta * cov).
#' @inheritParams CovariateModelEquation
#' @include CovariateModelEquation.R
#' @return An S7 object of class \code{Exponential}.
#' @export
Exponential = new_class(
  "Exponential", package = "PFIM", parent = CovariateModelEquation,
  constructor = function( beta = 0.0, combinedEffect = 0.0, value = 0.0 ) {
    # `value` is derived; the argument is kept for a stable constructor signature.
    new_object(
      CovariateModelEquation(
        beta           = beta,
        combinedEffect = combinedEffect,
        value          = beta * exp( combinedEffect )
      )
    )
  }
)

#' Apply the exponential covariate link to typical values.
#'
#' Builds a new \code{Exponential} whose \code{value} is
#' \code{beta * exp(combinedEffect)} (elementwise over parameters).
#' @param equation First argument of generic.
#' @param beta Named numeric vector of baseline parameter values (mu).
#' @param combinedEffect Named numeric vector of summed covariate effects.
#' @return \code{Exponential} object containing transformed parameter values.
#' @name computeCovariateValue
#' @keywords internal
method( computeCovariateValue, Exponential ) = function( equation, beta, combinedEffect ) {
  Exponential( beta = beta, combinedEffect = combinedEffect )
}
