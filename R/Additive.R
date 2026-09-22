#' @title Additive
#' @description
#' Additive covariate link: theta(cov) = mu * (1 + beta * cov).
#' @inheritParams CovariateModelEquation
#' @include CovariateModelEquation.R
#' @return An S7 object of class \code{Additive}.
#' @export
Additive = new_class(
  "Additive", package = "PFIM", parent = CovariateModelEquation,
  constructor = function( beta = 0.0, combinedEffect = 0.0, value = 0.0 ) {
    # `value` is derived; the argument is kept for a stable constructor signature.
    new_object(
      CovariateModelEquation(
        beta           = beta,
        combinedEffect = combinedEffect,
        value          = beta * ( 1 + combinedEffect )
      )
    )
  }
)

#' Apply the additive covariate link to typical values.
#'
#' Builds a new \code{Additive} whose \code{value} is
#' \code{beta * (1 + combinedEffect)} (elementwise over parameters).
#' @param equation First argument of generic.
#' @param beta Named numeric vector of baseline parameter values (mu).
#' @param combinedEffect Named numeric vector of summed covariate effects.
#' @return \code{Additive} object containing transformed parameter values.
#' @name computeCovariateValue
#' @keywords internal
method( computeCovariateValue, Additive ) = function( equation, beta, combinedEffect ) {
  Additive( beta = beta, combinedEffect = combinedEffect )
}
