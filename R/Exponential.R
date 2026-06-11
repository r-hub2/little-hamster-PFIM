#' @title Exponential
#' @description
#' Exponential covariate link: theta(cov) = mu * exp(beta * cov).
#' @inheritParams CovariateModelEquation
#' @include CovariateModelEquation.R
#' @export

Exponential = new_class("Exponential", package = "PFIM", parent = CovariateModelEquation,
                     properties = list(
                       beta = new_property( class_double, default = 0 ),
                       combinedEffect = new_property( class_double, default = 0 ),
                       value = new_property( class_double, default = 0 )
                     ),
                     constructor = function( beta = 0.0,
                                            combinedEffect = 0.0,
                                            value = 0.0) {
                       new_object( .parent = CovariateModelEquation,
                                   beta = beta,
                                   combinedEffect = combinedEffect,
                                   value = beta * exp( combinedEffect ) )
                     })

method(computeCovariateValue, Exponential) = function(equation, beta, combinedEffect) {
  Exponential(beta = beta, combinedEffect = combinedEffect)
}
