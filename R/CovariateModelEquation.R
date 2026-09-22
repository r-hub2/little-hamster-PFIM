#' @description
#' The class \code{CovariateModelEquation} is the base class for the two
#' covariate structural equations supported by PFIM:
#' \code{Additive} and \code{Exponential}.
#'
#' A covariate model equation adjusts the population mean for individual i
#' with covariate value
#' \itemize{
#'   \item \code{Additive}:    \eqn{\theta = \mu \cdot (1 + \beta \cdot \mathrm{cov})}
#'   \item \code{Exponential}: \eqn{\theta = \mu \cdot \exp(\beta \cdot \mathrm{cov})}
#' }
#' Concrete subclasses implement \code{computeCovariateValue()}.
#'
#' @title CovariateModelEquation
#' @param beta           Fixed-effect coefficient (default 0).
#' @param combinedEffect Combined covariate effect term (default 0).
#' @param value          Computed adjusted value (default 0).
#' @return An S7 object of class \code{CovariateModelEquation}.
#' @export

CovariateModelEquation = new_class(
  "CovariateModelEquation", package = "PFIM",
  properties = list(
    beta           = new_property( class_double, default = 0.0 ),
    combinedEffect = new_property( class_double, default = 0.0 ),
    value          = new_property( class_double, default = 0.0 )
  ),
  constructor = function( beta = 0.0, combinedEffect = 0.0, value = 0.0 ) {
    new_object(
      S7_object(),
      beta           = beta,
      combinedEffect = combinedEffect,
      value          = value
    )
  }
)

#' Apply \code{beta} and \code{combinedEffect} to produce the adjusted parameter value.
#'
#' Dispatches to \code{Additive} or \code{Exponential} methods (see those classes).
#' @name computeCovariateValue
#' @usage NULL
#' @param equation       A \code{CovariateModelEquation} subclass object.
#' @param beta           Numeric scalar -- fixed-effect coefficient.
#' @param combinedEffect Numeric scalar -- combined covariate term.
#' @return An updated \code{CovariateModelEquation} subclass object.
#' @keywords internal

computeCovariateValue = new_generic( "computeCovariateValue", "equation" )
