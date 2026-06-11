#' @description
#' The class \code{CovariateModelEquation} is the base class for the two
#' covariate structural equations supported by PFIM:
#' \code{Additive} and \code{Exponential}.
#'
#' A covariate model equation adjusts the population mean for individual i
#' with covariate value
#' \itemize{
#'   \item \code{Additive}:    additive model
#'   \item \code{Exponential}: exponential model
#' }
#'
#' @title CovariateModelEquation
#' @param beta           Fixed-effect coefficient Î² (default 0).
#' @param combinedEffect Combined effect term Î£ Î²â‚–Â·covâ‚– (default 0).
#' @param value          Computed adjusted value (default 0).
#' @export

CovariateModelEquation = new_class( "CovariateModelEquation", package = "PFIM",
                                    properties = list(
                                      beta           = new_property( class_double, default = 0.0 ),
                                      combinedEffect = new_property( class_double, default = 0.0 ),
                                      value          = new_property( class_double, default = 0.0 )
                                    ),
                                    constructor = function( beta           = 0.0,
                                                            combinedEffect = 0.0,
                                                            value          = 0.0 )
                                      new_object( CovariateModelEquation,
                                                  beta           = beta,
                                                  combinedEffect = combinedEffect,
                                                  value          = value )
)

#' computeCovariateValue: apply beta + combinedEffect to produce the adjusted value.
#' @name computeCovariateValue
#' @usage NULL
#' @param equation       A \code{CovariateModelEquation} subclass object.
#' @param beta           Numeric scalar -- fixed-effect coefficient.
#' @param combinedEffect Numeric scalar -- combined covariate term.
#' @return An updated \code{CovariateModelEquation} subclass object.
#' @export

computeCovariateValue = new_generic( "computeCovariateValue", "equation" )
