#' @title ModelError
#' @description
#' Base class for residual error models (\code{Constant}, \code{Combined1}, \code{Proportional}).
#' @param output Character string: outcome name.
#' @param equation Expression defining residual variance as a function of predictions.
#' @param derivatives List of partial derivatives w.r.t. error parameters.
#' @param sigmaInter Additive error standard deviation (inter-subject / residual).
#' @param sigmaSlope Proportional error standard deviation (slope).
#' @param sigmaInterFixed If TRUE, \code{sigmaInter} is fixed in the FIM.
#' @param sigmaSlopeFixed If TRUE, \code{sigmaSlope} is fixed in the FIM.
#' @param cError Power parameter in the error model (when applicable).
#' @export

ModelError = new_class("ModelError", package = "PFIM",
                       properties = list(
                         output          = new_property(class_character,  default = "output"),
                         equation        = new_property(class_expression, default = expression()),
                         derivatives     = new_property(class_list,       default = list()),
                         sigmaInter      = new_property(class_double,     default = 0.0),
                         sigmaSlope      = new_property(class_double,     default = 0.0),
                         sigmaInterFixed = new_property(class_logical,    default = FALSE),
                         sigmaSlopeFixed = new_property(class_logical,    default = FALSE),
                         cError          = new_property(class_double,     default = 1.0)
                       ),
                       constructor = function(output          = "output",
                                              equation        = expression(),
                                              derivatives     = list(),
                                              sigmaInter      = 0.0,
                                              sigmaSlope      = 0.0,
                                              sigmaInterFixed = FALSE,
                                              sigmaSlopeFixed = FALSE,
                                              cError          = 1.0) {
                         new_object(.parent         = ModelError,
                                    output          = output,
                                    equation        = equation,
                                    derivatives     = derivatives,
                                    sigmaInter      = sigmaInter,
                                    sigmaSlope      = sigmaSlope,
                                    sigmaInterFixed = sigmaInterFixed,
                                    sigmaSlopeFixed = sigmaSlopeFixed,
                                    cError          = cError)
                       })

evaluateErrorModelDerivatives = new_generic("evaluateErrorModelDerivatives", c("modelError"))
getModelErrorData              = new_generic("getModelErrorData",              c("modelError"))

#' Derivatives of the residual error model
#' @name evaluateErrorModelDerivatives
#' @export

method( evaluateErrorModelDerivatives, ModelError ) = function( modelError, evaluationModel ) {

  sigmaInter      = prop( modelError, "sigmaInter" )
  sigmaSlope      = prop( modelError, "sigmaSlope" )
  sigmaInterFixed = prop( modelError, "sigmaInterFixed" )
  sigmaSlopeFixed = prop( modelError, "sigmaSlopeFixed" )

  n = length( evaluationModel )
  I = diag( n )

  varianceExpr = expression((sigmaInter + sigmaSlope * evaluationModel)^2)

  # Only estimable (non-zero, non-fixed) sigma parameters
  sigmaSpecs = list(
    list(name = "sigmaInter", value = sigmaInter, fixed = sigmaInterFixed),
    list(name = "sigmaSlope", value = sigmaSlope, fixed = sigmaSlopeFixed)
  ) |> keep(~ .x$value != 0 && !.x$fixed)

  sigmaDerivatives = sigmaSpecs |>
    map(function(spec) {
      derivMat = eval(D(varianceExpr, spec$name)) * I
      setNames(list(derivMat), spec$name)
    }) |>
    list_flatten()

  list(
    sigmaDerivatives = sigmaDerivatives,
    errorVariance    = (sigmaInter + sigmaSlope * evaluationModel)^2 * I
  )
}

#' Residual error settings for reports
#' @name getModelErrorData
#' @export

method( getModelErrorData, ModelError ) = function( modelError ) {
  data.frame(
    output     = prop( modelError, "output"     ),
    type       = str_remove( class( modelError )[[1L]], "^PFIM::" ),
    sigmaSlope = as.character( prop( modelError, "sigmaSlope" ) ),
    sigmaInter = as.character( prop( modelError, "sigmaInter" ) ),
    stringsAsFactors = FALSE
  )
}
