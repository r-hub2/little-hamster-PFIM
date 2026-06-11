#' @title Constant
#' @description
#' Pure additive residual error: constant variance \code{sigmaInter^2}.
#' @inheritParams ModelError
#' @include ModelError.R
#' @export

Constant = new_class("Constant", package = "PFIM", parent = ModelError,
                     properties = list(
                       output = new_property(class_character, default = character(0)),
                       equation = new_property(class_expression, default = expression(sigmaInter)),
                       derivatives = new_property(class_list, default = list()),
                       sigmaInter = new_property(class_double, default = 0.0),
                       sigmaSlope = new_property(class_double, default = 0.0),
                       sigmaInterFixed = new_property(class_logical, default = FALSE),
                       sigmaSlopeFixed = new_property(class_logical, default = FALSE),
                       cError = new_property(class_double, default = 1.0)
                     ),
                     constructor = function(output = character(0),
                                            equation = expression(sigmaInter),
                                            derivatives = list(),
                                            sigmaInter = 0.0,
                                            sigmaSlope = 0.0,
                                            sigmaInterFixed = FALSE,
                                            sigmaSlopeFixed = FALSE,
                                            cError = 1.0) {
                       new_object(.parent = ModelError,
                                  output = output,
                                  equation = equation,
                                  derivatives = derivatives,
                                  sigmaInter = sigmaInter,
                                  sigmaSlope = sigmaSlope,
                                  sigmaInterFixed = sigmaInterFixed,
                                  sigmaSlopeFixed = sigmaSlopeFixed,
                                  cError = cError)
                     })






