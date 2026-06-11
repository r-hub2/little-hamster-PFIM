#' @title SamplingTimes
#' @description
#' Observation times for one outcome within an arm.
#' @param outcome Character string: outcome name.
#' @param samplings Numeric vector: sampling times (hours or model time unit).
#' @export

SamplingTimes = new_class("SamplingTimes", package = "PFIM",

                          properties = list(
                            outcome = new_property(class_character, default = character(0)),
                            samplings = new_property(class_double, default = numeric(0))
                          ))
