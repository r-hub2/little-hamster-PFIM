#' @title SamplingTimes
#' @description
#' Observation times for one outcome within an arm.
#'
#' Times are in the model's time unit (often hours). Multiple
#' \code{SamplingTimes} objects on an arm cover multi-response designs.
#' @param outcome Character string: outcome name.
#' @param samplings Numeric vector: sampling times (hours or model time unit).
#' @return An S7 object of class \code{SamplingTimes}.
#' @export

SamplingTimes = new_class("SamplingTimes", package = "PFIM",

                          properties = list(
                            outcome = new_property(class_character, default = character(0)),
                            samplings = new_property(class_double, default = numeric(0))
                          ),
                          validator = function( self ) {
                            s = prop( self, "samplings" )
                            if ( length( s ) && any( !is.finite( s ) ) )
                              return( "SamplingTimes: samplings must be finite." )
                            o = prop( self, "outcome" )
                            if ( length( o ) > 1L )
                              return( "SamplingTimes: outcome must be a single string." )
                            if ( .pfimIsBlankScalar( o ) )
                              return( "SamplingTimes: outcome must be a non-empty string." )
                            NULL
                          })
