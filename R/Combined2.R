#' @title Combined2
#' @description
#' Combined residual error with independent additive and proportional SDs:
#' variance \code{sigmaInter^2 + (sigmaSlope * f^{cError})^2}.
#' Matches PopED \code{y*(1+e_prop)+e_add} when both sigmas are residual SDs.
#' @inheritParams ModelError
#' @include ModelError.R
#' @return An S7 object of class \code{Combined2}.
#' @export
Combined2 = new_class(
  "Combined2", package = "PFIM", parent = ModelError,
  validator = function( self ) {
    if ( !identical( prop( self, "varianceForm" ), "combined2" ) )
      return( "Combined2: varianceForm must be 'combined2'." )
    NULL
  },
  constructor = function( output          = character( 0 ),
                          sigmaInter      = 0.0,
                          sigmaSlope      = 0.0,
                          sigmaInterFixed = FALSE,
                          sigmaSlopeFixed = FALSE,
                          cError          = 1.0,
                          ... ) {
    a = .pfimConstructModelErrorArgs(
      output, sigmaInter, sigmaSlope, sigmaInterFixed, sigmaSlopeFixed,
      cError, "combined2", ...
    )
    new_object(
      ModelError(
        output          = a$output,
        sigmaInter      = a$sigmaInter,
        sigmaSlope      = a$sigmaSlope,
        sigmaInterFixed = a$sigmaInterFixed,
        sigmaSlopeFixed = a$sigmaSlopeFixed,
        cError          = a$cError,
        varianceForm    = a$varianceForm
      )
    )
  }
)
