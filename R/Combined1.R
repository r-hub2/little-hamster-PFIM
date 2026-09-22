#' @title Combined1
#' @description
#' Combined residual error (additive + proportional SD):
#' variance \code{(sigmaInter + sigmaSlope * f^{cError})^2}.
#' @inheritParams ModelError
#' @include ModelError.R
#' @return An S7 object of class \code{Combined1}.
#' @export
Combined1 = new_class(
  "Combined1", package = "PFIM", parent = ModelError,
  # Lock form so mutation cannot silently switch to Combined2 math.
  validator = function( self ) {
    if ( !identical( prop( self, "varianceForm" ), "combined1" ) )
      return( "Combined1: varianceForm must be 'combined1'." )
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
      cError, "combined1", ...
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
