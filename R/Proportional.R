#' @title Proportional
#' @description
#' Pure proportional residual error: \code{V = (sigmaSlope * f^{cError})^2}
#' (\code{combined1} with \code{sigmaInter} locked to 0).
#' @inheritParams ModelError
#' @include ModelError.R
#' @return An S7 object of class \code{Proportional}.
#' @export
Proportional = new_class(
  "Proportional", package = "PFIM", parent = ModelError,
  validator = function( self ) {
    if ( length( prop( self, "sigmaInter" ) ) != 1L ||
         is.na( prop( self, "sigmaInter" ) ) ||
         prop( self, "sigmaInter" ) != 0 )
      return( "Proportional: sigmaInter must be 0." )
    if ( !identical( prop( self, "varianceForm" ), "combined1" ) )
      return( "Proportional: varianceForm must be 'combined1'." )
    NULL
  },
  constructor = function( output          = character( 0 ),
                          sigmaInter      = 0.0,
                          sigmaSlope      = 0.0,
                          sigmaInterFixed = FALSE,
                          sigmaSlopeFixed = FALSE,
                          cError          = 1.0,
                          ... ) {
    if ( length( sigmaInter ) && ( is.na( sigmaInter ) || sigmaInter != 0 ) )
      stop( "Proportional(): sigmaInter must be 0.", call. = FALSE )
    a = .pfimConstructModelErrorArgs(
      output, 0.0, sigmaSlope, sigmaInterFixed, sigmaSlopeFixed,
      cError, "combined1", ...,
      .forceSigmaInter = 0.0
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
