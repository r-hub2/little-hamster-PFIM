#' @title Constant
#' @description
#' Pure additive residual error: \code{V = sigmaInter^2}
#' (\code{combined1} with \code{sigmaSlope} locked to 0).
#' @inheritParams ModelError
#' @include ModelError.R
#' @return An S7 object of class \code{Constant}.
#' @export
Constant = new_class(
  "Constant", package = "PFIM", parent = ModelError,
  # Validator also blocks prop<- mutation after construction.
  validator = function( self ) {
    if ( length( prop( self, "sigmaSlope" ) ) != 1L ||
         is.na( prop( self, "sigmaSlope" ) ) ||
         prop( self, "sigmaSlope" ) != 0 )
      return( "Constant: sigmaSlope must be 0." )
    if ( !identical( prop( self, "varianceForm" ), "combined1" ) )
      return( "Constant: varianceForm must be 'combined1'." )
    NULL
  },
  constructor = function( output          = character( 0 ),
                          sigmaInter      = 0.0,
                          sigmaSlope      = 0.0,
                          sigmaInterFixed = FALSE,
                          sigmaSlopeFixed = FALSE,
                          cError          = 1.0,
                          ... ) {
    # Fail fast with a clear API message (validator is the lasting guard).
    if ( length( sigmaSlope ) && ( is.na( sigmaSlope ) || sigmaSlope != 0 ) )
      stop( "Constant(): sigmaSlope must be 0.", call. = FALSE )
    a = .pfimConstructModelErrorArgs(
      output, sigmaInter, 0.0, sigmaInterFixed, sigmaSlopeFixed,
      cError, "combined1", ...,
      .forceSigmaSlope = 0.0
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
