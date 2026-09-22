# =============================================================================
# Residual-error S7 hierarchy (ModelError and helpers)
# =============================================================================
#
# Estimable residual sigma rule (must stay identical in FIM column builders,
# Rcpp kernels, and R fallbacks):
#
#   estimable  <=>  value != 0  &&  !fixed
#
# Non-zero fixed sigmas still enter V (observation variance) but never become
# FIM parameters / SE-RSE rows. Zero sigmas are never estimable.
#
# Variance forms (property varianceForm):
#   combined1 : V = (a + b f^c)^2     Constant / Proportional / Combined1
#   combined2 : V = a^2 + (b f^c)^2   Combined2 (PopED additive + prop. SDs)
#
# Hot path: evaluateErrorModelDerivatives -> .pfimErrorModelDerivatives (Rcpp).
# =============================================================================

#' Validation message for residual SD constraints, or \code{NULL} if OK.
#' @noRd
#' @keywords internal
.modelErrorValidationMsg = function( sigmaInter, sigmaSlope, cError = 1.0 ) {
  if ( ( length( sigmaInter ) && ( is.na( sigmaInter ) || sigmaInter < 0 ) ) ||
       ( length( sigmaSlope ) && ( is.na( sigmaSlope ) || sigmaSlope < 0 ) ) )
    return( "ModelError: sigmaInter and sigmaSlope must be non-negative." )
  if ( length( cError ) && ( is.na( cError ) || cError <= 0 ) )
    return( "ModelError: cError must be positive." )
  NULL
}

#' Stop when residual SD constraints fail (constructors / helpers).
#' @noRd
#' @keywords internal
.checkModelErrorSigmas = function( sigmaInter, sigmaSlope, cError = 1.0 ) {
  msg = .modelErrorValidationMsg( sigmaInter, sigmaSlope, cError )
  if ( !is.null( msg ) ) stop( msg, call. = FALSE )
  invisible( NULL )
}

#' Whether a residual SD is an estimable FIM parameter.
#'
#' Single source of truth for names, values, and dV/dsigma assembly.
#' @noRd
#' @keywords internal
.pfimSigmaIsEstimable = function( value, fixed ) {
  isTRUE( value != 0 ) && !isTRUE( fixed )
}

#' Allowed \code{varianceForm} tokens.
#' @noRd
#' @keywords internal
.pfimModelErrorVarianceForms = function() c( "combined1", "combined2" )

#' Validation message for \code{varianceForm}, or \code{NULL} if OK.
#' @noRd
#' @keywords internal
.modelErrorVarianceFormMsg = function( varianceForm ) {
  if ( length( varianceForm ) != 1L || !is.character( varianceForm ) ||
       !varianceForm %in% .pfimModelErrorVarianceForms() )
    return( "ModelError: varianceForm must be 'combined1' or 'combined2'." )
  NULL
}

#' Absorb removed constructor args (\code{equation}, \code{derivatives}).
#' @noRd
#' @keywords internal
.pfimModelErrorLegacyDots = function( ... ) {
  dots = list( ... )
  if ( !length( dots ) ) return( invisible( NULL ) )
  nms = names( dots )
  if ( is.null( nms ) || any( is.na( nms ) | !nzchar( nms ) ) )
    stop( "ModelError: unnamed arguments are not allowed.", call. = FALSE )
  dead = c( "equation", "derivatives" )
  legacy = intersect( nms, dead )
  if ( length( legacy ) )
    warning(
      "ModelError: ", paste( legacy, collapse = ", " ),
      " ignored (removed from the class API).",
      call. = FALSE
    )
  extra = setdiff( nms, dead )
  if ( length( extra ) )
    stop(
      "ModelError: unused argument(s) ", paste( extra, collapse = ", " ), ".",
      call. = FALSE
    )
  invisible( NULL )
}

#' Normalize and validate constructor fields for every ModelError subclass.
#'
#' Call from subclass constructors only; then pass a parent \emph{instance}
#' into a literal \code{new_object(ModelError(...), ...)} so S7 keeps the
#' concrete class (S7 requires a direct \code{new_object()} call in the
#' constructor body - do not wrap it in a helper). S7 devel (#409) rejects a
#' class object as parent: use \code{ModelError(...)} or \code{S7_object()},
#' not \code{ModelError} itself. The parent is the first positional argument
#' (\code{.parent} on CRAN S7, \code{_parent} on S7 devel).
#'
#' @param .forceSigmaInter,.forceSigmaSlope Optional locks (Constant / Proportional).
#' @noRd
#' @keywords internal
.pfimConstructModelErrorArgs = function(
    output           = character( 0 ),
    sigmaInter       = 0.0,
    sigmaSlope       = 0.0,
    sigmaInterFixed  = FALSE,
    sigmaSlopeFixed  = FALSE,
    cError           = 1.0,
    varianceForm     = "combined1",
    ...,
    .forceSigmaInter = NULL,
    .forceSigmaSlope = NULL ) {
  .pfimModelErrorLegacyDots( ... )
  if ( !is.null( .forceSigmaInter ) ) sigmaInter = .forceSigmaInter
  if ( !is.null( .forceSigmaSlope ) ) sigmaSlope = .forceSigmaSlope
  .checkModelErrorSigmas( sigmaInter, sigmaSlope, cError )
  if ( length( output ) > 1L )
    stop( "ModelError: output must be a single string.", call. = FALSE )
  if ( .pfimIsBlankScalar( output ) )
    stop( "ModelError: output must be a non-empty string.", call. = FALSE )
  formMsg = .modelErrorVarianceFormMsg( varianceForm )
  if ( !is.null( formMsg ) ) stop( formMsg, call. = FALSE )
  list(
    output          = output,
    sigmaInter      = as.numeric( sigmaInter ),
    sigmaSlope      = as.numeric( sigmaSlope ),
    sigmaInterFixed = isTRUE( sigmaInterFixed ),
    sigmaSlopeFixed = isTRUE( sigmaSlopeFixed ),
    cError          = as.numeric( cError ),
    varianceForm    = varianceForm
  )
}

#' @title ModelError
#' @description
#' Base class for residual error models (\code{Constant}, \code{Combined1},
#' \code{Combined2}, \code{Proportional}).
#'
#' Variance is determined by \code{sigmaInter} (additive SD),
#' \code{sigmaSlope} (proportional SD), \code{cError}, and
#' \code{varianceForm}:
#' \itemize{
#'   \item \code{"combined1"}: \eqn{V=(a + b f^{c})^2}
#'   \item \code{"combined2"}: \eqn{V=a^2 + (b f^{c})^2} (PopED-style)
#' }
#' A sigma is estimable in the FIM iff it is non-zero and not marked fixed
#' (\code{sigmaInterFixed} / \code{sigmaSlopeFixed}).
#' @param output Outcome name.
#' @param sigmaInter Additive residual SD.
#' @param sigmaSlope Proportional residual SD.
#' @param sigmaInterFixed If \code{TRUE}, \code{sigmaInter} is not estimated.
#' @param sigmaSlopeFixed If \code{TRUE}, \code{sigmaSlope} is not estimated.
#' @param cError Power on the proportional prediction term.
#' @param varianceForm \code{"combined1"} or \code{"combined2"} (subclasses set this).
#' @param ... Legacy \code{equation}/\code{derivatives} warn and are ignored.
#' @return An S7 object of class \code{ModelError}.
#' @export
ModelError = new_class(
  "ModelError", package = "PFIM",
  properties = list(
    output          = new_property( class_character, default = "output" ),
    sigmaInter      = new_property( class_double,    default = 0.0 ),
    sigmaSlope      = new_property( class_double,    default = 0.0 ),
    sigmaInterFixed = new_property( class_logical,   default = FALSE ),
    sigmaSlopeFixed = new_property( class_logical,   default = FALSE ),
    cError          = new_property( class_double,    default = 1.0 ),
    varianceForm    = new_property( class_character, default = "combined1" )
  ),
  validator = function( self ) {
    .modelErrorValidationMsg(
      prop( self, "sigmaInter" ),
      prop( self, "sigmaSlope" ),
      prop( self, "cError" )
    ) %||% .modelErrorVarianceFormMsg( prop( self, "varianceForm" ) )
  },
  constructor = function( output          = "output",
                          sigmaInter      = 0.0,
                          sigmaSlope      = 0.0,
                          sigmaInterFixed = FALSE,
                          sigmaSlopeFixed = FALSE,
                          cError          = 1.0,
                          varianceForm    = "combined1",
                          ... ) {
    a = .pfimConstructModelErrorArgs(
      output, sigmaInter, sigmaSlope, sigmaInterFixed, sigmaSlopeFixed,
      cError, varianceForm, ...
    )
    new_object(
      S7_object(),
      output          = a$output,
      sigmaInter      = a$sigmaInter,
      sigmaSlope      = a$sigmaSlope,
      sigmaInterFixed = a$sigmaInterFixed,
      sigmaSlopeFixed = a$sigmaSlopeFixed,
      cError          = a$cError,
      varianceForm    = a$varianceForm
    )
  }
)

#' Evaluate residual variance and dV/dsigma diagonals
#' @param modelError A \code{ModelError} object.
#' @param ... For \code{ModelError}: \code{evaluationModel} (numeric predictions).
#' @name evaluateErrorModelDerivatives
#' @return List with \code{errorVariance} and \code{sigmaDerivatives}.
#' @keywords internal
evaluateErrorModelDerivatives = new_generic(
  "evaluateErrorModelDerivatives", c( "modelError" )
)

#' Extract residual error settings for reports
#' @param modelError A \code{ModelError} object.
#' @param ... Optional method arguments.
#' @name getModelErrorData
#' @return One-row data frame of residual-error settings.
#' @export
getModelErrorData = new_generic( "getModelErrorData", c( "modelError" ) )

#' @rdname evaluateErrorModelDerivatives
#' @name evaluateErrorModelDerivatives
#' @usage NULL
#' @keywords internal
method( evaluateErrorModelDerivatives, ModelError ) = function(
    modelError, evaluationModel ) {
  # Polymorphic via varianceForm - subclasses need no method override.
  .pfimErrorModelDerivatives(
    prop( modelError, "sigmaInter" ),
    prop( modelError, "sigmaSlope" ),
    prop( modelError, "sigmaInterFixed" ),
    prop( modelError, "sigmaSlopeFixed" ),
    prop( modelError, "cError" ),
    evaluationModel,
    form = prop( modelError, "varianceForm" )
  )
}

#' @rdname getModelErrorData
#' @name getModelErrorData
#' @keywords internal
method( getModelErrorData, ModelError ) = function( modelError ) {
  data.frame(
    output       = prop( modelError, "output" ),
    type         = str_remove( class( modelError )[[ 1L ]], "^PFIM::" ),
    varianceForm = prop( modelError, "varianceForm" ),
    sigmaSlope   = as.character( prop( modelError, "sigmaSlope" ) ),
    sigmaInter   = as.character( prop( modelError, "sigmaInter" ) ),
    cError       = as.character( prop( modelError, "cError" ) ),
    stringsAsFactors = FALSE
  )
}
