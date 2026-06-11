.pfimSession <- new.env( parent = emptyenv() )

.pfimOptionDefaults <- list(
  devPath               = NULL,
  fim.cache             = TRUE,
  fim.cache.hits        = 0L,
  fim.cache.scope       = NULL,
  model.cache.signature = TRUE,
  constraints.maxTasks  = NULL,
  verbose               = FALSE,
  eval.batch            = FALSE,
  perf.adminCache       = TRUE,
  perf.fdCache          = TRUE,
  perf.odeTimesCache    = TRUE
)

.pfimNormalizeOptionName <- function( name ) {
  sub( "^PFIM\\.", "", name )
}

.pfimInitSession <- function() {
  for ( nm in names( .pfimOptionDefaults ) )
    assign( nm, .pfimOptionDefaults[[ nm ]], envir = .pfimSession )
  invisible( NULL )
}

#' Read a PFIM session option
#'
#' Session options replace \code{options(PFIM.*)} for cache, performance, and
#' optimization toggles. Names may be given with or without the \code{PFIM.} prefix.
#' @param name Option name (e.g. \code{"fim.cache"} or \code{"PFIM.fim.cache"}).
#' @param default Value when the option is unset.
#' @return Option value.
#' @export
pfim_get_option <- function( name, default = NULL ) {
  key = .pfimNormalizeOptionName( name )
  if ( exists( key, envir = .pfimSession, inherits = FALSE ) )
    return( get( key, envir = .pfimSession ) )
  if ( key %in% names( .pfimOptionDefaults ) )
    return( .pfimOptionDefaults[[ key ]] )
  default
}

#' Set PFIM session options
#'
#' @param ... Named option values (\code{fim.cache = FALSE}, etc.).
#' @return Invisibly \code{NULL}.
#' @export
pfim_set_option <- function( ... ) {
  dots = list( ... )
  for ( nm in names( dots ) ) {
    key = .pfimNormalizeOptionName( nm )
    assign( key, dots[[ nm ]], envir = .pfimSession )
  }
  invisible( NULL )
}

#' Reset PFIM session options to package defaults
#' @return Invisibly \code{NULL}.
#' @export
pfim_reset_session <- function() {
  rm( list = ls( .pfimSession, all.names = TRUE ), envir = .pfimSession )
  .pfimInitSession()
  invisible( NULL )
}

.pfimInitSession()
