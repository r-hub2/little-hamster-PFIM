# Session-scoped PFIM options (not global R options()).
#
# Options live in `.pfimSession` for the life of the R process. Defaults are
# applied lazily via `.pfimDefaultOption()` when a key is unset. Use
# `pfim_set_option()` / `pfim_get_option()` from user code; call
# `pfim_reset_session()` between unrelated runs to clear options and caches.

.pfimSession = new.env( parent = emptyenv() )

# Strip optional "PFIM." prefix so callers can use either naming style.
.pfimNormalizeOptionName = function( name ) {
  sub( "^PFIM\\.", "", name )
}

# Names accepted by pfim_set_option() (with or without a PFIM. prefix).
.pfimKnownOptionNames = c(
  "fim.cache", "fim.cache.hits", "fim.cache.maxEntries",
  "fim.cache.scope", "fim.cache.evaluation",
  "covariate.occasion.cache", "covariate.occasion.cache.hits",
  "model.cache.signature", "model.detect.legacy_warn",
  "verbose", "eval.batch",
  "perf.adminCache", "perf.fdCache", "perf.fdLinearOnly",
  "perf.odeTimesCache", "perf.odeTimesCache.maxEntries",
  "devPath", "constraints.maxTasks"
)

# Built-in defaults when the session env has no stored value for `key`.
.pfimDefaultOption = function( key, missing = NULL ) {
  switch( key,
    fim.cache = TRUE,
    fim.cache.hits = 0L,
    fim.cache.maxEntries = 2048L,
    covariate.occasion.cache = TRUE,
    covariate.occasion.cache.hits = 0L,
    model.cache.signature = TRUE,
    verbose = FALSE,
    eval.batch = FALSE,
    perf.adminCache = TRUE,
    perf.fdCache = TRUE,
    perf.fdLinearOnly = TRUE,
    perf.odeTimesCache = TRUE,
    perf.odeTimesCache.maxEntries = 2048L,
    devPath = NULL,
    fim.cache.scope = NULL,
    constraints.maxTasks = NULL,
    missing )
}

# ggplot2 `facet_wrap(space=)` needs 4.0.0; cache the version check (DESCRIPTION
# is otherwise re-read on every SE/RSE facet).
.pfimGgplot2SupportsFacetSpace = function() {
  key = "ggplot2.supportsFacetSpace"
  if ( exists( key, envir = .pfimSession, inherits = FALSE ) )
    return( isTRUE( .pfimSession[[ key ]] ) )
  ok = utils::packageVersion( "ggplot2" ) >= "4.0.0"
  assign( key, ok, envir = .pfimSession )
  ok
}

# Warn at most once per (flag, id) in the session. `id` distinguishes callers
# (hash of the tiny mus, etc.) so a second badly scaled design still warns.
.pfimWarnOnce = function( flag, ..., id = NULL ) {
  key = paste0( "warned.", flag )
  if ( !is.null( id ) )
    key = paste0( key, ".", as.character( id ) )
  if ( isTRUE( .pfimSession[[ key ]] ) )
    return( invisible( NULL ) )
  assign( key, TRUE, envir = .pfimSession )
  warning( ..., call. = FALSE )
  invisible( NULL )
}

#' Read a PFIM session option
#'
#' @param name Option name, with or without \code{PFIM.} prefix.
#' @param default Value when the option is unset.
#' @return Option value.
#' @seealso \code{\link{pfim_set_option}}, \code{\link{pfim_reset_session}}
#' @export
pfim_get_option = function( name, default = NULL ) {
  key = .pfimNormalizeOptionName( name )
  if ( exists( key, envir = .pfimSession, inherits = FALSE ) )
    return( get( key, envir = .pfimSession ) )
  .pfimDefaultOption( key, default )
}

#' Set PFIM session options
#'
#' Options apply for the current R session until \code{\link{pfim_reset_session}}.
#'
#' @param ... Named option values. Common options:
#' \describe{
#'   \item{\code{fim.cache}}{Cache FIM results during optimization (default \code{TRUE}).}
#'   \item{\code{fim.cache.maxEntries}}{LRU cap on design and covariate/occasion cache
#'     entries (default \code{2048L}). Set \code{NULL} for no limit.}
#'   \item{\code{eval.batch}}{Batch model rebuilds on constraint grids (default \code{FALSE};
#'     set automatically during \code{\link{generateFimsFromConstraints}}).}
#'   \item{\code{constraints.maxTasks}}{\code{NULL} (default) evaluates every cell of the
#'     dose \eqn{\times} sampling grid in discrete optimization. Set to a positive integer
#'     to cap the number of FIM evaluations: cells are subsampled deterministically
#'     (evenly spaced within each dose stratum), so Fedorov-Wynn and Multiplicative runs
#'     remain reproducible without touching the global RNG.}
#'   \item{\code{perf.adminCache}, \code{perf.fdCache}, \code{perf.odeTimesCache}}{ODE / FD caches
#'     (default \code{TRUE}).}
#'   \item{\code{perf.odeTimesCache.maxEntries}}{LRU cap on the C++ ODE sim-time cache
#'     (default \code{2048L}). Set \code{NULL} for no limit.}
#'   \item{\code{verbose}}{Extra messages (default \code{FALSE}).}
#' }
#' @return Invisibly \code{NULL}.
#' @seealso \code{\link{pfim_get_option}}, \code{\link{pfim_reset_session}}
#' @export
pfim_set_option = function( ... ) {
  dots = list( ... )
  nms  = names( dots )
  if ( is.null( nms ) || any( is.na( nms ) | !nzchar( nms ) ) )
    .pfimStop( "pfim_set_option() requires named arguments." )
  iwalk( dots, function( val, nm ) {
    key = .pfimNormalizeOptionName( nm )
    if ( !key %in% .pfimKnownOptionNames )
      .pfimStop(
        "unknown option '", nm, "'. See ?pfim_set_option for valid names."
      )
    assign( key, val, envir = .pfimSession )
  } )
  invisible( NULL )
}

#' Reset PFIM session options and clear all computed-data caches.
#'
#' Clears \code{.pfimSession} options and empties the design FIM cache
#' (\code{.pfimFimDesignCache}), covariate/occasion cache, eval-model rebuild cache,
#' and gradient/ODE performance caches (FD scheme, admin, ODE time grids).
#' Call between unrelated runs in a long-lived R session (Shiny, services).
#' @param closeDevices Logical. If \code{TRUE}, close all open graphics devices
#'   except the null device. Default \code{FALSE}.
#' @return Invisibly \code{NULL}.
#' @export
pfim_reset_session = function( closeDevices = FALSE ) {
  rm( list = ls( .pfimSession, all.names = TRUE ), envir = .pfimSession )
  .pfimClearFimCaches()
  if ( isTRUE( closeDevices ) )
    .pfimCloseGraphicsDevices()
  invisible( NULL )
}

#' Temporarily set the RNG seed and restore \code{.Random.seed} on exit.
#'
#' If \code{seed} is \code{NULL}, does nothing (does not call \code{set.seed(NULL)}).
#' @param seed Integer seed, or \code{NULL}.
#' @return A zero-argument restore function for use with \code{on.exit()}.
#' @noRd
#' @keywords internal
.pfimLocalSeed = function( seed ) {
  if ( is.null( seed ) )
    return( function() invisible( NULL ) )
  has_old = exists( ".Random.seed", envir = .GlobalEnv, inherits = FALSE )
  old = if ( has_old ) get( ".Random.seed", envir = .GlobalEnv ) else NULL
  set.seed( seed )
  function() {
    if ( is.null( old ) ) {
      if ( exists( ".Random.seed", envir = .GlobalEnv, inherits = FALSE ) )
        rm( ".Random.seed", envir = .GlobalEnv )
    } else {
      assign( ".Random.seed", old, envir = .GlobalEnv )
    }
    invisible( NULL )
  }
}
