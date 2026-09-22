# User-facing vs developer stops.
#
# .pfimStop()  - reaches Evaluation / Optimization / run() callers.
# .pfimWarn()  - same audience.
# .pfimInternalStop() - package invariants; tests may hit these via :::.

#' User-facing abort (no call stack, no internal function name).
#' @noRd
#' @keywords internal
.pfimStop = function( ... ) {
  stop( paste0( "PFIM: ", paste0( ..., collapse = "" ) ), call. = FALSE )
}

#' User-facing warning (no call stack, no internal function name).
#' @noRd
#' @keywords internal
.pfimWarn = function( ... ) {
  warning( paste0( "PFIM: ", paste0( ..., collapse = "" ) ), call. = FALSE )
}

#' Developer abort for broken internals (still no \code{.fn:} leak).
#' @noRd
#' @keywords internal
.pfimInternalStop = function( ... ) {
  stop( paste0( ..., collapse = "" ), call. = FALSE )
}
