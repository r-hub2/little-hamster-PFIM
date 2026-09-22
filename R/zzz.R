# Package load hooks and S7 method registration.
#
# On load: optionally record the source tree as `devPath` (when developing from
# a checkout), then ensure S4 class tables exist and call `S7::methods_register()`.
#
# Top-level `S4_register()` next to each `new_class()` is still required so
# `method<-` works while sourcing the package. Lazy-load does not replay those
# top-level side effects on a fresh install load, so `.onLoad` re-registers
# every class that participates in S4 generics (e.g. `show`) before wiring
# the deferred S7 methods.

#' Package \code{.onLoad}: set \code{devPath} when present, register S7 methods.
#' @noRd
#' @keywords internal
.onLoad = function(libname, pkgname) {
  ns = asNamespace( pkgname )
  # Parent of the installed package path is the source tree during load_all().
  pkg_root = normalizePath( file.path( getNamespaceInfo( ns, "path" ), ".." ), winslash = "/" )
  if ( file.exists( file.path( pkg_root, "DESCRIPTION" ) ) &&
       is.null( pfim_get_option( "devPath" ) ) )
    pfim_set_option( devPath = pkg_root )

  # Replay S4 tables (lazy-load skips top-level S4_register side effects).
  S4_register( Evaluation )
  S4_register( Optimization )
  S4_register( CovariateTest )
  S4_register( PopulationFim )
  S4_register( IndividualFim )
  S4_register( BayesianFim )
  S4_register( MultiplicativeAlgorithm )
  S4_register( FedorovWynnAlgorithm )
  S4_register( PSOAlgorithm )
  S4_register( PGBOAlgorithm )
  S4_register( SimplexAlgorithm )

  S7::methods_register()
}

# Silence R CMD check notes for magrittr/purrr placeholder `.`.
utils::globalVariables(c("."))
