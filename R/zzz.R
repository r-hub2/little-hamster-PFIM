.registerS7WithS4 = function(pkgname) {
  ns = asNamespace(pkgname)
  for (nm in c(
    "BayesianFim", "CovariateTest", "Evaluation", "FedorovWynnAlgorithm",
    "IndividualFim", "MultiplicativeAlgorithm", "Optimization",
    "PGBOAlgorithm", "PopulationFim", "PSOAlgorithm", "SimplexAlgorithm"
  )) {
    if (exists(nm, envir = ns, inherits = FALSE)) {
      S4_register(get(nm, envir = ns))
    }
  }
}

.onLoad = function(libname, pkgname) {
  if ( getRversion() < "4.4.0" ) {
    `%||%` = function(x, y) if ( is.null( x ) ) y else x
    assign( "%||%", `%||%`, envir = asNamespace( pkgname ) )
  }
  ns = asNamespace( pkgname )
  pkg_root = normalizePath( file.path( getNamespaceInfo( ns, "path" ), ".." ), winslash = "/" )
  if ( file.exists( file.path( pkg_root, "DESCRIPTION" ) ) &&
       is.null( pfim_get_option( "devPath" ) ) )
    pfim_set_option( devPath = pkg_root )
  .registerS7WithS4( pkgname )
  S7::methods_register()
}

utils::globalVariables(c("."))

