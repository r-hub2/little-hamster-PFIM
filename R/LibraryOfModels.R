#' @title LibraryOfModels
#' @description
#' Registry of built-in PK and PD model definitions.
#'
#' Subclasses \code{LibraryOfPKModels} and \code{LibraryOfPDModels} store named
#' equation lists; \code{defineModelEquationsFromLibraryOfModel()} looks them up
#' via the package instances \code{pkModelLibrary} / \code{pdModelLibrary}.
#' @param models List of PK and PD model metadata objects.
#' @return An S7 object of class \code{LibraryOfModels}.
#' @examples
#' \donttest{
#' length(prop(pkModelLibrary, "models"))
#' length(prop(pdModelLibrary, "models"))
#' }
#' @export

LibraryOfModels = new_class("LibraryOfModels", package = "PFIM",

                       properties = list(
                         models = new_property( class_list, default = list())))

#' Look up one named equation list from a PK or PD library.
#' @noRd
#' @keywords internal
.pfimLibraryEquation = function( models, name, slot ) {
  if ( !is.character( name ) || !.pfimIsNonEmptyScalar( name ) )
    .pfimStop( "modelFromLibrary$", slot, " must be a single library name." )
  eq = models[[ name ]]
  if ( is.null( eq ) )
    .pfimStop( "Library model '", name, "' was not found in ", slot, "." )
  eq
}

