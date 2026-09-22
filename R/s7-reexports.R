#' S7 \code{prop} / \code{prop<-} reexports with Optimization project delegation.
#'
#' For \code{Optimization} objects, project-field names listed in
#' \code{.pfimProjectFieldNames()} read and write the nested \code{project} slot
#' via \code{projectOf()} / \code{projectProp<-}. All other properties use
#' \code{S7::prop} unchanged.
#' @name s7-reexports
#' @include PFIMProject.R
#' @include Optimization.R
#' @include pfim-project-access.R
#' @keywords internal
NULL

#' Property accessor (S7); project fields on \code{Optimization} read slot \code{project}.
#' @param object An S7 object.
#' @param name Property name.
#' @name prop
#' @return The requested property value, or the modified object for replacement.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' prop(ev, "name")
#' prop(ev, "name") = "renamed"
#' }
#' @export
prop = function( object, name ) {
  if ( .pfimDelegatesProjectProp( object, name ) )
    return( S7::prop( projectOf( object ), name ) )
  S7::prop( object, name )
}

#' @param object An S7 object.
#' @param name Property name.
#' @param value Value to assign.
#' @name prop
#' @export
`prop<-` = function( object, name, value ) {
  if ( .pfimDelegatesProjectProp( object, name ) ) {
    object = `projectProp<-`( object, name, value )
    return( object )
  }
  S7::prop( object, name ) = value
  object
}
