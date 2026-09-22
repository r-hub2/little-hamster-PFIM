# Project field access for Evaluation / Optimization / PFIMProject.
#
# Evaluation stores project fields on the object itself. Optimization nests a
# `project` slot; `projectProp()` / `prop()` delegate reads and writes there.
# Also owns the pfim#N cache-scope counter (Evaluation.R collates before
# pfim-fim-cache.R) and helpers to flatten sampling vectors onto arms.

# Monotonic generator for optimization-scoped cache ids ("pfim#1", ...).
.pfimCacheState = new.env( parent = emptyenv() )
.pfimCacheState$gen = 0L
# Maps project object address -> assigned cache scope string.
.pfimProjectScopeRegistry = new.env( parent = emptyenv() )

#' Stable key for a project's object address (used by the scope registry).
#' @noRd
#' @keywords internal
.pfimProjectScopeKey = function( project ) {
  paste0( "scope_", rlang::obj_address( projectOf( project ) ) )
}

#' Allocate the next unique FIM cache scope id for an optimization run.
#' @noRd
#' @keywords internal
.pfimNextCacheScope = function() {
  .pfimCacheState$gen = .pfimCacheState$gen + 1L
  paste0( "pfim#", .pfimCacheState$gen )
}

#' Names of fields that live on PFIMProject (and on Evaluation / nested project).
#' @noRd
#' @keywords internal
.pfimProjectFieldNames = function() {
  c(
    "name", "modelClass", "modelEquations", "modelCovariatesEquation",
    "modelFromLibrary", "modelParameters", "modelCovariates", "modelError",
    "optimizer", "optimizerParameters", "outputs", "designs", "fimType", "fim",
    "odeSolverParameters", "numberOfOccasions"
  )
}

#' TRUE when \code{prop()} on an Optimization should read/write the nested project.
#' @noRd
#' @keywords internal
.pfimDelegatesProjectProp = function( object, name ) {
  S7::S7_inherits( object, Optimization ) && name %in% .pfimProjectFieldNames()
}

#' Nested \code{PFIMProject} for an \code{Optimization}
#' @param x An \code{Evaluation}, \code{Optimization}, or \code{PFIMProject} object.
#' @return The underlying \code{PFIMProject}.
#' @export
projectOf = function( x ) {
  if ( S7::S7_inherits( x, Optimization ) )
    return( S7::prop( x, "project" ) )
  x
}

#' Read a project field
#' @param x An \code{Evaluation}, \code{Optimization}, or \code{PFIMProject} object.
#' @param name Character name of the project property.
#' @return The requested project property value, or the modified object for replacement.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' projectProp(ev, "name")
#' projectProp(ev, "name") = "renamed"
#' }
#' @export
projectProp = function( x, name ) {
  S7::prop( projectOf( x ), name )
}

#' Write a project field
#' @param x An \code{Evaluation}, \code{Optimization}, or \code{PFIMProject} object.
#' @param name Character name of the project property.
#' @param value Value to assign.
#' @rdname projectProp
#' @export
`projectProp<-` = function( x, name, value ) {
  # Optimization: mutate nested project, then write it back (S7 props are not refs).
  if ( S7::S7_inherits( x, Optimization ) ) {
    proj = S7::prop( x, "project" )
    S7::prop( proj, name ) = value
    S7::prop( x, "project" ) = proj
  } else {
    S7::prop( x, name ) = value
  }
  invisible( x )
}

#' Spread a flat sampling-time vector onto cloned arms (continuous optimizers).
#'
#' Layout matches the concatenation of all samplingTimes$samplings across arms.
#' @noRd
#' @keywords internal
.applyFlatToArms = function( flat_pos, arms ) {
  # Walk arms left-to-right, consuming contiguous slices of the flat vector.
  out = vector( "list", length( arms ) )
  idx = 1L
  for ( a in seq_along( arms ) ) {
    arm = .pfimCloneS7( arms[[ a ]] )
    sts = prop( arm, "samplingTimes" )
    for ( s in seq_along( sts ) ) {
      n = length( prop( sts[[ s ]], "samplings" ) )
      prop( sts[[ s ]], "samplings" ) = flat_pos[ idx:( idx + n - 1L ) ]
      idx = idx + n
    }
    prop( arm, "samplingTimes" ) = sts
    out[[ a ]] = arm
  }
  out
}

#' Build an Evaluation that copies model/FIM settings from a project.
#' @noRd
#' @keywords internal
.evaluationFromProject = function( project, designs = projectProp( project, "designs" ),
                                    name = "" ) {
  Evaluation(
    name                    = name,
    modelEquations          = projectProp( project, "modelEquations" ),
    modelClass              = projectProp( project, "modelClass" ),
    modelFromLibrary        = projectProp( project, "modelFromLibrary" ),
    modelParameters         = projectProp( project, "modelParameters" ),
    modelError              = projectProp( project, "modelError" ),
    modelCovariates         = projectProp( project, "modelCovariates" ),
    modelCovariatesEquation = projectProp( project, "modelCovariatesEquation" ),
    designs                 = designs,
    fimType                 = projectProp( project, "fimType" ),
    fim                     = projectProp( project, "fim" ),
    outputs                 = projectProp( project, "outputs" ),
    odeSolverParameters     = projectProp( project, "odeSolverParameters" ),
    numberOfOccasions       = projectProp( project, "numberOfOccasions" )
  )
}

#' Evaluation for one design taken from an Optimization (post-search reporting).
#' @noRd
#' @keywords internal
.evaluationFromOptimization = function( optimization, design, name = NULL ) {
  if ( is.null( name ) ) {
    opt_name = projectProp( optimization, "name" )
    name = if ( .pfimIsNonEmptyScalar( opt_name ) )
      paste0( opt_name, "_", prop( design, "name" ) ) else ""
  }
  .evaluationFromProject( optimization, designs = list( design ), name = name )
}
