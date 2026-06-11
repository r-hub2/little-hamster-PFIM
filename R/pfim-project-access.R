#' @keywords internal
projectOf = function( x ) {
  if ( S7::S7_inherits( x, Optimization ) )
    return( prop( x, "project" ) )
  x
}

#' @keywords internal
projectProp = function( x, name ) {
  prop( projectOf( x ), name )
}

#' @keywords internal
`projectProp<-` = function( x, name, value ) {
  if ( S7::S7_inherits( x, Optimization ) ) {
    proj <- prop( x, "project" )
    prop( proj, name ) <- value
    prop( x, "project" ) <- proj
  } else {
    prop( x, name ) <- value
  }
  invisible( x )
}

#' @keywords internal
.applyFlatToArms = function( flat_pos, arms ) {
  idx = 1L
  out = vector( "list", length( arms ) )
  for ( i in seq_along( arms ) ) {
    arm = arms[[ i ]]
    sts = prop( arm, "samplingTimes" )
    for ( j in seq_along( sts ) ) {
      n = length( prop( sts[[ j ]], "samplings" ) )
      prop( sts[[ j ]], "samplings" ) = flat_pos[ idx:( idx + n - 1L ) ]
      idx = idx + n
    }
    prop( arm, "samplingTimes" ) = sts
    out[[ i ]] = arm
  }
  out
}

#' @keywords internal
.evaluationFromProject = function( project, designs = projectProp( project, "designs" ),
                                    name = "...." ) {
  Evaluation(
    name                    = name,
    modelEquations          = projectProp( project, "modelEquations" ),
    modelFromLibrary        = projectProp( project, "modelFromLibrary" ),
    modelParameters         = projectProp( project, "modelParameters" ),
    modelError              = projectProp( project, "modelError" ),
    modelCovariates         = projectProp( project, "modelCovariates" ),
    modelCovariatesEquation = projectProp( project, "modelCovariatesEquation" ),
    designs                 = designs,
    fimType                 = projectProp( project, "fimType" ),
    fim                     = projectProp( project, "fim" ),
    outputs                 = projectProp( project, "outputs" ),
    odeSolverParameters     = projectProp( project, "odeSolverParameters" )
  )
}

#' Build an \code{Evaluation} from an \code{Optimization} and a single \code{Design}.
#' @keywords internal
.evaluationFromOptimization = function( optimization, design, name = NULL ) {
  if ( is.null( name ) ) {
    opt_name = projectProp( optimization, "name" )
    name = if ( length( opt_name ) == 1L && nzchar( opt_name ) )
      paste0( opt_name, "_", prop( design, "name" ) ) else ""
  }
  .evaluationFromProject( optimization, designs = list( design ), name = name )
}
