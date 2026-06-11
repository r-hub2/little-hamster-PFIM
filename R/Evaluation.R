#' @title Evaluation
#' @description
#' Design evaluation: build the model, compute the Fisher information matrix (FIM),
#' and store results per design.
#' @inheritParams PFIMProject
#' @param evaluationDesign List of evaluated designs (filled by \code{run()}).
#' @param modelCovariatesEquation See \code{\link{PFIMProject}} (\code{modelCovariatesEquation}).
#' @include PFIMProject.R
#' @include Model.R
#' @export

Evaluation = new_class("Evaluation", package = "PFIM", parent = PFIMProject,
                       properties = list(
                         evaluationDesign = new_property(class_list, default = list())
                       ),
                       constructor = function(evaluationDesign = list(),
                                              name = character(0),
                                              modelParameters = list(),
                                              modelCovariates = list(),
                                              modelCovariatesEquation = character(0),
                                              modelEquations = list(),
                                              modelFromLibrary = list(),
                                              modelError = list(),
                                              designs = list(),
                                              outputs = list(),
                                              fimType = character(0),
                                              odeSolverParameters = list(),
                                              fim = NULL ) {
                         if ( is.null( fim ) )
                           fim = .placeholderFim( fimType )
                         new_object(
                           .parent = PFIMProject(
                             name = name,
                             modelEquations = modelEquations,
                             modelCovariatesEquation = modelCovariatesEquation,
                             modelFromLibrary = modelFromLibrary,
                             modelParameters = modelParameters,
                             modelCovariates = modelCovariates,
                             modelError = modelError,
                             designs = designs,
                             outputs = outputs,
                             fimType = fimType,
                             fim = fim,
                             odeSolverParameters = odeSolverParameters
                           ),
                           evaluationDesign = evaluationDesign
                         )
                       })
S4_register( Evaluation )

.placeholderFim = function( fimType = character(0) ) {
  if ( length( fimType ) == 1L && nzchar( fimType ) ) {
    switch(
      tolower( fimType ),
      population = PopulationFim(),
      individual = IndividualFim(),
      bayesian   = BayesianFim(),
      PopulationFim()
    )
  } else {
    PopulationFim()
  }
}

getFim = new_generic( "getFim", c( "evaluation" ) )

#' Names of the deepest elements in a nested list (used for model-type detection).
#'
#' @keywords internal

.getListLastName = function( list ) {
  if ( is.list( list ) ) {
    result = map( list, .getListLastName )
    result = unlist( result, recursive = FALSE )
    if ( length( result ) == 0 ) {
      return( names( list ) )
    } else {
      return( result )
    }
  }
}

#' Run design optimization
#' @name run
#' @export

method( run, Evaluation ) = function( pfimproject )
{
  if ( !isTRUE( pfim_get_option( "eval.batch", FALSE ) ) ) {
    .invalidateEvalModelCache( pfimproject )
    .pfimClearGradientPerfCaches()
  }
  model = rebuildEvalModel( pfimproject, finiteDifference = TRUE )

  fim = defineFim( pfimproject )

  # evaluate the designs
  designs = prop( pfimproject, "designs" )
  evaluationDesign = map(
    designs,
    function( design ) evaluateDesign( design, model, fim )
  )
  prop( pfimproject, "evaluationDesign" ) = evaluationDesign

  # results of the evaluation
  designNumber = 1
  evaluationDesign = pluck( evaluationDesign, designNumber )
  prop( pfimproject, "fim" ) = setEvaluationFim(
    prop( evaluationDesign, "fim" ), pfimproject
  )

  return( pfimproject )
}

#' Extract the FIM object from an evaluation
#' @name getFim
#' @export

method( getFim, Evaluation ) = function( evaluation )
{
  fim = prop( evaluation, "fim" )
  fisherMatrix = prop( fim, "fisherMatrix" )
  fixedEffects = prop( fim, "fixedEffects" )
  varianceEffects = prop( fim, "varianceEffects" )

  return( list( fisherMatrix = fisherMatrix, fixedEffects = fixedEffects, varianceEffects = varianceEffects ) )
}
