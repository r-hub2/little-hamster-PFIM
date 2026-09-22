#' @title Evaluation
#' @description
#' Design evaluation: build the model, compute the Fisher information matrix (FIM),
#' and store results per design. Inherits \code{PFIMProject} fields directly.
#' Use \code{\link{prop}} for project fields on both \code{Evaluation} and
#' \code{Optimization}. After \code{run()}, \code{\link{getEvaluationDesign}} returns
#' per-design results when several designs were evaluated.
#'
#' Equation strings may use declared \code{\link{ModelParameter}} names,
#' \code{t}, \code{dose_*} / \code{Tinf_*} for administered outcomes, and other
#' outcome names. Free symbols resolve in \code{knitr::knit_global()} (while
#' knitting) or \code{globalenv()}. Names that already exist on the search path
#' (\code{beta}, \code{gamma}, \code{pi}) are not reported as missing; declare
#' every model parameter explicitly.
#' @return An \code{Evaluation} object with \code{evaluationDesign} and \code{fim} filled.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' }
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
                                              modelClass = character(0),
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
                                              numberOfOccasions = NA_real_,
                                              fim = NULL ) {
                         if ( is.null( fim ) )
                           fim = .placeholderFim( fimType )
                         new_object(
                           PFIMProject(
                             name = name,
                             modelClass = modelClass,
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
                             odeSolverParameters = odeSolverParameters,
                             numberOfOccasions = numberOfOccasions
                           ),
                           evaluationDesign = evaluationDesign
                         )
                       })
S4_register( Evaluation )

#' Build a default FIM instance from type label.
#'
#' Used by the \code{Evaluation} constructor so \code{fim} is never \code{NULL}
#' before \code{defineFim()} / \code{run()}. Empty labels default to
#' \code{PopulationFim}; unknown labels error via \code{.pfimInstantiateFimType}.
#' @param fimType \code{fimType} label (may be empty before \code{defineFim()}).
#' @return A \code{PopulationFim}, \code{IndividualFim}, or \code{BayesianFim} object.
#' @noRd
#' @keywords internal
.placeholderFim = function( fimType = character(0) ) {
  if ( !.pfimIsNonEmptyScalar( fimType ) )
    return( PopulationFim() )
  .pfimInstantiateFimType( fimType )
}

#' Extract FIM matrix blocks from an evaluation.
#'
#' Returns the named blocks of the project (or selected design) FIM rather than
#' the S7 object itself - convenient for scripts that only need matrices.
#' @param evaluation An \code{Evaluation} object.
#' @param ... Optional multi-design selector (\code{NULL} keeps the primary FIM).
#' @usage getFim(evaluation, ...)
#' @name getFim
#' @return List with \code{fisherMatrix}, \code{fixedEffects}, \code{varianceEffects}.
#' @export
getFim = new_generic( "getFim", c( "evaluation" ) )

#' Per-design evaluation result from a multi-design \code{run()}.
#' @param evaluation An \code{Evaluation} object after \code{run()}.
#' @param design Design index (integer) or name (character).
#' @return One element of \code{evaluationDesign} (evaluated design with \code{fim}).
#' @export
getEvaluationDesign = function( evaluation, design = 1L ) {
  .getEvaluationDesignEntry( evaluation, design )
}

#' Names of the deepest elements in a nested list (used for model-type detection).
#'
#' @noRd
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

#' Run design evaluation.
#'
#' Pipeline: set FIM cache scope -> rebuild model (with finite differences) ->
#' instantiate FIM type -> \code{evaluateDesign()} per design -> post-process the
#' primary FIM via \code{setEvaluationFim()}. Additional designs remain in
#' \code{evaluationDesign}; use \code{getEvaluationDesign()} to retrieve them.
#' @param pfimproject An \code{Evaluation} object.
#' @return The same \code{Evaluation} with \code{evaluationDesign} and \code{fim} updated.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' }
#' @name run
#' @export

method( run, Evaluation ) = function( pfimproject )
{
  # Reuse an outer cache scope when nested (e.g. constraint-grid batch eval).
  prev_scope = pfim_get_option( "fim.cache.scope", NULL )
  had_scope = .pfimIsNonEmptyScalar( prev_scope )
  eval_scope = if ( had_scope ) prev_scope else .pfimProjectScopeId( pfimproject )
  if ( !had_scope ) {
    pfim_set_option( fim.cache.scope = eval_scope )
    on.exit( pfim_set_option( fim.cache.scope = NULL ), add = TRUE )
  }

  # Batch mode skips invalidation so repeated cells share model/gradient caches.
  if ( !isTRUE( pfim_get_option( "eval.batch", FALSE ) ) ) {
    .invalidateEvalModelCache( pfimproject )
    .pfimClearGradientPerfCaches()
    .pfimCovOccasionCacheClear( eval_scope )
  }
  .pfimSetCacheEvaluation( pfimproject )
  on.exit( pfim_set_option( fim.cache.evaluation = NULL ), add = TRUE )
  model = rebuildEvalModel( pfimproject, finiteDifference = TRUE )

  fim = defineFim( pfimproject )

  designs = prop( pfimproject, "designs" )
  evaluationDesign = map(
    designs,
    function( design ) evaluateDesign( design, model, fim )
  )
  prop( pfimproject, "evaluationDesign" ) = evaluationDesign

  # Primary slot mirrors design 1 for getSE/getRSE/... convenience accessors.
  prop( pfimproject, "fim" ) = setEvaluationFim(
    prop( pluck( evaluationDesign, 1L ), "fim" ), pfimproject
  )

  return( pfimproject )
}

#' @name getFim
#' @keywords internal

method( getFim, Evaluation ) = function( evaluation, design = NULL )
{
  # design = NULL -> project-level fim; otherwise index/name into evaluationDesign.
  fim = if ( is.null( design ) ) {
    prop( evaluation, "fim" )
  } else {
    prop( .getEvaluationDesignEntry( evaluation, design ), "fim" )
  }
  fisherMatrix = prop( fim, "fisherMatrix" )
  fixedEffects = prop( fim, "fixedEffects" )
  varianceEffects = prop( fim, "varianceEffects" )

  return( list( fisherMatrix = fisherMatrix, fixedEffects = fixedEffects, varianceEffects = varianceEffects ) )
}

#' Resolve one evaluated design entry by index or name.
#'
#' Character selectors match against the *input* \code{designs} names (not the
#' evaluated copies), then index into \code{evaluationDesign}.
#' @noRd
#' @keywords internal
.getEvaluationDesignEntry = function( evaluation, design = 1L ) {
  evaluated = prop( evaluation, "evaluationDesign" )
  if ( !length( evaluated ) )
    .pfimStop( "no evaluationDesign: call run() first." )
  if ( is.character( design ) ) {
    idx = match(
      design,
      map_chr( prop( evaluation, "designs" ), ~ prop( .x, "name" ) )
    )
    if ( is.na( idx ) )
      .pfimStop( "unknown design name: ", design )
    design = idx
  }
  pluck( evaluated, design )
}
