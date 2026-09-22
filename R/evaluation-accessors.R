# Evaluation accessors: always refresh dimnames / SE-RSE via setEvaluationFim()
# so getSE/getRSE/... stay consistent after model rebuilds. Not a raw prop read.

#' Post-process the evaluation FIM (labels, SE/RSE) before accessors / show.
#' @noRd
#' @keywords internal
.fimForEvaluation = function( evaluation ) {
  setEvaluationFim( prop( evaluation, "fim" ), evaluation )
}

#' Fisher matrix blocks from an evaluated design.
#'
#' Runs \code{setEvaluationFim()} so dimnames and SE/RSE match the current model
#' metadata, then returns the matrix blocks used by report helpers, plus
#' \code{singularFim}.
#' @return List with \code{fisherMatrix}, \code{fixedEffects},
#'   \code{varianceEffects}, and \code{singularFim}.
#' @name getFisherMatrix
#' @export

method( getFisherMatrix, Evaluation ) = function( pfimproject ) {
  fim = .fimForEvaluation( pfimproject )
  list(
    fisherMatrix    = prop( fim, "fisherMatrix" ),
    fixedEffects    = prop( fim, "fixedEffects" ),
    varianceEffects = prop( fim, "varianceEffects" ),
    singularFim     = isTRUE( prop( fim, "singularFim" ) )
  )
}

#' Console summary of an evaluated design (FIM, SE, RSE, criteria).
#' @param object First argument of generic.
#' @return Invisibly returns the printed \code{Evaluation} object.
#' @name show-methods
#' @keywords internal
method( show, Evaluation ) = function( object ) {
  fim = .fimForEvaluation( object )
  prop( object, "fim" ) = fim
  showFIM( fim )
  invisible( object )
}

#' Standard errors from the evaluated design FIM.
#' @name getSE
#' @export
method( getSE, Evaluation ) = function( pfimproject ) {
  prop( .fimForEvaluation( pfimproject ), "SEAndRSE" )$SE
}

#' Relative standard errors (\%) from the evaluated design FIM.
#' @name getRSE
#' @export

method( getRSE, Evaluation ) = function( pfimproject ) {
  prop( .fimForEvaluation( pfimproject ), "SEAndRSE" )$RSE
}

#' Parameter shrinkage from a Bayesian evaluation FIM (empty for other types).
#' @name getShrinkage
#' @export

method( getShrinkage, Evaluation ) = function( pfimproject ) {
  prop( .fimForEvaluation( pfimproject ), "shrinkage" )
}

#' Determinant of the evaluated design FIM.
#' @name getDeterminant
#' @export

method( getDeterminant, Evaluation ) = function( pfimproject ) {
  .fimDeterminant( prop( .fimForEvaluation( pfimproject ), "fisherMatrix" ) )
}

#' D-criterion of the evaluated design.
#' @name getDcriterion
#' @export

method( getDcriterion, Evaluation ) = function( pfimproject ) {
  Dcriterion( .fimForEvaluation( pfimproject ) )
}

#' Parameter correlations from the evaluated design FIM.
#' @name getCorrelationMatrix
#' @export

method( getCorrelationMatrix, Evaluation ) = function( pfimproject ) {
  .fimCorrelationMatrix( prop( .fimForEvaluation( pfimproject ), "fisherMatrix" ) )
}
