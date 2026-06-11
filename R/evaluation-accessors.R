#' Fisher matrix blocks for the optimal design
#' @name getFisherMatrix
#' @export

method( getFisherMatrix, Evaluation ) = function( pfimproject )
{
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  fisherMatrix = prop( fim, "fisherMatrix" )
  fixedEffects = prop( fim, "fixedEffects" )
  varianceEffects = prop( fim, "varianceEffects" )

  return( list( fisherMatrix = fisherMatrix, fixedEffects = fixedEffects, varianceEffects = varianceEffects ) )
}


method( show, Evaluation ) = function( object )
{
  fim = setEvaluationFim( prop( object, "fim" ), object )
  prop( object, "fim" ) = fim
  .showPopulationFimConsole( fim, object )
}

#' Response plots for all designs in an evaluation
#' @name plotEvaluation
#' @export


method( getSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  SEAndRSE = prop( fim, "SEAndRSE" )
  SE = SEAndRSE$SE
  return( SE )
}

#' Relative standard errors from the optimal design FIM
#' @name getRSE
#' @export

method( getRSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  SEAndRSE = prop( fim, "SEAndRSE" )
  RSE = SEAndRSE$RSE
  return( RSE )
}

#' Parameter shrinkage from the Bayesian FIM
#' @name getShrinkage
#' @export

method( getShrinkage, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  shrinkage = prop( fim, "shrinkage" )
  return( shrinkage )
}

#' Determinant of the optimal design FIM
#' @name getDeterminant
#' @export

method( getDeterminant, Evaluation ) = function( pfimproject )
{
  fisherMatrix = getFisherMatrix( pfimproject )
  return( det( fisherMatrix$fisherMatrix ) )
}

#' D-criterion of the optimal design
#' @name getDcriterion
#' @export

method( getDcriterion, Evaluation ) = function( pfimproject )
{
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  return( Dcriterion( fim ) )
}

#' Parameter correlations from the optimal design FIM
#' @name getCorrelationMatrix
#' @export

method( getCorrelationMatrix, Evaluation ) = function( pfimproject )
{
  fim       = setEvaluationFim( prop( pfimproject, "fim" ), pfimproject )
  M         = prop( fim, "fisherMatrix" )
  covMatrix = .safeCholInv( M )
  cov2cor( covMatrix )
}
