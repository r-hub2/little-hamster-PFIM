#' Response plots for all designs in an evaluation
#' @name plotEvaluation
#' @export

method( plotEvaluation, Evaluation ) = function( pfimproject, plotOptions )
{
  designs = prop( pfimproject, "designs" )
  model = rebuildEvalModel( pfimproject, finiteDifference = TRUE )
  fim = defineFim( pfimproject )
  design = pluck( designs, 1 )
  designName = prop( design, "name" )
  arms = prop( design, "arms" )
  # generate and print all plots
  allPlots = map( arms, ~ processArmEvaluationResults( .x, model, fim, designName, plotOptions ) )
  allPlots = setNames( list( allPlots |> map( ~ .x[[designName]] ) |> list_flatten() ), designName )
  return( allPlots )
}

#' Sensitivity indices for the optimal design
#' @name plotSensitivityIndices
#' @export

method( plotSensitivityIndices, Evaluation ) = function( pfimproject, plotOptions )
{
  designs = prop( pfimproject, "designs" )
  model = rebuildEvalModel( pfimproject, finiteDifference = TRUE )
  fim = defineFim( pfimproject )
  design = pluck( designs, 1 )
  designName = prop( design, "name" )
  arms = prop( design, "arms" )

  # generate and print all plots
  allPlots = map( arms, ~ processArmEvaluationSI( .x, model, fim, designName, plotOptions ) )
  allPlots = setNames( list( allPlots |> map( ~ .x[[designName]] ) |> list_flatten() ), designName )
  return( allPlots )
}

#' SE barplot for the optimal design
#'
#' @description
#' Generates a bar plot showing the Standard Errors (SE) for the fixed effects
#' and variance components of the model. This visualization helps assess the
#' expected precision of the parameter estimates for the current design.
#'
#' @name plotSE
#' @param pfimproject An object of class \code{PFIMProject} containing the evaluation results.
#' @return A bar plot displaying the calculated SE for each model parameter.
#'
#' @examples
#' \dontrun{
#' # Assuming 'myPFIMproject' has been evaluated using run()
#'
#' # Generate the bar plot of Standard Errors
#' plotSE(myPFIMproject)
#' }
#' @export

# plot SE  from evaluation
method( plotSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  plotSE = plotSEFIM( fim, pfimproject )
  return( plotSE )
}

#' RSE barplot for the optimal design
#' @name plotRSE
#' @export

method( plotRSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot RSE
  fim = prop( pfimproject, "fim" )
  plotRSE = plotRSEFIM( fim, pfimproject )
  return( plotRSE )
}
