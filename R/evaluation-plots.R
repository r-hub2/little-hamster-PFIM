#' Nested response or SI plots for every design in an evaluation.
#'
#' Rebuilds the model once (no FD here - processors that need gradients build
#' their own stencil). Arms come from \code{.pfimEvaluatedArmsForDesign()} so
#' multi-design runs plot the evaluated schedules, not the unevaluated inputs.
#' @noRd
#' @keywords internal
.pfimPlotAllDesigns = function( pfimproject, plotOptions, processor ) {
  designs = prop( pfimproject, "designs" )
  # FD is deferred to SI processors; response path stays prediction-only.
  model   = rebuildEvalModel( pfimproject, finiteDifference = FALSE )
  fim     = defineFim( pfimproject )
  stats::setNames(
    map( designs, function( design ) {
      designName = prop( design, "name" )
      plots      = map(
        .pfimEvaluatedArmsForDesign( pfimproject, design ),
        ~ processor( .x, model, fim, designName, plotOptions )
      )
      list_flatten( map( plots, ~ .x[[ designName ]] ) )
    }),
    map_chr( designs, ~ prop( .x, "name" ) )
  )
}

#' Predicted responses with sampling markers for all designs/arms.
#' @name plotEvaluation
#' @usage NULL
#' @export
method( plotEvaluation, Evaluation ) = function( pfimproject, plotOptions = list() )
  .pfimPlotAllDesigns( pfimproject, plotOptions, processArmEvaluationResults )

#' Sensitivity indices over time for all designs/arms.
#' @name plotSensitivityIndices
#' @usage NULL
#' @export
method( plotSensitivityIndices, Evaluation ) = function( pfimproject, plotOptions = list() )
  .pfimPlotAllDesigns( pfimproject, plotOptions, processArmEvaluationSI )

#' SE bar chart for the primary evaluation FIM.
#' @name plotSE
#' @usage NULL
#' @export
method( plotSE, Evaluation ) = function( pfimproject )
  plotSEFIM( prop( pfimproject, "fim" ), pfimproject )

#' RSE bar chart for the primary evaluation FIM.
#' @name plotRSE
#' @usage NULL
#' @export
method( plotRSE, Evaluation ) = function( pfimproject )
  plotRSEFIM( prop( pfimproject, "fim" ), pfimproject )
