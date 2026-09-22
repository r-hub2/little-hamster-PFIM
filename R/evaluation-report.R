#' HTML evaluation report.
#'
#' Assembles model / design / FIM kables and response / SI / SE / RSE plots, then
#' dispatches to the FIM-type R Markdown template via
#' \code{generateReportEvaluation()}.
#' Administration / initial-design tables use design 1 only; response and SI
#' plots cover every design in \code{designs}.
#' @rdname Report
#' @name Report
#' @export

method( Report, Evaluation ) = function( pfimproject, outputPath, outputFile, plotOptions )
{
  projectName       = .pfimProjectNameOrDefault( pfimproject )
  evaluationOutputs = prop( pfimproject, "outputs" )

  # Same rebuild as run(), but without FD stencil (report needs predictions only).
  model          = rebuildEvalModel( pfimproject, finiteDifference = FALSE )
  modelEquations = prop( model, "modelEquations" )

  # Model error
  modelErrorTable = .buildModelErrorKable( prop( pfimproject, "modelError" ) )

  # Model parameters table (with IOV column)
  modelParametersTable = .buildModelParametersKable( prop( pfimproject, "modelParameters" ) )

  # Covariates: structure table + CovariateTest tables
  modelCovariates = prop( pfimproject, "modelCovariates" )
  hasCov          = length( modelCovariates ) > 0L

  covariatesTable     = .buildCovariatesKable( modelCovariates )          # NULL if no covariates
  covariateTestTables = if ( hasCov ) .buildCovariateTestSection( pfimproject ) else NULL

  # Arm tables (administration uses the first design; plots cover all designs)
  designs  = prop( pfimproject, "designs" )
  armsData = .armsDataFromEvaluation( pfimproject )

  administrationData = list_flatten( map(
    prop( designs[[ 1L ]], "arms" ),
    ~ armAdministration( .x, prop( designs[[ 1L ]], "name" ) )
  ) ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( administrationData ) = c( "Design name", "Arms name", "Number of subjects",
                                      "Outcome", "Dose", "Time dose", "$\\tau$", "$T_{inf}$" )
  administrationTable = .kblReportStyled(
    administrationData,
    align = .kblAlign( ncol( administrationData ), 3L )
  )

  # Initial design table
  initialDesignData = armsData |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( initialDesignData ) = c( "Arms name", "Number of subjects", "Outcome", "Dose", "Sampling times" )
  initialDesignTable = .kblReportStyled(
    initialDesignData,
    align = .kblAlign( ncol( initialDesignData ), 1L )
  )

  # Fisher matrix, SE/RSE tables
  fim                   = prop( pfimproject, "fim" )
  fim                   = setEvaluationFim( fim, pfimproject )
  fimInitialDesignTable = tablesForReport( fim, pfimproject )

  # Plots
  plotsEvaluation = plotEvaluation( pfimproject, plotOptions )
  plotsSI         = plotSensitivityIndices( pfimproject, plotOptions )
  plotSEBars      = plotSE( pfimproject )
  plotRSEBars     = plotRSE( pfimproject )

  # Assemble and render
  reportTables = c(
    list(
      evaluationOutputs      = evaluationOutputs,
      modelEquations         = modelEquations,
      modelErrorTable        = modelErrorTable,
      modelParametersTable   = modelParametersTable,
      covariatesTable        = covariatesTable,
      covariateTestTables    = covariateTestTables,
      administrationTable    = administrationTable,
      initialDesignTable     = initialDesignTable,
      fimInitialDesignTable  = fimInitialDesignTable,
      plotsEvaluation        = plotsEvaluation,
      plotSensitivityIndices = plotsSI,
      plotSE                 = plotSEBars,
      plotRSE                = plotRSEBars,
      fim                    = fim,
      pfimproject            = pfimproject,
      projectName            = projectName
    ),
    .covReportFlags( hasCov, covariateTestTables, covariatesTable )
  )

  generateReportEvaluation( fim, reportTables,
                            outputFile = outputFile,
                            outputPath = outputPath )
}
