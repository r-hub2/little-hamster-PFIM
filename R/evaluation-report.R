#' HTML report for an optimization run
#' @name Report
#' @export

method( Report, Evaluation ) = function( pfimproject, outputPath, outputFile, plotOptions )
{
  projectName       = prop( pfimproject, "name"    )
  evaluationOutputs = prop( pfimproject, "outputs" )

  # Model equations (same pipeline as run(): library models must be resolved first)
  model          = rebuildEvalModel( pfimproject, finiteDifference = FALSE )
  modelEquations = prop( model, "modelEquations" )

  # Model error
  modelErrorData = prop( pfimproject, "modelError" ) |>
    map( getModelErrorData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( modelErrorData ) = c( "Output", "Type", "$\\sigma_{slope}$", "$\\sigma_{inter}$" )
  modelErrorTable = kbl( modelErrorData, align = c( "c","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE, position = "center", font_size = 13 )

  # Model parameters table (with IOV column)
  modelParametersTable = .buildModelParametersKable( prop( pfimproject, "modelParameters" ) )

  # Covariates: structure table + CovariateTest tables
  modelCovariates = prop( pfimproject, "modelCovariates" )
  hasCov          = length( modelCovariates ) > 0L

  covariatesTable     = .buildCovariatesKable( modelCovariates )          # NULL if no covariates
  covariateTestTables = if ( hasCov ) .buildCovariateTestSection( pfimproject ) else NULL

  # Arm tables
  designs      = prop( pfimproject, "designs" )
  designsNames = map_chr( designs, "name" )
  arms         = map( designs, ~ prop( .x, "arms" ) )
  armsData     = list_flatten( map( pluck( arms, 1 ), getArmData ) )

  # Administration table
  administrationData = list_flatten( map( pluck( arms, 1 ), armAdministration ) ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( administrationData ) = c( "Design name", "Arms name", "Number of subjects",
                                      "Outcome", "Dose", "Time dose", "$\\tau$", "$T_{inf}$" )
  administrationTable = kbl( administrationData, align = c( "l","l","l","c","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE, position = "center", font_size = 13 )

  # Initial design table
  initialDesignData = armsData |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( initialDesignData ) = c( "Arms name", "Number of subjects", "Outcome", "Dose", "Sampling times" )
  initialDesignTable = kbl( initialDesignData, align = c( "l","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE, position = "center", font_size = 13 )

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
  reportTables = list(
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
  )

  generateReportEvaluation( fim, reportTables,
                            outputFile = outputFile,
                            outputPath = outputPath )
}
