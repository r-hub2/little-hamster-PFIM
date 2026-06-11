# Observation variance and sigma derivatives for the population FIM lambda block.

evaluateModelVarianceCore = function( model, evaluationModel ) {

  modelErrors = prop( model, "modelError"  )
  outputNames = prop( model, "outputNames" )

  # Error model derivatives, one entry per output (NULL entries dropped).
  errorDerivativesList = map( modelErrors, function( err ) {
    outcome = prop( err, "output" )
    if ( outcome %in% outputNames )
      evaluateErrorModelDerivatives( err, evaluationModel[[ outcome ]][ , outcome ] )
    else NULL
  }) |> compact() |> set_names( outputNames )

  errorVariance  = map( errorDerivativesList, ~ bdiag( .x$errorVariance ) ) |> bdiag()
  totalSamplings = map_int( evaluationModel, ~ length( .x[["time"]] ) ) |> sum()

  samplingOffsets = accumulate(
    outputNames,
    function( offset, outName ) offset + length( evaluationModel[[ outName ]][["time"]] ),
    .init = 0L
  ) |> head( -1L ) |> set_names( outputNames )

  sigmaDerivatives = outputNames |>
    map( function( outName ) {
      n   = length( evaluationModel[[ outName ]][["time"]] )
      rng = ( samplingOffsets[[ outName ]] + 1L ):( samplingOffsets[[ outName ]] + n )

      map( errorDerivativesList[[ outName ]]$sigmaDerivatives, function( derivComp ) {
        mat             = matrix( 0, nrow = totalSamplings, ncol = totalSamplings )
        mat[ rng, rng ] = derivComp
        mat
      })
    }) |> list_flatten()

  list( errorVariance = errorVariance, sigmaDerivatives = sigmaDerivatives )
}

method( evaluateModelVariance, Model ) = function( model, arm ) {
  evaluationModel = prop( arm, "evaluationModel" )
  if ( !usesCovariateOccasionStructure( model ) ) {
    evaluateModelVarianceCore( model, evaluationModel )
  } else {
    evaluateModelVarianceWithCovariates( model, evaluationModel )
  }
}

evaluateModelVarianceWithCovariates = function( model, evaluationModelWithCovariates ) {
  map( evaluationModelWithCovariates, function( combinationData ) {
    list(
      combination = combinationData$combination,
      proportion  = combinationData$proportion,
      variances   = map( combinationData$evaluations, function( occasionData ) {
        list(
          occasion = occasionData$occasion,
          variance = evaluateModelVarianceCore( model, occasionData$evaluation )
        )
      })
    )
  })
}
