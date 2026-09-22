# Observation variance + sigma derivatives for the population FIM lambda block.
#
# Multi-output residual blocks are stacked with Matrix::bdiag. Per-output
# dV/dsigma slices are embedded as sparse diagonals so large designs
# (many samples x outputs) never allocate dense zero fills during assembly.

#' Core residual variance for one (flat) evaluationModel list.
#' @noRd
#' @keywords internal
evaluateModelVarianceCore = function( model, evaluationModel ) {
  modelErrors = prop( model, "modelError" )
  outputNames = prop( model, "outputNames" )

  # One evaluateErrorModelDerivatives call per configured outcome.
  errs = keep( modelErrors, ~ prop( .x, "output" ) %in% outputNames )
  errorDerivativesList = set_names(
    map( errs, function( err ) {
      out = prop( err, "output" )
      evaluateErrorModelDerivatives(
        err,
        evaluationModel[[ out ]][ , out, drop = TRUE ]
      )
    } ),
    map_chr( errs, ~ prop( .x, "output" ) )
  )

  # Require a residual-error model for every outcome. A missing modelError used
  # to fall back to Diagonal(n) (= R = I), silently matching a unit residual.
  missingOut = setdiff( outputNames, map_chr( errs, ~ prop( .x, "output" ) ) )
  if ( length( missingOut ) )
    stop(
      "No modelError for output(s): ", paste( missingOut, collapse = ", " ),
      ". Provide a residual-error model (Constant, Proportional, Combined1/2).",
      call. = FALSE
    )
  if ( !length( errs ) && length( outputNames ) )
    stop(
      "modelError is empty but the model has outputs. ",
      "Provide at least one residual-error model.",
      call. = FALSE
    )

  blocks = map( outputNames, function( outName ) {
    as.matrix( errorDerivativesList[[ outName ]]$errorVariance )
  } )
  errorVariance = if ( length( blocks ) )
    bdiag( blocks ) else Matrix::Diagonal( n = 0L )

  n_by_out = map_int( outputNames, ~ length( evaluationModel[[ .x ]][[ "time" ]] ) )
  totalSamplings = sum( n_by_out )
  samplingOffsets = set_names(
    c( 0L, cumsum( n_by_out )[ -length( n_by_out ) ] ),
    outputNames
  )

  # Flatten per-output estimable sigma derivatives into FIM column order.
  sigmaDerivatives = outputNames |>
    map( function( outName ) {
      offset = samplingOffsets[[ outName ]]
      map(
        pluck( errorDerivativesList, outName, "sigmaDerivatives", .default = list() ),
        ~ .pfimEmbedDiagonalBlock( .x, offset, totalSamplings )
      )
    } ) |>
    list_flatten()

  list( errorVariance = errorVariance, sigmaDerivatives = sigmaDerivatives )
}

method( evaluateModelVariance, Model ) = function( model, arm ) {
  evaluationModel = prop( arm, "evaluationModel" )
  if ( !usesCovariateOccasionStructure( model ) )
    evaluateModelVarianceCore( model, evaluationModel )
  else
    evaluateModelVarianceWithCovariates( model, evaluationModel )
}

#' Residual variance under covariate x occasion nesting.
#' @noRd
#' @keywords internal
evaluateModelVarianceWithCovariates = function( model, evaluationModelWithCovariates ) {
  map( evaluationModelWithCovariates, ~ list(
    combination = .x$combination,
    proportion  = .x$proportion,
    variances   = map( .x$evaluations, ~ list(
      occasion = .x$occasion,
      variance = evaluateModelVarianceCore( model, .x$evaluation )
    ) )
  ) )
}
