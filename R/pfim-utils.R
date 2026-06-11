# Internal helpers shared across PFIM modules.

.pfimEvalModelCache = new.env( parent = emptyenv() )

#' @keywords internal
.invalidateEvalModelCache = function( pfimproject ) {
  cacheId = .pfimModelCacheId( pfimproject )
  .pfimEvalModelCache[[ cacheId ]] = NULL
  invisible( NULL )
}

#' @keywords internal
.evaluateFimConstraintsCell = function(
    fimIndex,
    iterDose,
    iterComb,
    totalIterations,
    show_progress,
    evaluation,
    design,
    arms,
    dosesForDesign,
    samplingsForFIMs,
    designName,
    combinationGrid,
    baseModel = NULL,
    baseFim   = NULL ) {

  armsWithDoses = purrr::map( arms, function( arm ) {
    armName         = prop( arm, "name" )
    administrations = purrr::map( prop( arm, "administrations" ), function( adm ) {
      prop( adm, "dose" ) = dosesForDesign[[ armName ]][[ prop( adm, "outcome" ) ]][ iterDose ]
      adm
    } )
    prop( arm, "administrations" ) = administrations
    arm
  } )

  armsUpdated = purrr::map( armsWithDoses, function( arm ) {
    armName       = prop( arm, "name" )
    idx           = combinationGrid[ iterComb, armName ]
    samplingEntry = purrr::pluck( samplingsForFIMs, designName, armName, idx )
    prop( arm, "samplingTimes" ) = samplingEntry
    list(
      arm = arm,
      samplingsForFW = unlist( purrr::map( samplingEntry, ~ prop( .x, "samplings" ) ), use.names = FALSE )
    )
  } )

  armResult  = armsUpdated[[ 1L ]]
  tempDesign = design
  prop( tempDesign, "arms" ) = list( armResult$arm )

  if ( !is.null( baseModel ) && !is.null( baseFim ) ) {
    evaluatedDesign = evaluateDesign( tempDesign, baseModel, baseFim )
    fisherMatrix    = prop( prop( evaluatedDesign, "fim" ), "fisherMatrix" )
    evalResult      = .pfimEvaluationFromDesign( evaluation, tempDesign, evaluatedDesign )
  } else {
    tempEval = evaluation
    prop( tempEval, "designs" ) = list( tempDesign )
    evalResult   = .pfimRunEvaluationCached( tempEval )
    fisherMatrix = getFim( evalResult )$fisherMatrix
  }

  dimFim = nrow( fisherMatrix )
  dimVec = dimFim * ( dimFim + 1L ) / 2L

  fisherMatrixForAlgoFW = matrix(
    fisherMatrix[ rev( lower.tri( t( fisherMatrix ), diag = TRUE ) ) ],
    ncol  = dimVec,
    byrow = TRUE
  )

  if ( show_progress )
    message( sprintf( "FIM evaluation: %d / %d", fimIndex, totalIterations ) )

  list(
    armResult             = armResult,
    samplingsForFW        = armResult$samplingsForFW,
    fisherMatrixForAlgoFW = fisherMatrixForAlgoFW,
    fisherMatrix          = fisherMatrix,
    dimFim                = dimFim,
    cachedEvaluation      = evalResult
  )
}

#' @keywords internal
.as_df_rows = function( lst ) {
  if ( !length( lst ) ) return( as.data.frame( list() ) )
  map( lst, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
}

#' @keywords internal
.constraintsArmsTable = function( armsConstraints )
  map( armsConstraints, .as_df_rows ) |> list_rbind()

#' @keywords internal
.armConstraintsContinuous = function( arm ) {
  armName = prop( arm, "name" )
  armSize = prop( arm, "size" )
  map( prop( arm, "samplingTimesConstraints" ), function( sc ) {
    fmt = function( x ) paste0( "(", paste( x, collapse = ", " ), ")" )
    list(
      "Arms name"                  = armName,
      "Number of subjects"         = armSize,
      "Outcome"                    = prop( sc, "outcome" ),
      "Initial samplings"          = fmt( prop( sc, "initialSamplings" ) ),
      "Samplings windows"          = paste( map_chr( prop( sc, "samplingsWindows" ), ~ paste0( "(", paste( .x, collapse = "," ), ")" ) ), collapse = ", " ),
      "Number of times by windows" = fmt( prop( sc, "numberOfTimesByWindows" ) ),
      "Min sampling"               = fmt( prop( sc, "minSampling" ) )
    )
  } )
}

#' @keywords internal
rebuildEvalModel = function( pfimproject, finiteDifference = FALSE ) {

  cacheId = .pfimModelCacheId( pfimproject )
  if ( is.null( .pfimEvalModelCache[[ cacheId ]] ) )
    .pfimEvalModelCache[[ cacheId ]] = list()
  cacheKey = if ( isTRUE( finiteDifference ) ) "fd" else "plain"
  if ( !is.null( .pfimEvalModelCache[[ cacheId ]][[ cacheKey ]] ) )
    return( .pfimEvalModelCache[[ cacheId ]][[ cacheKey ]] )

  if ( length( projectProp( pfimproject, "modelFromLibrary" ) ) != 0L ) {
    projectProp( pfimproject, "modelEquations" ) <-
      defineModelEquationsFromLibraryOfModel( pfimproject )
  }

  model = defineModelType( pfimproject )
  if ( isTRUE( finiteDifference ) )
    model = finiteDifferenceHessian( model )
  model = defineModelWrapper( model, pfimproject )
  if ( usesCovariateOccasionStructure( model ) )
    model = defineCovariatesData( model )
  ensureModelOutputNames( model, pfimproject )

  .pfimEvalModelCache[[ cacheId ]][[ cacheKey ]] = model
  model
}

#' @keywords internal
ensureModelOutputNames = function( model, pfimproject ) {
  if ( length( prop( model, "outputNames" ) ) > 0L ) return( model )
  outputs = prop( pfimproject, "outputs" )
  if ( length( outputs ) == 0L ) return( model )
  on = unname( unlist( outputs, use.names = FALSE ) )
  if ( length( on ) == 0L ) on = names( outputs )
  prop( model, "outputNames" ) = on
  model
}
