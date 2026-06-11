NULL
#' Evaluate a single constraint-grid cell and return the updated Evaluation.
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
    administrations = prop( arm, "administrations" )
    administrations = purrr::map( administrations, function( adm ) {
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
      samplingsForFW = unlist(
        purrr::map( samplingEntry, ~ prop( .x, "samplings" ) ),
        use.names = FALSE
      )
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
    tempEval   = evaluation
    prop( tempEval, "designs" ) = list( tempDesign )
    evalResult = .pfimRunEvaluationCached( tempEval )
    fim        = getFim( evalResult )
    fisherMatrix = fim$fisherMatrix
  }

  dimFim       = nrow( fisherMatrix )
  dimVec       = dimFim * ( dimFim + 1L ) / 2L

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

  if ( length( lst ) == 0L )

    return( as.data.frame( list() ) )

  map( lst, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()

}



#' @keywords internal

.constraintsArmsTable = function( armsConstraints ) {

  map( armsConstraints, .as_df_rows ) |> list_rbind()

}



#' @keywords internal

.armConstraintsContinuous = function( arm ) {

  armName = prop( arm, "name" )

  armSize = prop( arm, "size" )

  map( prop( arm, "samplingTimesConstraints" ), function( samplingConstraint ) {

    outcome = prop( samplingConstraint, "outcome" )

    initialSamplings = paste0(

      "(", paste( prop( samplingConstraint, "initialSamplings" ), collapse = ", " ), ")"

    )

    samplingsWindows = paste(

      map_chr(

        prop( samplingConstraint, "samplingsWindows" ),

        ~ paste0( "(", paste( .x, collapse = "," ), ")" )

      ),

      collapse = ", "

    )

    numberOfTimesByWindows = paste0(

      "(", paste( prop( samplingConstraint, "numberOfTimesByWindows" ), collapse = ", " ), ")"

    )

    minSampling = paste0(

      "(", paste( prop( samplingConstraint, "minSampling" ), collapse = ", " ), ")"

    )

    list(

      "Arms name"                  = armName,

      "Number of subjects"         = armSize,

      "Outcome"                    = outcome,

      "Initial samplings"          = initialSamplings,

      "Samplings windows"          = samplingsWindows,

      "Number of times by windows" = numberOfTimesByWindows,

      "Min sampling"               = minSampling

    )

  } )

}
