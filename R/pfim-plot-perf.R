# Plot helpers for evaluation reports (response and sensitivity curves).
#
# Re-evaluate on a dense time grid (step 0.05), overlay design
# sampling times as red markers. FIM evaluation itself stays at design times.

#' Cap for densify via \code{updateSamplingTimes} (long horizons stay tractable).
#' @noRd
#' @keywords internal
.pfimDefaultPlotMaxPoints = function() 400L

#' Time grid densifier for plot re-evaluation (keeps design times, adds steps).
#'
#' Uses \code{seq(0, tmax, by = 0.05)} when that stays under
#' \code{maxPoints}; otherwise the step grows so the curve remains smooth.
#' @noRd
#' @keywords internal
.pfimDensePlotTimes = function(
    times,
    maxPoints = .pfimDefaultPlotMaxPoints(),
    minStep   = 0.05 ) {
  times = sort( unique( as.numeric( times ) ) )
  if ( !length( times ) ) return( times )
  tmax = max( times, 0 )
  nGrid = as.integer( floor( tmax / minStep ) ) + 1L
  if ( nGrid > maxPoints )
    minStep = tmax / max( maxPoints - 1L, 1L )
  sort( unique( c( times, seq( 0, tmax, by = minStep ) ) ) )
}

#' Clone an arm and densify its sampling grids for response / SI plots.
#' @noRd
#' @keywords internal
.pfimDensifyArmForPlots = function( arm ) {
  updateSamplingTimes( .pfimCloneS7( arm ), getSamplingData( arm ) )
}

#' Model bound to a (possibly densified) plot arm; optional FD stencil for SI.
#' @noRd
#' @keywords internal
.pfimPreparePlotModel = function( model, plotArm, needFd = FALSE ) {
  modelPlot = .pfimPrepareModelForEvaluation( .pfimCloneS7( model ), plotArm )
  if ( isTRUE( needFd ) ) {
    gp = prop( modelPlot, "parametersForComputingGradient" )
    if ( is.null( gp ) || is.null( gp$shifted ) )
      modelPlot = finiteDifferenceHessian( modelPlot )
  }
  modelPlot
}

#' Sampling grids for one arm (response or SI plots).
#' @noRd
#' @keywords internal
.pfimArmPlotSampling = function( arm, model, outputNames, plotOptions ) {
  outNames = unlist( outputNames, use.names = FALSE )
  list(
    samplingsByResponse = .pfimSamplingsByResponse(
      getSamplingData( arm ), outNames, model
    )
  )
}

#' Evaluated arms from \code{run()} when available (carry cached model/gradient results).
#' @noRd
#' @keywords internal
.pfimEvaluatedArmsForDesign = function( pfimproject, design ) {
  evaluated = prop( pfimproject, "evaluationDesign" )
  if ( length( evaluated ) ) {
    idx = match( prop( design, "name" ), map_chr( evaluated, ~ prop( .x, "name" ) ) )
    if ( !is.na( idx ) ) {
      arms = prop( evaluated[[ idx ]], "evaluationArms" )
      if ( length( arms ) ) return( arms )
    }
  }
  prop( design, "arms" )
}

#' Minimum number of distinct sampling times required for response/SI plots.
#' @noRd
#' @keywords internal
.pfimMinSamplingTimesForPlots = function() 2L

#' Whether an arm has enough sampling times to build response/SI plots.
#' @noRd
#' @keywords internal
.pfimArmPlotsEnabled = function( arm, model, outputNames ) {
  outNames = unlist( outputNames, use.names = FALSE )
  if ( !length( outNames ) ) return( FALSE )
  samplings = .pfimSamplingsByResponse( getSamplingData( arm ), outNames, model )
  counts = map_int(
    samplings,
    ~ length( unique( stats::na.omit( as.numeric( .x ) ) ) )
  )
  all( counts >= .pfimMinSamplingTimesForPlots() )
}

#' Empty nested plot list for one arm (skip plotting).
#' @noRd
#' @keywords internal
.pfimEmptyArmPlotResult = function( designName, arm ) {
  stats::setNames(
    list( stats::setNames( list( list() ), prop( arm, "name" ) ) ),
    designName
  )
}

#' Sampling grids keyed by response names (RespPK) when outputs map to state names (Cc).
#' @noRd
#' @keywords internal
.pfimSamplingsByResponse = function( samplingData, responseNames, model = NULL ) {
  outputFormula = if ( !is.null( model ) ) .getModelOutputFormulas( model ) else list()
  .odeSamplingsByOutput( samplingData$samplingTimes, responseNames, outputFormula )
}

#' Stored FIM response curves at design sampling times.
#' @noRd
#' @keywords internal
.pfimResponseCurvesForPlot = function( arm, outputNames ) {
  stored = prop( arm, "evaluationModel" )
  if ( !length( stored ) ) return( NULL )
  em = if ( .isNestedEvaluationModel( stored ) )
    .aggregateEvaluationModelForPlot( stored ) else stored
  outNames = unlist( outputNames, use.names = FALSE )
  out = map( outNames, function( out_name ) {
    df = em[[ out_name ]]
    if ( is.null( df ) || !nrow( df ) || !out_name %in% names( df ) )
      return( NULL )
    df[ , c( "time", out_name ), drop = FALSE ]
  } )
  names( out ) = outNames
  if ( !any( lengths( out ) ) ) NULL else out
}

#' Drop covariate beta columns from sensitivity-index parameter names.
#' @noRd
#' @keywords internal
.pfimSiParamNamesNoBeta = function( names ) {
  if ( !length( names ) ) return( names )
  names[
    !startsWith( names, "beta_" ) &
      !startsWith( names, "\u03b2_" )
  ]
}

#' Proportion-weighted gradient matrix for one model output (covariate structure).
#' @noRd
#' @keywords internal
.pfimAggregateGradientForOutput = function( model, arm, outName ) {
  allGradientsData = prop( arm, "evaluationGradients" )
  zeroGrad = pluck( allGradientsData, 1L, "gradients", 1L, "gradient", outName ) * 0
  matList = map( allGradientsData, function( combinationData ) {
    nOcc     = length( combinationData$gradients )
    meanGrad = reduce(
      map( combinationData$gradients, ~ .x$gradient[[ outName ]] ),
      `+`
    ) / nOcc
    combinationData$proportion * meanGrad
  } )
  reduce( matList, `+`, .init = zeroGrad )
}

#' Time column plus gradient columns for one SI plot frame.
#' @noRd
#' @keywords internal
.pfimSiFrame = function( sparse_t, grad ) {
  if ( is.null( sparse_t ) ) return( NULL )
  grad_cols = as.data.frame( grad, stringsAsFactors = FALSE )
  if ( is.null( colnames( grad_cols ) ) || !ncol( grad_cols ) ||
       nrow( grad_cols ) != length( sparse_t ) )
    return( NULL )
  data.frame( time = sparse_t, grad_cols, check.names = FALSE )
}

#' Stored FIM gradients at design sampling times.
#' @noRd
#' @keywords internal
.pfimSiCurvesForPlot = function( arm, model, outputNames ) {
  raw = prop( arm, "evaluationGradients" )
  if ( !length( raw ) ) return( NULL )
  outNames = unlist( outputNames, use.names = FALSE )
  em = prop( arm, "evaluationModel" )
  emAgg = if ( length( em ) && .isNestedEvaluationModel( em ) )
    .aggregateEvaluationModelForPlot( em ) else em

  frames = if ( .isNestedArmEvaluation( raw ) ) {
    map( outNames, function( out ) {
      .pfimSiFrame(
        emAgg[[ out ]]$time,
        .pfimAggregateGradientForOutput( model, arm, out )
      )
    } )
  } else {
    map2( outNames, raw, function( out, grad ) {
      .pfimSiFrame( emAgg[[ out ]]$time, grad )
    } )
  }
  names( frames ) = outNames
  if ( !any( lengths( frames ) ) ) NULL else frames
}

#' Convert freshly computed gradients to per-output data frames with time.
#' @noRd
#' @keywords internal
.pfimSiFramesFromGradients = function(
    evaluationModelGradient, model, arm, outputNames, parametersNames ) {
  outNames = unlist( outputNames, use.names = FALSE )
  stList = .pfimSamplingsByResponse( getSamplingData( arm ), outNames, model )
  if ( .isNestedArmEvaluation( evaluationModelGradient ) ) {
    armForGrad = arm
    prop( armForGrad, "evaluationGradients" ) = evaluationModelGradient
    colNames = colnames( .pfimAggregateGradientForOutput( model, armForGrad, outNames[[ 1L ]] ) )
    if ( is.null( colNames ) || !length( colNames ) )
      colNames = parametersNames
    frames = map( outNames, function( out ) {
      gmat = .pfimAggregateGradientForOutput( model, armForGrad, out )
      st = stList[[ out ]][ seq_len( nrow( gmat ) ) ]
      df = as.data.frame( gmat, stringsAsFactors = FALSE )
      colnames( df ) = colNames
      df$time = st
      df
    } )
    set_names( frames, outNames )
  } else {
    map( outNames, function( out ) {
      grad = evaluationModelGradient[[ out ]]
      df   = as.data.frame( grad, stringsAsFactors = FALSE )
      df$time = stList[[ out ]][ seq_len( nrow( df ) ) ]
      df
    } ) |> set_names( outNames )
  }
}
