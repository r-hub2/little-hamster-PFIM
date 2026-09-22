# Shared analytic evaluation / administration helpers.
#
# The four ModelAnalytic* cores share a time x outcome grid, superposition of
# doses, and per-outcome sampling restriction. Keep the physics (bolus vs
# infusion, tau) in the callers; this file is the scaffolding.

#' Sorted unique observation grid across all arm sampling times.
#' @noRd
#' @keywords internal
.analyticSamplingGrid = function( arm ) {
  map( prop( arm, "samplingTimes" ), ~ prop( .x, "samplings" ) ) |>
    unlist() |> sort() |> unique()
}

#' Re-bind analytic wrappers to the current user / knit environment.
#' @noRd
#' @keywords internal
.analyticBindUserEnv = function( ... ) {
  fns = list( ... )
  user_env = .pfimUserSymbolEnv()
  lapply( fns, function( fn ) {
    if ( is.function( fn ) )
      environment( fn ) = user_env
    fn
  } )
}

#' Fill unused \code{t_*} formals on a passive analytic wrapper.
#'
#' Administered wrappers receive relative times since each dose. Passive
#' (no-admin) equations use calendar time at the current grid point, named
#' \code{t_<that outcome>} - not \code{t_<administered outcome>}.
#' @noRd
#' @keywords internal
.pfimFillMissingTimeFormals = function( fn, args, t ) {
  if ( !is.function( fn ) ) return( args )
  nms = names( formals( fn ) )
  if ( is.null( nms ) || !length( nms ) ) return( args )
  need = nms[ startsWith( nms, "t_" ) & !nms %in% names( args ) ]
  if ( !length( need ) ) return( args )
  for ( nm in need )
    args[[ nm ]] = t
  args
}

#' One outcome column group: administered prediction, optional passive/PD.
#' @noRd
#' @keywords internal
.analyticOutcomeFrame = function( evalAdmin, evalNoAdmin ) {
  if ( is.null( evalNoAdmin ) || length( evalNoAdmin ) == 0L )
    data.frame( Admin = evalAdmin )
  else
    data.frame( Admin = evalAdmin, NoAdmin = evalNoAdmin )
}

#' Restrict a dense timexoutcome table to each output's own sampling times.
#' @noRd
#' @keywords internal
.analyticRestrictBySampling = function( evaluationModelTmp, arm, outputNames ) {
  samplings_by_output = set_names(
    map( prop( arm, "samplingTimes" ), ~ prop( .x, "samplings" ) ),
    outputNames
  )
  set_names(
    map(
      outputNames,
      ~ evaluationModelTmp[
        evaluationModelTmp$time %in% samplings_by_output[[ .x ]],
        c( "time", .x )
      ]
    ),
    outputNames
  )
}

#' Time x outcome grid: \code{evalOutcome(iterTime, outcome, args)} returns
#' \code{list(admin, noAdmin, args)}.
#' @noRd
#' @keywords internal
.analyticEvalGrid = function( samplings, outcomes, args0, evalOutcome ) {
  rows = map( seq_along( samplings ), function( iterTime ) {
    evalAcc = reduce(
      outcomes,
      function( acc, outcome ) {
        piece = evalOutcome( iterTime, outcome, acc$args )
        acc$args = piece$args
        acc$results[[ length( acc$results ) + 1L ]] =
          .analyticOutcomeFrame( piece$admin, piece$noAdmin )
        acc
      },
      .init = list( args = args0, results = list() )
    )
    data.frame(
      time = samplings[[ iterTime ]],
      do.call( cbind, evalAcc$results )
    )
  } )
  list_rbind( rows )
}

#' Label columns and restrict to per-outcome sampling schedules.
#' @noRd
#' @keywords internal
.analyticFinishEvaluation = function( tmp, outputNames, arm ) {
  colnames( tmp ) = c( "time", outputNames )
  .analyticRestrictBySampling( tmp, arm, outputNames )
}

#' Covariate/occasion dispatch shared by every analytic \code{evaluateModel}.
#' @noRd
#' @keywords internal
.analyticDispatchEvaluate = function( model, arm, core ) {
  if ( usesCovariateOccasionStructure( model ) )
    evaluateModelWithCovariates( model, arm, core )
  else
    core( model, arm )
}

#' Relative-time dose table for bolus / oral superposition (optional tau grid).
#' @noRd
#' @keywords internal
.analyticRelativeDoseTable = function( samplings, timeDose, dose, tau ) {
  if ( tau != 0 ) {
    timeDose = seq( 0, max( samplings ), tau )
    dose     = rep( dose, length( timeDose ) )
  }
  timeRel = timeDose |>
    map( ~ ifelse( samplings - .x > 0, samplings - .x, samplings ) ) |>
    reduce( cbind )
  indicesDoses = if ( is.null( dim( timeRel ) ) ) {
    1L
  } else {
    map_int( seq_len( nrow( timeRel ) ), ~ length( unique( timeRel[ .x, ] ) ) )
  }
  list(
    data = data.frame( timeRel, indicesDoses ),
    dose = dose,
    tau  = tau
  )
}

#' Infusion window table: during/after flags, relative times, active dose index.
#' @noRd
#' @keywords internal
.analyticInfusionWindowTable = function( samplings, timeDose, dose, Tinf, tau ) {
  if ( tau != 0 ) {
    timeDose = seq( 0, max( samplings ), tau )
    n        = length( timeDose )
    dose     = rep( dose, n )
    Tinf     = rep( Tinf, n )
  }
  duringAndAfter = rep( "afterInfusion", length( samplings ) )
  Tinfs = map2( timeDose, timeDose + Tinf, c )
  samplingsDuringInfusion = map( Tinfs, function( iv ) {
    samplings |> keep( ~ .x >= min( iv ) & .x < max( iv ) )
  } ) |> unlist() |> unique()
  duringAndAfter[ samplings %in% samplingsDuringInfusion ] = "duringInfusion"

  samplingTimeDoses = timeDose |>
    map( ~ ifelse( samplings - .x > 0, samplings - .x, 0 ) )
  indicesDoses = map_int( samplings, function( s ) {
    idx = which( s >= timeDose )
    idx[ length( idx ) ]
  } )

  data = data.frame( duringAndAfter, indicesDoses, samplings, samplingTimeDoses )
  colnames( data ) = c(
    "duringAndAfter", "indicesDoses", "samplings",
    paste0( "samplingTimeDoses", seq_along( dose ) )
  )
  list( data = data, dose = dose, Tinf = Tinf, tau = tau )
}

#' Administered infusion prediction: current during-window plus finished after.
#' @noRd
#' @keywords internal
.analyticEvalInfusionAdmin = function( duringAfter, indicesDoses, sampTimes,
                                       dose, Tinf, argsCall,
                                       tName, doseName, TinfName,
                                       fnDuring, fnAfter ) {
  set3 = function( args, t, d, Tinf_i ) {
    args[[ tName ]]    = t
    args[[ doseName ]] = d
    args[[ TinfName ]] = Tinf_i
    args
  }
  if ( duringAfter == "duringInfusion" ) {
    if ( indicesDoses == 1L ) {
      argsCall = set3( argsCall, sampTimes[[ 1L ]], dose[[ 1L ]], Tinf[[ 1L ]] )
      return( do.call( fnDuring, argsCall )[[ 1L ]] )
    }
    idxBefore = seq_len( indicesDoses - 1L )
    argsCall  = set3(
      argsCall, sampTimes[ idxBefore ], dose[ idxBefore ], Tinf[ idxBefore ]
    )
    sumAfter = sum( do.call( fnAfter, argsCall )[[ 1L ]] )
    argsCall = set3(
      argsCall,
      sampTimes[[ indicesDoses ]],
      dose[[ indicesDoses ]],
      Tinf[[ indicesDoses ]]
    )
    return( sumAfter + do.call( fnDuring, argsCall )[[ 1L ]] )
  }
  idxAll   = seq_len( indicesDoses )
  argsCall = set3( argsCall, sampTimes[ idxAll ], dose[ idxAll ], Tinf[ idxAll ] )
  sum( do.call( fnAfter, argsCall )[[ 1L ]] )
}

#' Relative sampling times stored in an infusion solver-input data frame.
#' @noRd
#' @keywords internal
.analyticInfusionSampleTimes = function( data, iterTime ) {
  sampCols = grep( "samplingTimeDoses", colnames( data ), value = TRUE )
  as.numeric( data[ iterTime, sampCols ] )
}
