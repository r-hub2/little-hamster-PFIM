# Eval-model cache and constraint-grid helpers.
#
# Supports discrete design optimization (`generateFimsFromConstraints`):
# clone arms per grid cell, assign doses/samplings, evaluate the FIM, and pack
# the lower triangle for Fedorov-Wynn. Also rebuilds/caches the project model
# used by Evaluation (`rebuildEvalModel`).

# Per-project cache of rebuilt evaluation models (plain vs finite-difference).
.pfimEvalModelCache = new.env( parent = emptyenv() )

#' Remove every binding from an environment (cache flush helper).
#' @noRd
#' @keywords internal
.pfimEnvClear = function( env ) {
  rm( list = ls( env, all.names = TRUE ), envir = env )
  invisible( NULL )
}

#' Scalar non-empty string (length 1, not \code{NA}, not \code{""}).
#'
#' \code{nzchar(NA)} is TRUE, so callers must not use \code{nzchar()} alone
#' inside \code{if()}.
#' @noRd
#' @keywords internal
.pfimIsNonEmptyScalar = function( x ) {
  length( x ) == 1L && !is.na( x ) && nzchar( x )
}

#' Length-1 string that is \code{NA} or empty (a provided blank, not missing).
#' @noRd
#' @keywords internal
.pfimIsBlankScalar = function( x ) {
  length( x ) == 1L && !.pfimIsNonEmptyScalar( x )
}

#' Safe project display name for HTML reports.
#'
#' Handles \code{NULL}, \code{character(0)}, and \code{NA} - all of which make
#' \code{nzchar()} unsafe inside \code{if()}.
#' @noRd
#' @keywords internal
.pfimProjectNameOrDefault = function( pfimproject, default = "PFIM Report" ) {
  nm = tryCatch( prop( pfimproject, "name" ), error = function( e ) NULL )
  if ( !.pfimIsNonEmptyScalar( nm ) ) default else as.character( nm )
}

#' Extract atol/rtol from an ODE solver parameter list (defaults 1e-8).
#'
#' Defaults are deliberately stricter than a typical FD step based only on
#' \code{.Machine$double.eps^(1/3)}, so ODE truncation does not dominate
#' finite-difference gradient noise.
#' @noRd
#' @keywords internal
.pfimDeSolveTolerances = function( odeSolverParameters ) {
  p = odeSolverParameters %||% list()
  list( atol = p$atol %||% 1e-8, rtol = p$rtol %||% 1e-8 )
}

#' Drop the cached eval-model for one project (e.g. after model equations change).
#' @noRd
#' @keywords internal
.invalidateEvalModelCache = function( pfimproject ) {
  cacheId = .pfimModelCacheId( pfimproject )
  .pfimEvalModelCache[[ cacheId ]] = NULL
  invisible( NULL )
}

#' Arms configured in one constraint-grid cell (single- or multi-arm design).
#' @noRd
#' @keywords internal
.constraintCellArms = function( cell ) {
  # Multi-arm cells store a list of entries; single-arm cells have one `arm`.
  if ( is.null( cell ) )
    return( list() )
  if ( !is.null( cell$entries ) )
    return( compact( map( cell$entries, "arm" ) ) )
  if ( is.null( cell$arm ) )
    return( list() )
  list( cell$arm )
}

#' Coerce a constraint-grid entry or Arm to a bare \code{Arm}.
#'
#' Older FW outputs (and some RDS caches) stored
#' \code{list(arm = <Arm>, samplingsForFW = ...)} in \code{optimalArms}.
#' Plots and reports expect a bare \code{Arm}.
#' @noRd
#' @keywords internal
.pfimAsArm = function( x ) {
  if ( S7::S7_inherits( x, Arm ) )
    return( x )
  if ( is.list( x ) && !is.null( x$arm ) && S7::S7_inherits( x$arm, Arm ) )
    return( x$arm )
  .pfimStop( paste0(
    "Expected an Arm (or list(arm=...)); got ",
    paste( class( x ), collapse = "/" ), "."
  ) )
}

#' Deep copy of the first arm in a constraint-grid cell (avoid shared mutable state).
#' @noRd
#' @keywords internal
.cloneConstraintArm = function( cell ) {
  .pfimCloneS7( .constraintCellArms( cell )[[ 1L ]] )
}

#' Deep copies of all arms in a constraint-grid cell.
#' @noRd
#' @keywords internal
.cloneConstraintCellArms = function( cell ) {
  map( .constraintCellArms( cell ), .pfimCloneS7 )
}

#' TRUE when the first cell of a grid is a joint (multi-arm) design cell.
#' @noRd
#' @keywords internal
.constraintCellJoint = function( cells ) {
  if ( !length( cells ) )
    return( FALSE )
  length( .constraintCellArms( cells[[ 1L ]] ) ) > 1L
}

#' Map multiplicative weights (grid-cell order) onto sorted optimal arms.
#'
#' Optimal arms are named \code{Arm1}, \code{Arm2}, ... matching \code{weightsIndex}.
#' @noRd
#' @keywords internal
.alignOptimalWeightsToArms = function( optimalArms, weightsIndex, optimalWeights ) {
  wByIdx = stats::setNames( as.numeric( optimalWeights ), as.character( weightsIndex ) )
  vapply( optimalArms, function( arm ) {
    idx = sub( "^Arm", "", prop( arm, "name" ) )
    unname( wByIdx[[ idx ]] )
  }, numeric( 1 ) )
}

#' After joint Mult, keep a single winning protocol (weight 1).
#' @noRd
#' @keywords internal
.pfimShrinkJointMultMixture = function( thinOut ) {
  cells = thinOut$listArms
  if ( !length( cells ) || !.constraintCellJoint( cells ) )
    return( thinOut )
  w = as.numeric( thinOut$optimalWeights )
  if ( !length( w ) )
    return( thinOut )
  if ( length( w ) == 1L ) {
    thinOut$optimalWeights = 1
    return( thinOut )
  }
  best = which.max( w )
  thinOut$optimalWeights = 1
  thinOut$weightsIndex   = thinOut$weightsIndex[ best ]
  thinOut$listArms       = cells[ best ]
  thinOut
}

.pfimOptimizerLabel = function( algo ) {
  if ( is.null( algo ) ) return( "NULL" )
  cls = S7::S7_class( algo )
  if ( is.null( cls ) ) class( algo )[[ 1L ]] else cls@name
}

#' Multiplicative outputs nested on the algorithm object (canonical store).
#' @noRd
#' @keywords internal
.pfimMultAlgorithmOutputs = function( optimization ) {
  algo = .getOptimizationAlgorithm( optimization )
  if ( is.null( algo ) || !S7::S7_inherits( algo, MultiplicativeAlgorithm ) )
    .pfimStop(
      paste0(
        "Expected MultiplicativeAlgorithm results; got ",
        .pfimOptimizerLabel( algo ), "."
      )
    )
  out = prop( algo, "multiplicativeAlgorithmOutputs" )
  if ( !length( out ) )
    .pfimStop( "MultiplicativeAlgorithm outputs are empty; run optimizeDesign() first." )
  out
}

#' Bar-plot data for Mult weights / FW frequencies (protocol simplex).
#'
#' When \code{optimalArms} was expanded (joint multi-arm cell), labels are
#' \code{Protocol1}, \code{Protocol2}, ... and one bar per mixture weight.
#' @noRd
#' @keywords internal
.pfimDiscreteMixturePlotData = function( optimization ) {
  out = .pfimAlgoOutputs( optimization )
  w   = out$optimalWeights
  if ( is.null( w ) )
    w = out$frequencies
  w = as.numeric( w )
  if ( !length( w ) )
    .pfimStop( "No mixture weights or frequencies to plot." )
  arms = map( out$optimalArms %||% list(), .pfimAsArm )
  labels = if ( length( arms ) == length( w ) )
    map_chr( arms, ~ prop( .x, "name" ) )
  else
    paste0( "Protocol", seq_along( w ) )
  list(
    data  = data.frame( label = labels, value = w, stringsAsFactors = FALSE ),
    xlab  = if ( length( arms ) == length( w ) ) "Arm" else "Protocol"
  )
}

#' Require a specific optimizer class on an Optimization result.
#' @noRd
#' @keywords internal
.pfimRequireOptimizerClass = function( optimization, algoClass, fn ) {
  algo = .getOptimizationAlgorithm( optimization )
  if ( is.null( algo ) || !S7::S7_inherits( algo, algoClass ) ) {
    .pfimStop( paste0(
      fn, "() requires ", algoClass@name, " results; got ",
      .pfimOptimizerLabel( algo ), "."
    ) )
  }
  algo
}

#' Console summary of discrete mixture weights and proportional N.
#'
#' Only prints for MultiplicativeAlgorithm when more than one weight is active.
#' @noRd
#' @keywords internal
.pfimShowOptimalMixtureWeights = function( optimization, armsDataDf ) {
  algoOut = .pfimAlgoOutputs( optimization )
  algo    = algoOut$optimizationAlgorithm
  if ( !S7::S7_inherits( algo, MultiplicativeAlgorithm ) ) return( invisible( NULL ) )

  ma = prop( algo, "multiplicativeAlgorithmOutputs" )
  if ( length( ma$weightsIndex ) < 2L ) return( invisible( NULL ) )

  # Prefer the algorithm's total N; otherwise sum unique arm sizes from the table.
  N_total = {
    n = ma$numberOfSubjects
    if ( is.null( n ) )
      n = ma$numberOfArms
    if ( length( n ) == 1L && is.finite( n ) && n > 0 ) {
      as.double( n )
    } else {
      sum( tapply(
        as.numeric( armsDataDf[[ "Number of subjects" ]] ),
        armsDataDf[[ "Arms name" ]],
        function( x ) x[[ 1L ]]
      ) )
    }
  }
  nAlloc  = .allocProportionalSubjects( N_total, ma$optimalWeights )

  cat( "\n--- Optimal mixture weights ---\n\n" )
  wtDf = data.frame(
    `Grid cell`  = ma$weightsIndex,
    Weight       = round( ma$optimalWeights, 4 ),
    `N subjects` = nAlloc,
    check.names  = FALSE
  )
  print( wtDf, row.names = FALSE )
  cat( sprintf(
    "\n  (N subjects = Hamilton / largest-remainder of %g * weight; sum(N) = %g)\n",
    N_total, sum( nAlloc )
  ) )

  invisible( NULL )
}

#' Evaluate one (dose x sampling) cell of the discrete constraint grid.
#'
#' Clones arms, applies the dose index and sampling combination, runs the FIM
#' evaluation (optionally reusing a pre-built model), and packs the lower
#' triangle for Fedorov-Wynn.
#' @noRd
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

  # Fresh arm copies per grid cell (design template arms are reused across cells).
  arms = map( arms, .pfimCloneS7 )

  # Assign the dose for this dose-stratum index to each administration.
  armsWithDoses = map( arms, function( arm ) {
    armName         = prop( arm, "name" )
    administrations = map( prop( arm, "administrations" ), function( adm ) {
      adm = .pfimCloneS7( adm )
      prop( adm, "dose" ) = dosesForDesign[[ armName ]][[ prop( adm, "outcome" ) ]][ iterDose ]
      adm
    } )
    prop( arm, "administrations" ) = administrations
    arm
  } )

  # Assign sampling times from the combination grid column for this arm.
  armsUpdated = map( armsWithDoses, function( arm ) {
    armName       = prop( arm, "name" )
    idx           = combinationGrid[ iterComb, armName ]
    samplingEntry = pluck( samplingsForFIMs, designName, armName, idx )
    prop( arm, "samplingTimes" ) = map( samplingEntry, .pfimCloneS7 )
    list(
      arm            = .pfimCloneS7( arm ),
      samplingsForFW = unlist( map( samplingEntry, ~ prop( .x, "samplings" ) ), use.names = FALSE )
    )
  } )

  samplingsForFW = unlist( map( armsUpdated, "samplingsForFW" ), use.names = FALSE )
  armResult      = list( entries = armsUpdated, samplingsForFW = samplingsForFW )
  tempDesign     = .pfimCloneS7( design )
  prop( tempDesign, "arms" ) = map( armsUpdated, "arm" )

  # Fast path: reuse a prepared model/FIM; evaluateArm clones model per arm.
  if ( !is.null( baseModel ) && !is.null( baseFim ) ) {
    evaluatedDesign = evaluateDesign( tempDesign, baseModel, baseFim )
    fisherMatrix    = prop( prop( evaluatedDesign, "fim" ), "fisherMatrix" )
    evalResult      = .pfimEvaluationFromDesign( evaluation, tempDesign, evaluatedDesign )
  } else {
    tempEval = .pfimCloneS7( evaluation )
    prop( tempEval, "designs" ) = list( tempDesign )
    evalResult   = .pfimRunEvaluationCached( tempEval )
    fisherMatrix = getFim( evalResult )$fisherMatrix
  }

  dimFim = nrow( fisherMatrix )

  fisherMatrixForAlgoFW = .packFisherLowerTriangle( fisherMatrix )

  if ( show_progress )
    message( sprintf( "FIM evaluation: %d / %d", fimIndex, totalIterations ) )

  list(
    armResult             = armResult,
    samplingsForFW        = samplingsForFW,
    fisherMatrixForAlgoFW = fisherMatrixForAlgoFW,
    fisherMatrix          = fisherMatrix,
    dimFim                = dimFim,
    cachedEvaluation      = evalResult
  )
}

#' Bind a list of row-lists into one data.frame (empty list -> empty frame).
#' @noRd
#' @keywords internal
.as_df_rows = function( lst ) {
  if ( !length( lst ) ) return( as.data.frame( list() ) )
  map( lst, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
}

#' Flatten nested arm-constraint row lists into a single table.
#' @noRd
#' @keywords internal
.constraintsArmsTable = function( armsConstraints )
  map( armsConstraints, .as_df_rows ) |> list_rbind()

#' Named rows describing continuous (window-based) sampling constraints for reports.
#' @noRd
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

#' Build project-level model spec (no arm binding yet).
#'
#' Resolves library equations, picks the model class, optionally attaches an FD
#' Hessian scheme, wraps the RHS, and applies covariate data when needed.
#' @noRd
#' @keywords internal
.buildEvalModelSpec = function( pfimproject, finiteDifference = FALSE ) {
  if ( length( projectProp( pfimproject, "modelFromLibrary" ) ) != 0L ) {
    projectProp( pfimproject, "modelEquations" ) =
      defineModelEquationsFromLibraryOfModel( pfimproject )
  }
  .pfimValidateProjectOutcomes( pfimproject )
  model = defineModelType( pfimproject )
  if ( isTRUE( finiteDifference ) )
    model = finiteDifferenceHessian( model )
  model = defineModelWrapper( model, pfimproject )
  if ( usesCovariateOccasionStructure( model ) )
    model = defineCovariatesData( model )
  ensureModelOutputNames( model, pfimproject )
  model
}

#' Cached project model; arm administration applied later.
#'
#' Keys on model signature + plain/FD mode. Returns a deep clone so callers can
#' mutate without contaminating the cache.
#' @return The updated project object with rebuilt evaluation model.
#' @keywords internal
rebuildEvalModel = function( pfimproject, finiteDifference = FALSE ) {
  cacheId = .pfimModelCacheId( pfimproject )
  if ( is.null( .pfimEvalModelCache[[ cacheId ]] ) )
    .pfimEvalModelCache[[ cacheId ]] = list()
  cacheKey = if ( isTRUE( finiteDifference ) ) "fd" else "plain"
  if ( !is.null( .pfimEvalModelCache[[ cacheId ]][[ cacheKey ]] ) )
    return( .pfimCloneS7( .pfimEvalModelCache[[ cacheId ]][[ cacheKey ]] ) )
  model = .buildEvalModelSpec( pfimproject, finiteDifference = finiteDifference )
  .pfimEvalModelCache[[ cacheId ]][[ cacheKey ]] = model
  .pfimCloneS7( model )
}

#' Ensure \code{outputNames} is filled from project \code{outputs} when empty.
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
