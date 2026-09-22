# Continuous-optimizer helpers (prepare design, seed samplings, D-criterion /
# metaheuristic fitness). Used by Simplex / PSO / PGBO.
#
# Fitness contract shared with the C++ kernels:
#   - PSO / Simplex minimise cost = 1/D (invalid -> .metaheuristicFitnessPenalty).
#   - PGBO's eval_d returns raw D (invalid -> tiny sentinel); C++ minimises 1/D.
# Flat positions always go through .applyFlatToArms before FIM evaluation.

# --- continuous sampling optimizers (PSO / PGBO / Simplex)

#' Validate constraints, attach optimisation windows, seed from user
#' \code{initialSamplings}, return working design + arms.
#' Shared by Simplex / PSO / PGBO so the reported initial design matches the
#' optimizer start.
#' @noRd
#' @keywords internal
.pfimPrepareContinuousDesign = function( optimizationObject ) {
  design = pluck( projectProp( optimizationObject, "designs" ), 1L )
  checkValiditySamplingConstraint( design )
  design = setSamplingConstraintForOptimization( design )
  arms   = .pfimSeedArmsFromInitialSamplings( prop( design, "arms" ) )
  prop( design, "arms" ) = arms
  list(
    design = design,
    arms   = arms
  )
}

#' Seed arm sampling times from user \code{SamplingTimeConstraints$initialSamplings}.
#'
#' Ensures report initial FIM and C++ \code{initial_pos} / Simplex start agree.
#' Exploration RNG (\code{seed}) applies only after this step for PSO/PGBO.
#' @noRd
#' @keywords internal
.pfimSeedArmsFromInitialSamplings = function( arms ) {
  map( arms, function( arm ) {
    arm = .pfimCloneS7( arm )
    constraints = prop( arm, "samplingTimesConstraints" )
    armName = prop( arm, "name" )
    by_outcome = set_names(
      map( constraints, function( sc ) {
        samp = as.numeric( prop( sc, "initialSamplings" ) )
        outcome = prop( sc, "outcome" )
        if ( !length( samp ) )
          stop(
            "Continuous optimizer: initialSamplings is empty for outcome '",
            outcome, "'.",
            call. = FALSE
          )
        n_by_win = prop( sc, "numberOfTimesByWindows" )
        if ( length( n_by_win ) && sum( n_by_win ) != length( samp ) )
          stop(
            "Continuous optimizer: length(initialSamplings) for '", outcome,
            "' (", length( samp ), ") must equal sum(numberOfTimesByWindows) (",
            sum( n_by_win ), ").",
            call. = FALSE
          )
        sort( samp )
      } ),
      map_chr( constraints, ~ prop( .x, "outcome" ) )
    )
    times = map( prop( arm, "samplingTimes" ), function( st ) {
      outcome = prop( st, "outcome" )
      if ( outcome %in% names( by_outcome ) ) {
        samp = by_outcome[[ outcome ]]
        arm_samp = as.numeric( prop( st, "samplings" ) )
        if ( length( arm_samp ) &&
             length( arm_samp ) != length( samp ) )
          stop(
            "Continuous optimizer: initialSamplings length for '", outcome,
            "' (", length( samp ), ") does not match SamplingTimes length (",
            length( arm_samp ), ").",
            call. = FALSE
          )
        if ( length( arm_samp ) &&
             !isTRUE( all.equal(
               sort( unname( arm_samp ) ), sort( unname( samp ) ),
               tolerance = 1e-8, check.attributes = FALSE
             ) ) )
          .pfimWarn(
            "Arm '", armName, "' outcome '", outcome,
            "': samplingTimes differ from SamplingTimeConstraints$initialSamplings; ",
            "the optimizer starts from initialSamplings."
          )
        prop( st, "samplings" ) = samp
      }
      st
    } )
    prop( arm, "samplingTimes" ) = times
    arm
  } )
}

#' Alias kept for older call sites.
#' @noRd
#' @keywords internal
.pfimSeedArmsFromConstraints = function( arms ) {
  .pfimSeedArmsFromInitialSamplings( arms )
}

#' Outcome names per arm (used by Simplex column keys and constraint checks).
#' @noRd
#' @keywords internal
.pfimArmOutcomesMap = function( arms ) {
  map(
    set_names( arms, map_chr( arms, ~ prop( .x, "name" ) ) ),
    ~ map_chr( prop( .x, "samplingTimesConstraints" ), ~ prop( .x, "outcome" ) )
  )
}

#' Optimal design = frozen initial snapshot with optimised arms.
#' @noRd
#' @keywords internal
.pfimOptimalDesignFrom = function( initialDesign, finalArms ) {
  optimal = .pfimCloneS7( initialDesign )
  prop( optimal, "arms" ) = finalArms
  optimal
}

#' D-criterion at a flat sampling position (PGBO fitness; PSO uses cost = 1/D).
#' @noRd
#' @keywords internal
.pfimFlatDesignDcriterion = function( evalTemplate, initialDesign, arms, flat_pos,
                                      ctx = NULL ) {
  if ( is.null( ctx ) )
    ctx = .pfimMetaheuristicEvalContext( evalTemplate, initialDesign, arms )
  d = as.numeric( Dcriterion( ctx( flat_pos ) ) )
  if ( !is.finite( d ) || d <= 0 ) .metaheuristicInvalidDcriterion else d
}

#' Store initial/optimal evaluations and algorithm outputs (shared by all algos).
#'
#' Never aliases the two design objects. Flushes the FIM cache before re-eval so
#' report FIMs are not stale from the search loop (metaheuristics hammer the
#' cache with many flat positions; the final pair must be recomputed cleanly).
#' @noRd
#' @keywords internal
.pfimStoreOptimization = function( optimizationObject, initialDesign, optimalDesign,
                                   algorithmOutputs, optimalEvalName = "" ) {
  .pfimOptimizationCacheFlush()
  prop( optimizationObject, "optimisationDesign" ) = list(
    evaluationInitialDesign = .pfimRunOptimizationEvaluation(
      optimizationObject, initialDesign
    ),
    evaluationOptimalDesign = .pfimRunOptimizationEvaluation(
      optimizationObject, optimalDesign, name = optimalEvalName
    )
  )
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = algorithmOutputs
  optimizationObject
}

#' Continuous store: set \code{optimizerOutputs} then shared store.
#' @noRd
#' @keywords internal
.pfimStoreContinuousOptimization = function( optimizationObject, optimizationAlgorithm,
                                             initialDesign, optimalDesign,
                                             optimizerOutputs, finalArms ) {
  prop( optimizationAlgorithm, "optimizerOutputs" ) = optimizerOutputs
  .pfimStoreOptimization(
    optimizationObject, initialDesign, optimalDesign,
    algorithmOutputs = list(
      optimizationAlgorithm = optimizationAlgorithm,
      optimalArms           = finalArms
    )
  )
}

#' Stall-stop \code{converged} flag for PSO/PGBO.
#'
#' When \code{tolerance <= 0}, stall stopping is disabled: return \code{NA}
#' (not applicable) instead of \code{FALSE}, so reports do not look like failure.
#' @noRd
#' @keywords internal
.pfimStallConvergedFlag = function( converged, tolerance ) {
  if ( ( tolerance %||% 0 ) <= 0 )
    return( NA )
  isTRUE( converged )
}

#' Shared continuous \code{algorithmOutput} blob.
#' @noRd
#' @keywords internal
.pfimContinuousAlgoOutputs = function( algorithmOutput, tolerance ) {
  algorithmOutput$converged = .pfimStallConvergedFlag(
    algorithmOutput$converged, tolerance
  )
  list( algorithmOutput = algorithmOutput )
}

#' Warn when a continuous optimizer had a stall tolerance but did not stop on it.
#' @noRd
#' @keywords internal
.pfimWarnIfNotConverged = function( algoName, tolerance, converged, extra = "" ) {
  if ( ( tolerance %||% 0 ) > 0 && !isTRUE( converged ) )
    warning( algoName, ": did not converge", extra, ".", call. = FALSE )
  invisible( NULL )
}

# --- metaheuristic objective (Simplex / PSO / PGBO)
# cost = 1/D; singular FIM -> large penalty (.metaheuristicFitnessPenalty in Optimization.R)

#' Reusable evaluation context for flat sampling fitness (one eval/design clone).
#' @noRd
#' @keywords internal
.pfimMetaheuristicEvalContext = function( evalTemplate, design, arms ) {
  evalWork   = .pfimCloneS7( evalTemplate )
  designWork = .pfimCloneS7( design )
  function( flat_pos ) {
    prop( designWork, "arms" ) = .applyFlatToArms( flat_pos, arms )
    prop( evalWork, "designs" ) = list( designWork )
    prop( .pfimRunEvaluationCached( evalWork ), "fim" )
  }
}

#' Scalar fitness (cost) at one flat sampling position.
#' @noRd
#' @keywords internal
.pfimMetaheuristicFitness = function( evalTemplate, design, arms, flat_pos,
                                      ctx = NULL ) {
  if ( is.null( ctx ) )
    ctx = .pfimMetaheuristicEvalContext( evalTemplate, design, arms )
  d = 1 / Dcriterion( ctx( flat_pos ) )
  if ( !is.finite( d ) ) .metaheuristicFitnessPenalty else d
}

#' Batch fitness for a matrix of flat positions (dedupes identical rows).
#' @noRd
#' @keywords internal
.pfimMetaheuristicFitnessBatch = function( evalTemplate, design, arms, pos_matrix,
                                           layout = NULL, ctx = NULL ) {
  pos_matrix = as.matrix( pos_matrix )
  n = nrow( pos_matrix )
  if ( n == 0L ) return( numeric( 0 ) )
  if ( is.null( ctx ) )
    ctx = .pfimMetaheuristicEvalContext( evalTemplate, design, arms )

  # Single validity gate when layout is supplied (scalar and batch share this).
  eval_flat = function( flat ) {
    if ( !is.null( layout ) && !.isFlatValid( layout, flat ) )
      return( .metaheuristicFitnessPenalty )
    .pfimMetaheuristicFitness( evalTemplate, design, arms, flat, ctx = ctx )
  }

  if ( n == 1L )
    return( eval_flat( pos_matrix[ 1L, , drop = TRUE ] ) )

  # Dedupe identical rows before calling run() (swarm / population batches).
  row_key = apply( pos_matrix, 1L, paste, collapse = "\r" )
  ukeys = unique( row_key )
  uniq_idx = match( ukeys, row_key )

  uniq_costs = map(
    uniq_idx,
    function( i ) eval_flat( pos_matrix[ i, , drop = TRUE ] )
  ) |> unlist( use.names = FALSE )
  stats::setNames( uniq_costs, ukeys )[ row_key ] |> unname()
}

#' Run a list of Evaluations through the design FIM cache.
#' @noRd
#' @keywords internal
.pfimRunEvaluations = function( evaluations ) {
  if ( !length( evaluations ) ) return( evaluations )
  map( seq_along( evaluations ), function( i ) .pfimRunEvaluationCached( evaluations[[ i ]] ) )
}

#' Cached Evaluation of one design under an Optimization's model settings.
#' @noRd
#' @keywords internal
.pfimRunOptimizationEvaluation = function( optimization, design, name = "" ) {
  .pfimRunEvaluationCached(
    .evaluationFromOptimization( optimization, design, name = name )
  )
}
