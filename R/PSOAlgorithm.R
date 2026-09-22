#' @title PSOAlgorithm
#' @description
#' Particle Swarm Optimization for sampling-time search (Rcpp kernel in \code{src/PSOKernel.cpp}).
#' Pass \code{maxIteration}, \code{populationSize}, \code{seed}, and learning coefficients
#' via \code{optimizerParameters} on \code{\link{Optimization}}.
#' Search starts from each outcome's \code{initialSamplings} (reported initial design);
#' \code{seed} controls swarm exploration after that start.
#' Optional \code{tolerance} (default \code{0}, stall stop disabled ->
#' \code{converged = NA}) and \code{stallIterations}
#' (default \code{5}) enable early stopping when \code{tolerance > 0}.
#' @param optimizerOutputs List filled by \code{optimizeDesign()} (optimal arms, etc.).
#' @return A \code{PSOAlgorithm} specification object.
#' @examples
#' \dontrun{
#' vignette("Example02")
#' }
#' @include Optimization.R
#' @export

PSOAlgorithm = new_class( "PSOAlgorithm", package = "PFIM",
                          properties = list(
                            optimizerOutputs = new_property( class_list, default = list() )
                          ) )
S4_register( PSOAlgorithm )

#' Particle Swarm Optimization kernel (Rcpp).
#'
#' Compiled implementation in \code{src/PSOKernel.cpp}. R<->C++ contract:
#' \code{windows_list} length equals \code{length(initial_pos)};
#' \code{sorting_groups} are 1-based index vectors; FIM evaluation uses
#' \code{eval_fitness} (scalar fallback) and \code{eval_fitness_batch}
#' (preferred: one matrix per swarm step). \code{sample_valid_pos} reseeds
#' particles that leave the feasible set.
#'
#' @name pso_optimize_Rcpp
#' @return A list of optimization results.
#' @keywords internal
NULL

#' PSO continuous D-optimal design (delegates to multi-design driver).
#'
#' @param optimizationObject An \code{\link{Optimization}} project.
#' @param optimizationAlgorithm A \code{PSOAlgorithm} instance.
#' @return Updated \code{Optimization} after optimizing each design in turn.
#' @name optimizeDesign
#' @keywords internal

method( optimizeDesign, list( Optimization, PSOAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {
  .pfimContinuousOptimizeDesigns(
    optimizationObject,
    optimizationAlgorithm,
    .optimizePSOOneDesign
  )
}

#' Run PSO on a single design: flatten sampling windows, call C++ swarm, rebuild arms.
#'
#' Sets the RNG from \code{optimizerParameters$seed} and restores
#' \code{.Random.seed} on exit (\code{NULL} seed leaves the stream unchanged).
#' Fitness is \code{1/D} via R callbacks; invalid flat positions receive a large penalty.
#' @param optimizationObject Parent \code{Optimization}.
#' @param optimizationAlgorithm \code{PSOAlgorithm} instance (outputs filled here).
#' @return List with optimal arms / design evaluation pieces for the driver.
#' @noRd
#' @keywords internal
.optimizePSOOneDesign = function( optimizationObject, optimizationAlgorithm ) {

  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  restore_seed = .pfimLocalSeed( optimizerParameters$seed )
  on.exit( restore_seed(), add = TRUE )

  prep = .pfimPrepareContinuousDesign( optimizationObject )
  design = prep$design
  arms   = prep$arms

  initialDesign = .pfimCloneS7( design )

  layout = .buildFlatSamplingLayout( design )
  evalTemplate = .evaluationFromOptimization( optimizationObject, initialDesign, name = "" )
  evalCtx      = .pfimMetaheuristicEvalContext( evalTemplate, design, arms )

  # Single validity gate: .pfimMetaheuristicFitnessBatch(layout=).
  eval_fitness = function( flat_pos ) {
    .pfimMetaheuristicFitnessBatch(
      evalTemplate, design, arms,
      matrix( as.numeric( flat_pos ), nrow = 1L ),
      layout = layout, ctx = evalCtx
    )
  }

  eval_fitness_batch = function( pos_matrix ) {
    .pfimMetaheuristicFitnessBatch(
      evalTemplate, design, arms, as.matrix( pos_matrix ),
      layout = layout, ctx = evalCtx
    )
  }

  phi1 = optimizerParameters$personalLearningCoefficient
  phi2 = optimizerParameters$globalLearningCoefficient
  phi  = phi1 + phi2
  # phi > 4 validated at Optimization() construction; constriction needs phi*(phi-4)>0.
  # Clerc-Kennedy constriction keeps velocities from exploding on the flat vector.
  phi_inner = phi * ( phi - 4 )
  constriction = 2 / abs( 2 - phi - sqrt( phi_inner ) )

  # Pass flat layout + R fitness callbacks; C++ owns the swarm loop only.
  res_cpp = pso_optimize_Rcpp(
    n_pop_in           = as.integer( optimizerParameters$populationSize ),
    max_iter           = as.integer( optimizerParameters$maxIteration ),
    initial_pos        = layout$initial_flat,
    windows_list       = layout$windows_list,
    sorting_groups     = layout$sorting_groups,
    phi1               = phi1,
    phi2               = phi2,
    constriction       = constriction,
    show_process       = optimizerParameters$showProcess,
    eval_fitness       = eval_fitness,
    eval_fitness_batch = eval_fitness_batch,
    sample_valid_pos   = function() .sampleFlatFromConstraints( layout ),
    ftol               = optimizerParameters$tolerance %||% 0,
    stall_iterations   = as.integer( optimizerParameters$stallIterations %||% 5L )
  )

  finalArms     = .applyFlatToArms( res_cpp$globalBestDesign, arms )
  optimalDesign = .pfimOptimalDesignFrom( initialDesign, finalArms )

  tol = optimizerParameters$tolerance %||% 0
  .pfimWarnIfNotConverged(
    "PSOAlgorithm", tol, res_cpp$converged,
    extra = paste0(
      " within maxIteration = ", optimizerParameters$maxIteration,
      " (tolerance = ", tol, ")"
    )
  )
  algoOut = .pfimContinuousAlgoOutputs(
    list(
      globalBestCost = res_cpp$globalBestCost,
      converged      = res_cpp$converged,
      iterations     = res_cpp$iterations,
      improved       = res_cpp$improved
    ),
    tolerance = tol
  )

  .pfimStoreContinuousOptimization(
    optimizationObject, optimizationAlgorithm,
    initialDesign, optimalDesign,
    optimizerOutputs = c( list( optimalArms = finalArms ), algoOut ),
    finalArms = finalArms
  )
}

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @keywords internal

method( constraintsTableForReport, PSOAlgorithm ) = function( optimizationAlgorithm, arms ) {
  .pfimConstraintsTableContinuous( optimizationAlgorithm, arms )
}
