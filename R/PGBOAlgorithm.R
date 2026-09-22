#' @title PGBOAlgorithm
#' @description
#' Population-based global optimization for sampling times (Rcpp kernel).
#' Pass \code{N}, \code{muteEffect}, \code{maxIteration}, \code{purgeIteration}, and
#' \code{seed} via \code{optimizerParameters} on \code{\link{Optimization}}.
#' \code{muteEffect} is an additive mutation scale in the same time unit as the
#' sampling windows (e.g. hours), not a ratio in \eqn{[0,1]}.
#' Search starts from each outcome's \code{initialSamplings}; \code{seed} controls
#' mutation draws only.
#' Optional \code{fitBase} (default \code{0.03}) and \code{cauchyProb} (default \code{0.8},
#' probability of a Cauchy jump vs Gaussian group mutation).
#' Optional \code{tolerance} (default \code{0}, stall stop disabled ->
#' \code{converged = NA}) and \code{stallIterations}
#' (default \code{5}) enable early stopping when \code{tolerance > 0}.
#' @param optimizerOutputs List filled by \code{optimizeDesign()} (optimal arms, etc.).
#' @return A \code{PGBOAlgorithm} specification object.
#' @examples
#' \dontrun{
#' vignette("Example02")
#' }
#' @include Optimization.R
#' @export

PGBOAlgorithm = new_class( "PGBOAlgorithm", package = "PFIM",
                           properties = list(
                             optimizerOutputs = new_property( class_list, default = list() )
                           ) )
S4_register( PGBOAlgorithm )

#' Population Genetics Based Optimization kernel (Rcpp).
#'
#' Compiled implementation in \code{src/PGBOAlgorithm.cpp}. R<->C++ contract:
#' \code{sorting_groups} are 1-based; mutations are per group.
#' \code{check_valid_group(trial, group_1based)} is the feasibility gate
#' (no \code{windows_list} clamp - unlike PSO). \code{eval_d} returns D
#' (not 1/D); the kernel minimises \code{1/D}.
#'
#' @name pgbo_optimize_Rcpp
#' @return A list of optimization results.
#' @keywords internal
NULL

#' PGBO continuous D-optimal design (delegates to multi-design driver).
#'
#' @param optimizationObject An \code{\link{Optimization}} project.
#' @param optimizationAlgorithm A \code{PGBOAlgorithm} instance.
#' @return Updated \code{Optimization} after optimizing each design in turn.
#' @name optimizeDesign
#' @keywords internal

method( optimizeDesign, list( Optimization, PGBOAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {
  .pfimContinuousOptimizeDesigns(
    optimizationObject,
    optimizationAlgorithm,
    .optimizePGBOOneDesign
  )
}

#' Run PGBO on a single design: flat sampling vector, C++ population genetics loop.
#'
#' Sets the RNG from \code{optimizerParameters$seed} and restores
#' \code{.Random.seed} on exit (\code{NULL} seed leaves the stream unchanged).
#' \code{eval_d} returns the D-criterion (invalid designs map to a tiny sentinel
#' inside the flat-design helper); the C++ kernel minimises \code{1/D}.
#' @param optimizationObject Parent \code{Optimization}.
#' @param optimizationAlgorithm \code{PGBOAlgorithm} instance.
#' @return List with optimal arms / design pieces for the continuous driver.
#' @noRd
#' @keywords internal
.optimizePGBOOneDesign = function( optimizationObject, optimizationAlgorithm ) {

  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  restore_seed = .pfimLocalSeed( optimizerParameters$seed )
  on.exit( restore_seed(), add = TRUE )

  prep = .pfimPrepareContinuousDesign( optimizationObject )
  design = prep$design
  arms   = prep$arms

  initialDesign = .pfimCloneS7( design )

  layout       = .buildFlatSamplingLayout( design )
  evalTemplate = .evaluationFromOptimization( optimizationObject, initialDesign, name = "" )
  evalCtx      = .pfimMetaheuristicEvalContext( evalTemplate, initialDesign, arms )

  eval_d = function( flat_pos ) {
    .pfimFlatDesignDcriterion( evalTemplate, initialDesign, arms, flat_pos, ctx = evalCtx )
  }

  cauchyProb = optimizerParameters$cauchyProb %||% 0.8

  # muteEffect is in the same time unit as the windows (hours, etc.), not a fraction.
  # check_valid_group: C++ passes 1-based group_id matching sorting_groups order.
  res_cpp = pgbo_optimize_Rcpp(
    initial_pos       = layout$initial_flat,
    sorting_groups    = layout$sorting_groups,
    max_iteration     = as.integer( optimizerParameters$maxIteration ),
    N                 = as.integer( optimizerParameters$N ),
    mute_effect       = optimizerParameters$muteEffect,
    purge_iteration   = as.integer( optimizerParameters$purgeIteration ),
    fit_base          = optimizerParameters$fitBase %||% 0.03,
    max_attempts      = 10000L,
    cauchy_prob       = cauchyProb,
    show_process      = optimizerParameters$showProcess,
    eval_d            = eval_d,
    check_valid_group = function( flat_pos, group_id ) {
      .checkFlatValidGroup( layout, flat_pos, group_id )
    },
    ftol             = optimizerParameters$tolerance %||% 0,
    stall_iterations = as.integer( optimizerParameters$stallIterations %||% 5L )
  )

  finalArms     = .applyFlatToArms( res_cpp$bestDesign, arms )
  optimalDesign = .pfimOptimalDesignFrom( initialDesign, finalArms )

  tol = optimizerParameters$tolerance %||% 0
  .pfimWarnIfNotConverged(
    "PGBOAlgorithm", tol, res_cpp$converged,
    extra = paste0(
      " within maxIteration = ", optimizerParameters$maxIteration,
      " (tolerance = ", tol, ")"
    )
  )
  algoOut = .pfimContinuousAlgoOutputs(
    list(
      bestD      = res_cpp$bestD,
      initialD   = res_cpp$initialD,
      converged  = res_cpp$converged,
      iterations = res_cpp$iterations,
      improved   = res_cpp$improved
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

method( constraintsTableForReport, PGBOAlgorithm ) = function( optimizationAlgorithm, arms ) {
  .pfimConstraintsTableContinuous( optimizationAlgorithm, arms )
}
