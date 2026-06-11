#' @title PSOAlgorithm
#' @description
#' Particle Swarm Optimization for sampling-time search (Rcpp kernel in \code{src/Psokernel.cpp}).
#' @param maxIteration Numeric: the maxIteration.
#' @param populationSize Numeric: the populationSize.
#' @param seed Numeric: the seed.
#' @param personalLearningCoefficient Numeric: the personalLearningCoefficient.
#' @param globalLearningCoefficient Numeric: the globalLearningCoefficient.
#' @param showProcess Logical: the showProcess.
#' @include Optimization.R
#' @export

PSOAlgorithm = new_class( "PSOAlgorithm", package = "PFIM",

                          properties = list( maxIteration = new_property(class_double, default = numeric(0)),
                                             populationSize = new_property(class_double, default = numeric(0)),
                                             seed = new_property(class_double, default = numeric(0)),
                                             personalLearningCoefficient = new_property(class_double, default = numeric(0)),
                                             globalLearningCoefficient = new_property(class_double, default = numeric(0)),
                                             showProcess = new_property(class_logical, default = FALSE ) ) )
S4_register( PSOAlgorithm )

#' Particle Swarm Optimization kernel (Rcpp).
#'
#' Compiled implementation in \code{src/Psokernel.cpp}. FIM evaluation uses
#' \code{eval_fitness} (scalar fallback) and \code{eval_fitness_batch}
#' (one matrix per swarm step).
#'
#' @name pso_optimize_Rcpp
#' @keywords internal
NULL

#' Search for a D-optimal sampling design
#' @name optimizeDesign
#' @export

method( optimizeDesign, list( Optimization, PSOAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {

  # Cache scope tied to this optimization object for the whole PSO run.
  .pfimFimCacheBegin( optimizationObject )
  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  set.seed( optimizerParameters$seed )

  design = pluck( projectProp( optimizationObject, "designs" ), 1L )
  checkValiditySamplingConstraint( design )
  design = setSamplingConstraintForOptimization( design )

  arms = prop( design, "arms" )

  initialArms = map( arms, function( arm ) {
    samplingConstraints = prop( arm, "samplingTimesConstraints" )
    newSamplings = map( samplingConstraints, ~ generateSamplingsFromSamplingConstraints( .x ) ) |>
      set_names( map_chr( samplingConstraints, ~ prop( .x, "outcome" ) ) )
    updatedTimes = map( prop( arm, "samplingTimes" ), function( st ) {
      outcome = prop( st, "outcome" )
      if ( outcome %in% names( newSamplings ) )
        prop( st, "samplings" ) = newSamplings[[ outcome ]]
      st
    })
    prop( arm, "samplingTimes" ) = updatedTimes
    arm
  })
  prop( design, "arms" ) = initialArms

  layout = .buildFlatSamplingLayout( design )

  evalObjTemplate = .evaluationFromOptimization( optimizationObject, design, name = "" )

  eval_fitness = function( flat_pos ) {
    if ( !.isFlatValid( layout, flat_pos ) )
      return( .metaheuristicFitnessPenalty )
    .pfimMetaheuristicFitness( evalObjTemplate, design, arms, flat_pos )
  }

  # Batch path: C++ passes the full swarm matrix after each velocity update.
  eval_fitness_batch = function( pos_matrix ) {
    pos_matrix = as.matrix( pos_matrix )
    n = nrow( pos_matrix )
    if ( n == 0L ) return( numeric( 0 ) )
    costs = .pfimMetaheuristicFitnessBatch(
      evalObjTemplate, design, arms, pos_matrix, layout = layout
    )
    invalid = !vapply(
      seq_len( n ),
      function( i ) .isFlatValid( layout, pos_matrix[ i, ] ),
      logical( 1L )
    )
    costs[ invalid ] = .metaheuristicFitnessPenalty
    costs
  }

  sample_valid_pos = function() .sampleFlatFromConstraints( layout )

  phi1 = optimizerParameters$personalLearningCoefficient
  phi2 = optimizerParameters$globalLearningCoefficient
  phi  = phi1 + phi2

  if ( phi <= 4 )
    stop( sprintf(
      paste0( "PSOAlgorithm: personalLearningCoefficient + globalLearningCoefficient must be > 4\n",
              "for the constriction factor to guarantee convergence (Clerc & Kennedy 2002).\n",
              "Current value: %.4f. Typical default: phi1 = phi2 = 2.05 (phi = 4.1)." ),
      phi
    ), call. = FALSE )

  phi_inner = phi * ( phi - 4 )
  if ( phi_inner <= 0 )
    stop(
      "PSOAlgorithm: phi * (phi - 4) must be positive for the constriction factor.",
      call. = FALSE
    )

  # Clerc & Kennedy constriction; phi1 + phi2 must exceed 4 (checked above).
  constrictionFactor = 2 / abs( 2 - phi - sqrt( phi_inner ) )

  res_cpp = pso_optimize_Rcpp(
    n_pop_in         = as.integer( optimizerParameters$populationSize ),
    max_iter         = as.integer( optimizerParameters$maxIteration ),
    initial_pos      = layout$initial_flat,
    windows_list     = layout$windows_list,
    sorting_groups   = layout$sorting_groups,
    phi1             = phi1,
    phi2             = phi2,
    constriction     = constrictionFactor,
    show_process     = optimizerParameters$showProcess,
    eval_fitness       = eval_fitness,
    eval_fitness_batch = eval_fitness_batch,
    sample_valid_pos   = sample_valid_pos
  )

  best_flat_pos = res_cpp$globalBestDesign

  finalArms = .applyFlatToArms( best_flat_pos, arms )

  optimalDesign = design
  prop( optimalDesign, "arms" ) = finalArms

  prop( optimizationObject, "optimisationDesign" ) = list(
    evaluationInitialDesign = .pfimRunOptimizationEvaluation( optimizationObject, design ),
    evaluationOptimalDesign = .pfimRunOptimizationEvaluation( optimizationObject, optimalDesign )
  )
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list(
    "optimizationAlgorithm" = optimizationAlgorithm,
    "optimalArms"           = finalArms
  )

  optimizationObject
}

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @export

method( constraintsTableForReport, PSOAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = .constraintsArmsTable( armsConstraints )
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Samplings windows", "Number of times by windows","Min sampling" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) |>
    kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}
