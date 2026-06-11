#' @title PGBOAlgorithm
#' @description
#' Population-based global optimization for sampling times (Rcpp kernel).
#' @param N Numeric: the parameter N.
#' @param muteEffect Numeric: the parameter muteEffect.
#' @param maxIteration Numeric: the parameter maxIteration.
#' @param purgeIteration Numeric: the parameter purgeIteration.
#' @param seed Numeric: the parameter seed.
#' @param showProcess Logical: showProcess.
#' @param fitBase Numeric: baseline fitness for PGBO (default 0.03).
#' @include Optimization.R
#' @export

PGBOAlgorithm = new_class( "PGBOAlgorithm", package = "PFIM",

                           properties = list( N = new_property(class_double, default = numeric(0)),
                                              muteEffect = new_property(class_double, default = numeric(0)),
                                              maxIteration = new_property(class_double, default = numeric(0)),
                                              purgeIteration = new_property(class_double, default = numeric(0)),
                                              fitBase = new_property(class_double, default = 0.03),
                                              seed = new_property(class_double, default = numeric(0)),
                                              showProcess = new_property(class_logical, default = FALSE )))
S4_register( PGBOAlgorithm )

#' Population Genetics Based Optimization kernel (Rcpp).
#'
#' Compiled implementation in \code{src/Pgboalgorithm.cpp}. FIM evaluation and
#' per-group constraint checks are R callbacks.
#'
#' @name pgbo_optimize_Rcpp
#' @keywords internal
NULL

#' Search for a D-optimal sampling design
#' @name optimizeDesign
#' @export

method( optimizeDesign, list( Optimization, PGBOAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {

  .pfimFimCacheBegin( optimizationObject )
  design = pluck( projectProp( optimizationObject, "designs" ), 1L )
  checkValiditySamplingConstraint( design )
  design = setSamplingConstraintForOptimization( design )

  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  set.seed( optimizerParameters$seed )

  arms     = prop( design, "arms" )
  outcomes = map( set_names( arms, map_chr( arms, ~ prop( .x, 'name' ) ) ), ~ map_chr( prop( .x, "samplingTimesConstraints" ), ~ prop( .x, "outcome") ) )

  initialArms = map( arms, function( arm ) {
    samplingConstraints = prop( arm, "samplingTimesConstraints" )

    newSamplings = map( samplingConstraints, ~ generateSamplingsFromSamplingConstraints( .x ) ) |>
      set_names( map_chr( samplingConstraints, ~ prop( .x, "outcome" ) ) )

    updatedTimes = map( prop( arm, "samplingTimes" ), function( st ) {
      outcome = prop( st, "outcome" )
      if ( outcome %in% names( newSamplings ) ) prop( st, "samplings" ) = newSamplings[[ outcome ]]
      st
    })

    prop( arm, "samplingTimes" ) = updatedTimes
    arm
  })

  prop( design, "arms" ) = initialArms
  initialDesign = design

  layout = .buildFlatSamplingLayout( design )

  evalObjTemplate = .evaluationFromOptimization( optimizationObject, design, name = "" )

  flat_from_design = function( tgt_design ) {
    .buildFlatSamplingLayout( tgt_design )$initial_flat
  }

  eval_d = function( flat_pos ) {
    .pfimMetaheuristicFitness( evalObjTemplate, design, arms, flat_pos )
  }

  check_valid_group = function( flat_pos, group_id ) {
    .checkFlatValidGroup( layout, flat_pos, group_id )
  }

  fitBase = prop( optimizationAlgorithm, "fitBase" )

  res_cpp = pgbo_optimize_Rcpp(
    initial_pos        = flat_from_design( initialDesign ),
    sorting_groups     = layout$sorting_groups,
    max_iteration      = as.integer( optimizerParameters$maxIteration ),
    N                  = as.integer( optimizerParameters$N ),
    mute_effect        = optimizerParameters$muteEffect,
    purge_iteration    = as.integer( optimizerParameters$purgeIteration ),
    fit_base           = fitBase,
    max_attempts       = 10000L,
    show_process       = optimizerParameters$showProcess,
    eval_d             = eval_d,
    check_valid_group  = check_valid_group
  )

  best_flat = res_cpp$bestDesign

  finalArms = .applyFlatToArms( best_flat, arms )

  optimalDesign = design
  prop( optimalDesign, "arms" ) = finalArms

  prop( optimizationObject, "optimisationDesign" ) = list(
    evaluationInitialDesign = .pfimRunOptimizationEvaluation( optimizationObject, initialDesign ),
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

method( constraintsTableForReport, PGBOAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = .constraintsArmsTable( armsConstraints )
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Samplings windows", "Number of times by windows","Min sampling" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) |>
    kable_styling( bootstrap_options = c( "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}
