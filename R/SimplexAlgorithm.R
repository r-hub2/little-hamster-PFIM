#' @title SimplexAlgorithm
#' @description
#' Nelder-Mead (amoeba) optimization for sampling times (Rcpp kernel).
#' @param pctInitialSimplexBuilding Numeric: the pctInitialSimplexBuilding.
#' @param maxIteration  Numeric: the maxIteration.
#' @param tolerance  Numeric: the tolerance.
#' @param seed  Numeric: the seed.
#' @param showProcess Logical: the showProcess.
#' @include Optimization.R
#' @export

SimplexAlgorithm = new_class( "SimplexAlgorithm", package = "PFIM",
                              properties = list( pctInitialSimplexBuilding = new_property(class_double, default = numeric(0)),
                                                 maxIteration = new_property(class_double, default = numeric(0)),
                                                 seed = new_property(class_double, default = numeric(0)),
                                                 tolerance = new_property(class_double, default = numeric(0)),
                                                 showProcess = new_property(class_logical, default = FALSE ) ) )
S4_register( SimplexAlgorithm )

fisherSimplex = new_generic( "fisherSimplex", c( "optimizationObject" ) )

#' Nelder-Mead simplex minimiser (Rcpp).
#'
#' Compiled implementation in \code{src/Simplexalgorithm.cpp}.
#' The objective \code{funk} remains an R callback
#' (\code{fisherSimplex}) because FIM evaluation stays in R.
#'
#' @param p Numeric matrix of simplex vertices (rows).
#' @param y Numeric vector of objective values at vertices.
#' @param ftol Relative tolerance on the spread of \code{y}.
#' @param itmax Maximum number of iterations.
#' @param funk R function \code{funk(data, pr, outcomes)}.
#' @param outcomes Outcome structure passed to \code{funk}.
#' @param data Optimization object passed to \code{funk}.
#' @param show_process Logical; print progress.
#' @return List with \code{p}, \code{y}, \code{iter}, \code{converge}, \code{results}.
#' @name fun_amoeba_Rcpp
#' @keywords internal
NULL

#' Nelder-Mead simplex minimiser (R interface).
#'
#' @param p parameter p
#' @param y parameter y
#' @param ftol parameter ftol
#' @param itmax parameter itmax
#' @param funk parameter funk
#' @param outcomes The model outcomes.
#' @param data parameter data
#' @param showProcess Boolean.
#' @return List with simplex results from the compiled kernel.
#' @keywords internal

fun_amoeba = function( p, y, ftol, itmax, funk, outcomes, data, showProcess ) {
  fun_amoeba_Rcpp( p, y, ftol, itmax, funk, data, outcomes, showProcess )
}

#' Simplex objective for continuous sampling optimization
#' @name fisherSimplex
#' @export

method( fisherSimplex, Optimization ) = function( optimizationObject, simplex, outcomes )
{
  designs = projectProp( optimizationObject, "designs" )

  # Hash-index simplex positions by outcome key (O(1) lookup per outcome).
  simplex_hash = split(unname(simplex), names(simplex))

  dCrit_list = map_dbl(designs, function(design) {

    designName = prop( design, "name" )
    arms = prop( design, "arms" )

    arms_processed = map(arms, function(arm) {
      armName = prop( arm, "name" )
      samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
      samplingTimes = prop( arm, "samplingTimes" )

      constraints_met = map_lgl(outcomes[[armName]], function(outcome) {
        namesSamplings = toString( c( designName, armName, outcome ) )

        # O(1) lookup from the hash table
        newSamplings = simplex_hash[[namesSamplings]]
        newSamplings = newSamplings[!is.na(newSamplings)]

        samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop( .x,"outcome" ) == outcome ) |> pluck(1)

        evaluation_result = checkSamplingTimeConstraintsForMetaheuristic( samplingTimesConstraint, arm, newSamplings, outcome )

        return( all(unlist(evaluation_result, use.names = FALSE)) )
      })

      updated_samplingTimes = map(samplingTimes, function(st) {
        outcome = prop( st, "outcome" )
        if (outcome %in% outcomes[[armName]]) {
          namesSamplings = toString( c( designName, armName, outcome ) )
          # Extraction O(1)
          prop( st, "samplings" ) = simplex_hash[[namesSamplings]]
        }
        return(st)
      })

      prop( arm, "samplingTimes" ) = updated_samplingTimes
      return(list(arm = arm, valid = all(constraints_met)))
    })

    prop( design, "arms" ) = map(arms_processed, "arm")
    constraints_global_valid = all(map_lgl(arms_processed, "valid"))

    if( constraints_global_valid ) {
      evaluationFIM = .pfimRunOptimizationEvaluation( optimizationObject, design, name = "" )
      fim = prop( evaluationFIM, "fim" )
      return( 1 / Dcriterion( fim ) )

    } else {
      return( 1.0 )
    }
  })

  dCrit = dCrit_list[1]
  dCrit[is.infinite(dCrit)] = 1.0

  return( dCrit )
}

#' Search for a D-optimal sampling design
#' @name optimizeDesign
#' @export

method( optimizeDesign, list( Optimization, SimplexAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm )
{
  .pfimFimCacheBegin( optimizationObject )
  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  showProcess = optimizerParameters$showProcess
  pctInitialSimplexBuilding = optimizerParameters$pctInitialSimplexBuilding
  tolerance = optimizerParameters$tolerance
  maxIteration = optimizerParameters$maxIteration

  designs = projectProp( optimizationObject, "designs" )
  design = pluck( designs, 1 )
  optimalDesign = pluck( designs, 1 )
  designName = prop( design, "name" )

  checkValiditySamplingConstraint( design )
  design = setSamplingConstraintForOptimization( design )
  arms = prop( design, "arms" )
  projectProp( optimizationObject, "designs" ) <- list( design )

  outcomes = map( set_names( arms, map_chr( arms, ~ prop( .x, 'name' ) ) ), ~ {
    map_chr( prop( .x, "samplingTimesConstraints" ), ~ prop( .x, "outcome"))
  })

  simplex_data_list = map(arms, function(arm) {
    armName = prop( arm, "name" )
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

    map(outcomes[[armName]], function(outcome) {
      samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop(.x,"outcome") == outcome ) |> pluck(1)
      samplings = prop( samplingTimesConstraint, "initialSamplings" )

      data.frame(
        name = rep( toString( c( designName, armName, outcome ) ), length( samplings ) ),
        value = samplings,
        stringsAsFactors = FALSE
      )
    }) |> list_rbind()
  })
  simplex_data = list_rbind(simplex_data_list)

  vectorizedSamplings = simplex_data$value
  namesSamplingsSimplex = simplex_data$name
  n_vars = length(vectorizedSamplings)

  samplingsSimplex = matrix( rep( vectorizedSamplings, n_vars + 1 ), ncol = n_vars, byrow = TRUE )
  colnames( samplingsSimplex ) = namesSamplingsSimplex
  samplingsSimplex = t( apply( samplingsSimplex, 1, sort ) )

  diag_indices = matrix(c(2:(n_vars + 1), 1:n_vars), ncol = 2)
  samplingsSimplex[diag_indices] = samplingsSimplex[diag_indices] * ( 1 - pctInitialSimplexBuilding / 100 )

  # Extraction de la fonction pure en amont pour Ã©viter l'Ã©valuation de signature Ã  chaque cycle
  pure_fisher_func = method(fisherSimplex, Optimization)

  y = map_dbl(seq_len(nrow(samplingsSimplex)), ~ pure_fisher_func( optimizationObject, samplingsSimplex[.x, ], outcomes ))

  # Injection de la fonction pure directement dans l'algorithme
  opti = fun_amoeba_Rcpp( samplingsSimplex, y, tolerance, maxIteration, pure_fisher_func, optimizationObject, outcomes, showProcess )

  indexOptimalDCriteria = which( opti$y == min( opti$y ) )[1]
  optimalsamplingTimes = opti$p[indexOptimalDCriteria, ]

  # Map optimized flat positions back onto arm sampling times
  optimalsamplingTimes_split = split(unname(optimalsamplingTimes), names(optimalsamplingTimes))

  armsList = map(prop(optimalDesign, "arms"), function(arm) {
    armName = prop( arm, "name" )

    listOfSamplingTimes = map(outcomes[[armName]], function(outcome) {
      namesSamplings = toString(c( designName, armName, outcome ))
      samplings = optimalsamplingTimes_split[[namesSamplings]]
      samplings = samplings[!is.na(samplings)]
      SamplingTimes( outcome, samplings = samplings )
    })

    prop( arm, "samplingTimes" ) = listOfSamplingTimes
    return(arm)
  })

  prop( optimalDesign, "arms" ) = armsList

  evaluationInitialDesign = .pfimRunOptimizationEvaluation( optimizationObject, design )
  evaluationOptimalDesign = .pfimRunOptimizationEvaluation( optimizationObject, optimalDesign )

  prop( optimizationObject, "optimisationDesign" ) = list( evaluationInitialDesign = evaluationInitialDesign, evaluationOptimalDesign = evaluationOptimalDesign )
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list( "optimizationAlgorithm" = optimizationAlgorithm, "optimalArms" = armsList )

  optimizationObject
}

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @export

method( constraintsTableForReport, SimplexAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = .constraintsArmsTable( armsConstraints )
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Samplings windows", "Number of times by windows","Min sampling" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) |>
    kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}
