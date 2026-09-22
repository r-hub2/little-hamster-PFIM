#' @title MultiplicativeAlgorithm
#' @description
#' Multiplicative weight algorithm for optimal sampling-time allocation (Rcpp kernel).
#' Pass \code{lambda}, \code{delta}, \code{numberOfIterations}, and
#' \code{weightThreshold} via \code{optimizerParameters} on \code{\link{Optimization}}.
#'
#' Outputs store two D-values: \code{mixtureDcriterion} is the D of the continuous
#' mixture over retained candidate cells; \code{realisedDcriterion} is the D of
#' the single protocol actually implemented. For joint multi-outcome cells (several
#' arms sharing one protocol in the constraint grid), only the best-weighted cell
#' is kept - mixture weights describe the optimisation simplex, not a multi-protocol
#' subject allocation across joint arms.
#' @param multiplicativeAlgorithmOutputs Compact outputs after optimize
#'   (retained \code{listArms}, weights, D-criteria, \code{algorithmOutput}).
#'   Filled by \code{optimizeDesign()}.
#' @return A \code{MultiplicativeAlgorithm} specification object.
#' @examples
#' \dontrun{
#' vignette("Example01")
#' }
#' @include Optimization.R
#' @export

MultiplicativeAlgorithm = new_class("MultiplicativeAlgorithm",
                                    package = "PFIM",

                                    properties = list(
                                      multiplicativeAlgorithmOutputs = new_property(class_list, default = list())
                                    ))
S4_register( MultiplicativeAlgorithm )

plotWeightsMultiplicativeAlgorithm = new_generic( "plotWeightsMultiplicativeAlgorithm", c( "optimization", "optimizationAlgorithm" ) )

#' MultiplicativeAlgorithm_Rcpp
#'
#' Calls the compiled Rcpp implementation of the multiplicative weight-update
#' algorithm. Registered via \code{Rcpp::compileAttributes()}
#' in \code{R/RcppExports.R}.
#'
#' @param fisherMatrices  List of FIM matrices.
#' @param n_fim           Integer number of FIMs.
#' @param weights         Numeric vector of initial weights.
#' @param p               Integer FIM dimension.
#' @param lambda          Numeric exponent parameter lambda.
#' @param delta           Numeric convergence tolerance delta.
#' @param iteration_init  Integer maximum iterations.
#' @param show_process    Logical; print iteration progress to the console.
#' @return Named list: \code{weights}, \code{iterations}, \code{converged},
#'   \code{singularFim}.
#' @name MultiplicativeAlgorithm_Rcpp
#' @keywords internal
NULL

#' Multiplicative discrete D-optimal design (multi-design driver).
#'
#' @param optimizationObject An \code{\link{Optimization}} project.
#' @param optimizationAlgorithm A \code{MultiplicativeAlgorithm} instance.
#' @return The same \code{Optimization} with \code{optimisationDesign} and
#'   \code{optimisationAlgorithmOutputs} filled.
#' @name optimizeDesign
#' @keywords internal

method( optimizeDesign, list( Optimization, MultiplicativeAlgorithm ) ) = function(
    optimizationObject, optimizationAlgorithm ) {
  .pfimDiscreteOptimizeDesigns(
    optimizationObject, optimizationAlgorithm, .optimizeMultiplicativeOneDesign
  )
}

#' @noRd
#' @keywords internal
.optimizeMultiplicativeOneDesign = function( optimizationObject, optimizationAlgorithm ) {

  p                  = projectProp( optimizationObject, "optimizerParameters" )
  lambda             = p$lambda
  delta              = p$delta
  numberOfIterations = p$numberOfIterations
  weightThreshold    = p$weightThreshold
  showProcess        = isTRUE( p$showProcess )

  # Candidate FIMs / arms for every dose x sampling cell under constraints.
  fimsFromConstraints = generateFimsFromConstraints( optimizationObject )
  design              = pluck( projectProp( optimizationObject, "designs" ), 1 )
  initialDesign       = .pfimCloneS7( design )
  optimalDesign       = .pfimCloneS7( design )
  designName          = prop( design, "name" )
  # Total study size (FW uses numberOfSubjects; Design$size / sum of arm sizes).
  numberOfSubjects    = {
    sz = prop( design, "size" )
    if ( length( sz ) == 1L && is.finite( sz ) && sz > 0 )
      as.double( sz )
    else
      sum( map_dbl( prop( design, "arms" ), ~ prop( .x, "size" ) ) )
  }
  armFims             = fimsFromConstraints$listArms[[ designName ]]
  # Dense pxp elementary FIMs (R<->C++ Mult contract). FW uses packed triangles.
  fisherMatrices      = fimsFromConstraints$listFimsAlgoMult[[ designName ]]
  fim                 = projectProp( optimizationObject, "fim" )

  # Individual / Bayesian: D-criterion is convex in covariance-mixture weights ->
  # optimum is a single protocol (vertex). Exhaustive max over candidates.
  if ( .pfimIsSubjectLevelFim( fim ) ) {
    best           = .pfimBestSubjectProtocol( fisherMatrices )
    weightsIndex   = best$index
    optimalWeights = 1
    mixtureDcriterion = best$Dcriterion
    thinOut = list(
      listArms          = armFims[ weightsIndex ],
      algorithmOutput   = list(
        iterations  = 0L,
        converged   = TRUE,
        singularFim = FALSE,
        method      = "bestSubjectProtocol"
      ),
      numberOfSubjects  = numberOfSubjects,
      weightThreshold   = weightThreshold,
      weightsIndex      = weightsIndex,
      optimalWeights    = optimalWeights,
      mixtureDcriterion = mixtureDcriterion
    )
    optimalArms = setOptimalArms( fim, optimizationAlgorithm, out = thinOut )
    prop( optimalDesign, "arms" ) = optimalArms
    thinOut$realisedDcriterion = mixtureDcriterion
    prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" ) = thinOut
    return( .pfimStoreOptimization(
      optimizationObject, initialDesign, optimalDesign,
      algorithmOutputs = list(
        optimizationAlgorithm = optimizationAlgorithm,
        optimalArms           = optimalArms,
        optimalWeights        = optimalWeights,
        mixtureDcriterion     = thinOut$mixtureDcriterion,
        realisedDcriterion    = thinOut$realisedDcriterion
      )
    ) )
  }

  # Uniform start on the simplex; C++ updates w <- w * m^lambda / ||.||_1 in place.
  # Returned weights are the last *certified* simplex point (test before update).
  numberOfFisherMatrices = length( fisherMatrices )
  weights                = rep( 1 / numberOfFisherMatrices, numberOfFisherMatrices )
  fimDim = dim( pluck( fisherMatrices, 1 ) )[ 1L ]

  multiplicativeAlgorithmOutput = MultiplicativeAlgorithm_Rcpp(
    fisherMatrices, numberOfFisherMatrices, weights, fimDim,
    lambda, delta, numberOfIterations, showProcess
  )

  if ( isTRUE( multiplicativeAlgorithmOutput$singularFim ) )
    stop(
      "MultiplicativeAlgorithm: the weighted mixture FIM became singular ",
      "(or weights collapsed). Check candidate protocols / FIMs, or relax ",
      "constraints / weightThreshold.",
      call. = FALSE
    )

  if ( !isTRUE( multiplicativeAlgorithmOutput$converged ) )
    warning(
      "MultiplicativeAlgorithm: did not converge within numberOfIterations = ",
      numberOfIterations, " (delta = ", delta, ").",
      call. = FALSE
    )

  weights        = multiplicativeAlgorithmOutput$weights
  # Drop near-zero weights before building optimal arms (affects realised D).
  weightsIndex   = which( weights > weightThreshold )
  if ( !length( weightsIndex ) )
    stop(
      "MultiplicativeAlgorithm: no protocol weight exceeds weightThreshold = ",
      weightThreshold, ".",
      call. = FALSE
    )
  optimalWeights = weights[ weightsIndex ]
  optimalWeights = optimalWeights / sum( optimalWeights )
  # D-criterion of the continuous mixture (kept instead of full FIM list -> smaller RDS).
    mixtureDcriterion = .fimDcriterionFromMatrix(
      reduce( map2( fisherMatrices, weights, `*` ), `+` )
    )
  # Compact algo status (shared key algorithmOutput with continuous optimizers).
  algoStatus = list(
    iterations  = multiplicativeAlgorithmOutput$iterations,
    converged   = multiplicativeAlgorithmOutput$converged,
    singularFim = multiplicativeAlgorithmOutput$singularFim
  )

  # Retained cells only - enough to re-call setOptimalArms; not the full grid.
  thinOut = list(
    listArms          = armFims[ weightsIndex ],
    algorithmOutput   = algoStatus,
    numberOfSubjects  = numberOfSubjects,
    weightThreshold   = weightThreshold,
    weightsIndex      = weightsIndex,
    optimalWeights    = optimalWeights,
    mixtureDcriterion = mixtureDcriterion
  )

  # Dispatch by FIM type (population: Hamilton N from weights).
  optimalArms = setOptimalArms( fim, optimizationAlgorithm, out = thinOut )
  # Joint multi-outcome cells: mixture weights are not a multi-protocol allocation -
  # keep one winning protocol (weight 1). mixtureDcriterion stays the continuous D;
  # realisedDcriterion below is the implemented single-protocol D.
  if ( length( thinOut$listArms ) && .constraintCellJoint( thinOut$listArms ) ) {
    thinOut        = .pfimShrinkJointMultMixture( thinOut )
    optimalWeights = as.numeric( thinOut$optimalWeights )
  } else if ( length( optimalArms ) == length( optimalWeights ) ) {
    optimalWeights = .alignOptimalWeightsToArms(
      optimalArms, weightsIndex, optimalWeights
    )
    optimalWeights = optimalWeights / sum( optimalWeights )
    thinOut$optimalWeights = optimalWeights
  }
  prop( optimalDesign, "arms" ) = optimalArms

  # Realised D of the allocated design (before store, so it is persisted).
  thinOut$realisedDcriterion = as.numeric(
    Dcriterion( prop(
      .pfimRunOptimizationEvaluation( optimizationObject, optimalDesign, name = "" ),
      "fim"
    ) )
  )
  prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" ) = thinOut

  .pfimStoreOptimization(
    optimizationObject, initialDesign, optimalDesign,
    algorithmOutputs = list(
      optimizationAlgorithm = optimizationAlgorithm,
      optimalArms           = optimalArms,
      optimalWeights        = optimalWeights,
      mixtureDcriterion     = thinOut$mixtureDcriterion,
      realisedDcriterion    = thinOut$realisedDcriterion
    )
  )
}

#' Weight trajectories of the multiplicative algorithm
#' @name plotWeightsMultiplicativeAlgorithm
#' @return A \code{ggplot2} plot object.
#' @keywords internal

method( plotWeightsMultiplicativeAlgorithm, list( Optimization, MultiplicativeAlgorithm ) ) = function( optimization, optimizationAlgorithm )
{
  plotDat = .pfimDiscreteMixturePlotData( optimization )
  .algoBars(
    data.frame(
      optimalArmsName = plotDat$data$label,
      optimalWeights  = plotDat$data$value
    ),
    x = "optimalArmsName", y = "optimalWeights",
    xlab = plotDat$xlab, ylab = "Weights",
    ylim = c( 0, 1 ), breaks = seq( 0, 1, by = 0.1 )
  )
}

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @keywords internal

method( constraintsTableForReport, MultiplicativeAlgorithm ) = function( optimizationAlgorithm, arms ) {
  .pfimConstraintsTableDiscrete( optimizationAlgorithm, arms )
}
