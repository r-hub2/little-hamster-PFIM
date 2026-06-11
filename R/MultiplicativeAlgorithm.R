#' @title MultiplicativeAlgorithm
#' @description
#' Multiplicative weight algorithm for optimal sampling-time allocation (Seurat et al., 2021).
#' @param lambda Step-size parameter for weight updates.
#' @param delta Convergence threshold on weight changes.
#' @param numberOfIterations Maximum number of iterations.
#' @param weightThreshold Minimum weight below which a time point is removed.
#' @param showProcess If TRUE, print FIM enumeration progress and a short summary.
#' @param multiplicativeAlgorithmOutputs Raw algorithm outputs (weights, history).
#' @include Optimization.R
#' @export

MultiplicativeAlgorithm = new_class("MultiplicativeAlgorithm",
                                    package = "PFIM",

                                    properties = list(
                                      lambda = new_property(class_numeric, default = 0.0),
                                      delta = new_property(class_numeric, default = 0.0),
                                      numberOfIterations = new_property(class_numeric, default = 0),
                                      weightThreshold = new_property(class_numeric, default = 0.0),
                                      showProcess = new_property(class_logical, default = FALSE),
                                      multiplicativeAlgorithmOutputs = new_property(class_list, default = list())
                                    ))
S4_register( MultiplicativeAlgorithm )

plotWeightsMultiplicativeAlgorithm = new_generic( "plotWeightsMultiplicativeAlgorithm", c( "optimization", "optimizationAlgorithm" ) )

#' MultiplicativeAlgorithm_Rcpp
#'
#' Calls the compiled Rcpp implementation of the multiplicative weight-update
#' algorithm (Seurat et al., 2021). Registered via \code{Rcpp::compileAttributes()}
#' in \code{R/RcppExports.R}.
#'
#' @param fisherMatrices         List of FIM matrices.
#' @param numberOfFisherMatrices Integer number of FIMs.
#' @param weights                Numeric vector of initial weights.
#' @param numberOfParameters     Integer number of parameters.
#' @param dim                    Integer FIM dimension (p).
#' @param lambda                 Numeric exponent parameter lambda.
#' @param delta                  Numeric convergence tolerance delta.
#' @param iterationInit          Integer maximum iterations.
#' @return Named list: \code{weights} (final weights), \code{iterationEnd}.
#' @name MultiplicativeAlgorithm_Rcpp
#' @keywords internal
NULL

#' Search for a D-optimal sampling design
#' @name optimizeDesign
#' @export

method( optimizeDesign, list( Optimization, MultiplicativeAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {

  # parameters of the optimization algorithm
  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  lambda = optimizerParameters$lambda
  delta = optimizerParameters$delta
  numberOfIterations = optimizerParameters$numberOfIterations
  weightThreshold = optimizerParameters$weightThreshold
  showProcess = isTRUE( optimizerParameters$showProcess )

  # generate the Fims from administration and sampling times constraints
  fimsFromConstraints = generateFimsFromConstraints( optimizationObject )

  # run the multiplicative algorithm
  designs = projectProp( optimizationObject, "designs" )
  design = pluck( designs, 1 )
  optimalDesign = pluck( designs, 1 )

  # list for the evaluation of the optimal design
  evaluationOptimalDesignList = list()
  evaluationInitialDesignList = list()

  # design name
  designName = prop( design, "name" )

  # number of arms in the design
  numberOfArms = prop( design, "numberOfArms" )

  # set arms and fims from the evaluation of the constraints
  armFims = fimsFromConstraints$listArms[[designName]]
  fisherMatrices = fimsFromConstraints$listFimsAlgoMult[[designName]]

  # multiplicative algorithm parameters
  numberOfFisherMatrices = length( fisherMatrices )
  weights = rep( 1/numberOfFisherMatrices, numberOfFisherMatrices )
  dim = dim( pluck( fisherMatrices,1 ) )[1]

  # run the multiplicative algorithm
  multiplicativeAlgorithmOutput = MultiplicativeAlgorithm_Rcpp(
    fisherMatrices, numberOfFisherMatrices, weights, dim,
    lambda, delta, numberOfIterations, showProcess
  )

  if ( showProcess )
    message( sprintf(
      "Multiplicative algorithm: %d weight updates",
      multiplicativeAlgorithmOutput$iterationEnd
    ) )

  #get the optimal weights
  weights = multiplicativeAlgorithmOutput[["weights"]]
  weightsIndex = which( weights > weightThreshold )
  optimalWeights = weights[ weightsIndex ]

  # set the multiplicativeAlgorithmOutputs
  prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" ) = list( armFims = armFims,
                                                                          multiplicativeAlgorithmOutput = multiplicativeAlgorithmOutput,
                                                                          numberOfArms = numberOfArms,
                                                                          weightThreshold = weightThreshold,
                                                                          weightsIndex = weightsIndex,
                                                                          optimalWeights = optimalWeights )

  # set the optimal arms to the optimal design
  fim =  projectProp( optimizationObject, "fim" )
  optimalArms = setOptimalArms( fim, optimizationAlgorithm )

  # set optimal arms
  prop( optimalDesign, "arms" ) = optimalArms

  evaluationsToRun = list(
    .evaluationFromOptimization( optimizationObject, optimalDesign, name = "...." ),
    .evaluationFromOptimization( optimizationObject, design, name = "" )
  )
  evaluationResults = .pfimRunEvaluations( evaluationsToRun )
  evaluationOptimalDesign = evaluationResults[[ 1L ]]
  evaluationInitialDesign = evaluationResults[[ 2L ]]

  # set the results in evaluation
  prop( optimizationObject, "optimisationDesign" ) = list( evaluationInitialDesign = evaluationInitialDesign,
                                                           evaluationOptimalDesign = evaluationOptimalDesign )

  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list( "optimizationAlgorithm" = optimizationAlgorithm,
                                                                     "optimalArms" = optimalArms, optimalWeights = optimalWeights )

  optimizationObject
}

#' Weight trajectories of the multiplicative algorithm
#' @name plotWeightsMultiplicativeAlgorithm
#' @export

method( plotWeightsMultiplicativeAlgorithm, list( Optimization, MultiplicativeAlgorithm ) ) = function( optimization, optimizationAlgorithm )
{
  optimisationAlgorithmOutputs = prop( optimization, "optimisationAlgorithmOutputs" )
  optimalArms = optimisationAlgorithmOutputs$optimalArms
  optimalArmsName = map( optimalArms, ~ prop(.x,"name" ) ) |> unlist()
  optimalWeights = optimisationAlgorithmOutputs$optimalWeights
  optimalArms = data.frame( optimalArmsName, optimalWeights )

  weightPlot = ggplot(optimalArms, aes(x = reorder(optimalArmsName, optimalWeights), y = optimalWeights)) +
    geom_bar(stat = "identity", fill = "gray50") +
    scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, by = 0.1),minor_breaks = seq(0, 1, by = 0.05),expand = c(0, 0)  ) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(  x = "Arms",  y = "Weights" ) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(color = "black", margin = margin(r = 10)),
      axis.text.x = element_text(color = "black", margin = margin(t = 5)),
      axis.text.y = element_text(color = "black", margin = margin(r = 5)),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor.x = element_line(color = "gray95", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      plot.margin = margin(10, 10, 10, 10) )
  return( weightPlot )
}

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @export

method( constraintsTableForReport, MultiplicativeAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = .as_df_rows( pluck( armsConstraints, 1L ) )
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Fixed times", "Number of samplings optimisable","Dose constraints" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) |>
    kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}
