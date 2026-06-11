#' @keywords internal
.fedorovWynnBufferSizes = function( ndimFim, nProtocols ) {
  ndimFim    = as.integer( ndimFim )
  nProtocols = as.integer( nProtocols )
  nMaxPop    = ndimFim * ( ndimFim + 1L ) / 2L + 1L
  nFisher    = ndimFim * ( ndimFim + 1L ) / 2L
  nBuf       = max( nProtocols * 2L, nMaxPop )
  list( nMaxPop = nMaxPop, nFisher = nFisher, nBuf = nBuf )
}

#' Fedorov-Wynn algorithm in Rcpp.
#' Wrapper around the compiled routine in \code{src/Fedorovwynnalgorithm.cpp}
#'
#' @param protocols List of elementary protocols.
#' @param ndimen    Integer vector of dimensions.
#' @param nbprot    Integer vector of protocol counts.
#' @param numprot   Integer vector of protocol indices.
#' @param freq      Numeric vector of frequencies.
#' @param nbdata    Integer vector of data counts.
#' @param vectps    Numeric vector of sampling times.
#' @param fisher    Numeric vector of Fisher information values.
#' @param error     Integer error/status code.
#' @param protdep   Integer protocol-dependence flags.
#' @param freqdep   Numeric frequency-dependence values.
#' @return A list with the results of the Fedorov-Wynn algorithm.
#' @name FedorovWynnAlgorithm_Rcpp
#' @keywords internal
NULL

#' .algoBars - shared ggplot bar chart for algorithm output frequencies/weights.
#' @param data   data.frame with columns named by \code{x} and \code{y}.
#' @param x      Name of the x-axis column (character).
#' @param y      Name of the y-axis column (character).
#' @param xlab   x-axis label.
#' @param ylab   y-axis label.
#' @keywords internal
.algoBars = function( data, x, y, xlab = x, ylab = y ) {
  ggplot( data, aes( x = stats::reorder( .data[[ x ]], .data[[ y ]] ), y = .data[[ y ]] ) ) +
    geom_bar( stat = "identity", fill = "gray50" ) +
    scale_y_continuous( limits = c( 0, max( data[[ y ]], 1 ) ),
                        expand = c( 0, 0 ) ) +
    scale_x_discrete( expand = c( 0, 0 ) ) +
    labs( x = xlab, y = ylab ) +
    coord_flip() +
    theme_minimal( base_size = 14 ) +
    theme(
      plot.title   = element_text( hjust = 0.5, face = "bold" ),
      axis.title   = element_text( color = "black" ),
      axis.text    = element_text( color = "black" ),
      panel.grid.major.x = element_line( color = "gray90", linewidth = 0.5 ),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect( color = "gray80", fill = NA, linewidth = 0.5 )
    )
}

#' .constraintsKbl - build a kableExtra table from arm constraints.
#' @param arms                  List-of-arms
#' @param optimizationAlgorithm An Optimization subclass object.
#' @param col_names             Character vector of column names for the table.
#' @keywords internal
.constraintsKbl = function( arms, optimizationAlgorithm, col_names ) {
  armsConstraints = map( pluck( arms, 1L ), ~ getArmConstraints( .x, optimizationAlgorithm ) )
  df = map( armsConstraints, ~ map( .x, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) ) |>
    list_flatten() |>
    list_rbind()
  colnames( df ) = col_names
  kbl( df, align = c( "l", rep( "c", ncol( df ) - 1L ) ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
}

#' FedorovWynnAlgorithm: S7 class
#' @title FedorovWynnAlgorithm
#' @description
#' Fedorov-Wynn exchange algorithm for D-optimal sampling-time design (Rcpp kernel).
#' @param elementaryProtocols       List of elementary protocols.
#' @param numberOfSubjects          Total number of subjects.
#' @param proportionsOfSubjects     Initial proportion vector for each protocol.
#' @param showProcess               Logical; print FIM enumeration progress.
#' @param FedorovWynnAlgorithmOutputs Output list from the Rcpp routine.
#' @export

FedorovWynnAlgorithm = new_class( "FedorovWynnAlgorithm",
                                  package    = "PFIM",
                                  properties = list(
                                    elementaryProtocols         = new_property( class_list,    default = list()  ),
                                    numberOfSubjects            = new_property( class_vector,  default = 0.0     ),
                                    proportionsOfSubjects       = new_property( class_vector,  default = 0.0     ),
                                    showProcess                 = new_property( class_logical, default = FALSE   ),
                                    FedorovWynnAlgorithmOutputs = new_property( class_list,    default = list()  )
                                  )
)
S4_register( FedorovWynnAlgorithm )

plotFrequenciesFedorovWynnAlgorithm = new_generic(
  "plotFrequenciesFedorovWynnAlgorithm",
  c( "optimization", "optimizationAlgorithm" )
)

#' Search for a D-optimal sampling design
#' @name optimizeDesign
#' @export

method( optimizeDesign, list( Optimization, FedorovWynnAlgorithm ) ) =
  function( optimizationObject, optimizationAlgorithm ) {

    p                          = projectProp( optimizationObject, "optimizerParameters" )
    showProcess                = isTRUE( p$showProcess )
    initialSamplings           = p$elementaryProtocols
    totalNumberOfIndividuals   = p$numberOfSubjects
    proportionsOfSubjects      = p$proportionsOfSubjects
    # Flatten initial elementary protocols into a single vector for row matching
    initialElementaryProtocols = unlist( initialSamplings )
    totalCost                  = sum( lengths( initialSamplings ) * totalNumberOfIndividuals )

    design      = projectProp( optimizationObject, "designs" )[[ 1L ]]
    optimalDesign = design
    designName  = prop( design, "name" )

    fimsFromConstraints = generateFimsFromConstraints( optimizationObject )

    # Reduce list-of-row-vectors into matrices expected by the C routine
    samplingsForFedorovWynn = reduce( fimsFromConstraints$samplingsForFedorovWynnAlgo[[ designName ]], rbind )
    fisherMatrices          = reduce( fimsFromConstraints$listFimsAlgoFW[[            designName ]], rbind )
    listArms                = fimsFromConstraints$listArms[[ designName ]]

    # Assemble input structure
    nProtocols = nrow( samplingsForFedorovWynn )
    nTimes     = ncol( samplingsForFedorovWynn )
    ndimFim    = fimsFromConstraints$dimFim

    elementaryProtocolsFW = list(
      numberOfprotocols = nProtocols,
      numberOfTimes     = nTimes,
      nbOfDimensions    = ndimFim,
      totalCost         = totalCost,
      samplingTimes     = samplingsForFedorovWynn,
      fisherMatrices    = fisherMatrices
    )

    ndimen  = c( nProtocols, ndimFim, totalCost )
    npInit  = nProtocols
    buf     = .fedorovWynnBufferSizes( ndimFim, nProtocols )
    nMaxPop = buf$nMaxPop
    nFisher = buf$nFisher
    nBuf    = buf$nBuf

    numprot = rep( 0, nBuf )
    freq    = rep( 0, nBuf )
    nbdata  = rep( 0, nBuf )
    vectps  = rep( 0, nBuf )
    fisher  = rep( 0, nFisher )
    nok     = 0L

    # Identify which rows of the sampling grid match the user-specified protocols
    indexElemProt = which(
      rowSums( samplingsForFedorovWynn ==
                 initialElementaryProtocols[ col( samplingsForFedorovWynn ) ] ) == nTimes
    )
    nElemProt = length( indexElemProt )

    zeprot = c( nElemProt, indexElemProt )
    zefreq = rep( 0, nBuf )
    zefreq[ seq_len( length( proportionsOfSubjects ) ) ] = proportionsOfSubjects

    # Run C routine
    prop( optimizationAlgorithm, "showProcess" ) = showProcess

    if ( showProcess )
      message( "Fedorov-Wynn exchange algorithm started" )

    output = FedorovWynnAlgorithm_Rcpp(
      elementaryProtocolsFW, ndimen, npInit,
      numprot, freq, nbdata, vectps, fisher, nok,
      zeprot, zefreq, showProcess
    )

    # Parse outputs
    if ( !any( rowSums( output$optimal_sampling_times ) > 0 ) )
      stop( "FedorovWynn algorithm did not converge." )

    if ( showProcess )
      message( "Fedorov-Wynn exchange algorithm finished" )

    optimalFrequencies       = output$freq[   output$freq    > 0 ]
    indexOptimalArms         = output$numprot[ output$numprot > 0 ]
    listArms                 = listArms[ indexOptimalArms ]
    numberOfIndividuals      = as.double( p$numberOfSubjects * optimalFrequencies )

    prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" ) = list(
      listArms            = listArms,
      optimalFrequencies  = optimalFrequencies,
      numberOfIndividuals = numberOfIndividuals
    )

    fim         = projectProp( optimizationObject, "fim" )
    optimalArms = setOptimalArms( fim, optimizationAlgorithm )
    prop( optimalDesign, "arms" ) = map( optimalArms, ~ .x$arm )

    evaluationsToRun = list(
      .evaluationFromOptimization( optimizationObject, optimalDesign ),
      .evaluationFromOptimization( optimizationObject, design )
    )
    evaluationResults = .pfimRunEvaluations( evaluationsToRun )

    prop( optimizationObject, "optimisationDesign" ) = list(
      evaluationInitialDesign = evaluationResults[[ 2L ]],
      evaluationOptimalDesign = evaluationResults[[ 1L ]]
    )
    prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list(
      optimizationAlgorithm = optimizationAlgorithm,
      optimalArms           = optimalArms,
      frequencies           = optimalFrequencies
    )

    optimizationObject
  }

#' Frequency trajectories of the Fedorov-Wynn algorithm
#' @name plotFrequenciesFedorovWynnAlgorithm
#' @export

method( plotFrequenciesFedorovWynnAlgorithm,
        list( Optimization, FedorovWynnAlgorithm ) ) =
  function( optimization, optimizationAlgorithm ) {

    out  = prop( optimization, "optimisationAlgorithmOutputs" )
    data = data.frame(
      arm       = map_chr( out$optimalArms, ~ prop( .x$arm, "name" ) ),
      frequency = out$frequencies
    )
    .algoBars( data, x = "arm", y = "frequency", xlab = "Arm", ylab = "Frequency" )
  }

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @export

method( constraintsTableForReport, FedorovWynnAlgorithm ) =
  function( optimizationAlgorithm, arms ) {
    .constraintsKbl(
      arms, optimizationAlgorithm,
      c( "Arms name", "Number of subjects", "Outcome",
         "Initial samplings", "Fixed times",
         "Number of samplings optimisable", "Dose constraints" )
    )
  }
