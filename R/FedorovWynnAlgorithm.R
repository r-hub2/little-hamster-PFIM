#' Buffer lengths for the Fedorov-Wynn C++ interface.
#'
#' R<->C++ contract (inout vectors filled by \code{FedorovWynnAlgorithm_Rcpp}):
#' \itemize{
#'   \item \code{nFisher} = packed triangle size \eqn{p(p+1)/2} (matches C++
#'     \code{PackedFim}; Mult uses dense \eqn{p \times p} instead).
#'   \item \code{nMaxPop} = Carathéodory bound \eqn{p(p+1)/2 + 1} on support size.
#'   \item \code{nBuf} sizes \code{numprot}/\code{freq}/\code{nbdata}/\code{zefreq};
#'     \code{2 * nProtocols} leaves headroom while the exchange temporarily grows
#'     support before pruning.
#'   \item \code{nVectps} must cover one protocol's sampling-time row
#'     (\code{nTimes}), which can exceed \code{nBuf} on tiny grids.
#' }
#' @param ndimFim FIM dimension \eqn{p}.
#' @param nProtocols Number of candidate protocols on the discrete grid.
#' @param nTimes Number of sampling columns per protocol row.
#' @return List with \code{nMaxPop}, \code{nFisher}, \code{nBuf}, \code{nVectps}.
#' @noRd
#' @keywords internal
.fedorovWynnBufferSizes = function( ndimFim, nProtocols, nTimes = 0L ) {
  ndimFim    = as.integer( ndimFim )
  nProtocols = as.integer( nProtocols )
  nTimes    = as.integer( nTimes )
  nMaxPop    = ndimFim * ( ndimFim + 1L ) / 2L + 1L
  nFisher    = ndimFim * ( ndimFim + 1L ) / 2L
  nBuf       = max( nProtocols * 2L, nMaxPop )
  nVectps    = max( nBuf, nTimes )
  list( nMaxPop = nMaxPop, nFisher = nFisher, nBuf = nBuf, nVectps = nVectps )
}

#' Map Fedorov-Wynn C++ status to stop / warn / continue.
#' @return \code{TRUE} if stationary (\code{success}), else \code{FALSE} (\code{incomplete}).
#' @noRd
#' @keywords internal
.pfimFedorovWynnCheckStatus = function( status ) {
  if ( identical( status, "success" ) )
    return( TRUE )
  if ( identical( status, "incomplete" ) ) {
    warning(
      "Fedorov-Wynn: outer cycle budget reached without a certified optimum ",
      "(status = incomplete; max Tr(F^{-1} F_i) <= p(1+delta)). ",
      "Frequencies are usable but not certified optimal.",
      call. = FALSE
    )
    return( FALSE )
  }
  detail = switch(
    as.character( status ),
    singular_fim = paste0(
      "mixture FIM became singular. Check elementaryProtocols, ",
      "candidate grid, and proportionsOfSubjects."
    ),
    paste0( "status = ", status, "." )
  )
  .pfimStop( "Fedorov-Wynn: ", detail )
}

#' Flatten one FW protocol entry to a numeric grid-row vector.
#' Accepts a numeric vector or a list of per-outcome numeric vectors
#' (each outcome block is sorted to match the sampling grid).
#' @noRd
#' @keywords internal
.pfimFlattenFedorovProtocol = function( protocol ) {
  if ( is.list( protocol ) && !is.numeric( protocol ) ) {
    parts = lapply( protocol, function( block ) {
      sort( as.numeric( unlist( block, use.names = FALSE ) ) )
    } )
    return( as.numeric( unlist( parts, use.names = FALSE ) ) )
  }
  as.numeric( unlist( protocol, use.names = FALSE ) )
}

#' Normalize \code{elementaryProtocols} for Fedorov-Wynn grid matching.
#'
#' Supported forms per support point: flat vector of length \code{nTimes}, or a
#' list of outcome vectors that concatenate to \code{nTimes}. If the user passes
#' \code{list(pk, pd, ...)} with a single proportion and the concatenated length
#' equals \code{nTimes}, treat that list as one multi-outcome protocol (vignette
#' PK/PD convenience).
#' @noRd
#' @keywords internal
.pfimNormalizeFedorovElementaryProtocols = function( protocols, proportions, nTimes ) {
  if ( !is.list( protocols ) || !length( protocols ) )
    stop( "Fedorov-Wynn: elementaryProtocols must be a non-empty list.", call. = FALSE )
  if ( !is.numeric( proportions ) || !length( proportions ) )
    stop( "Fedorov-Wynn: proportionsOfSubjects must be a non-empty numeric vector.", call. = FALSE )

  nProt = length( protocols )
  nProp = length( proportions )

  # Vignette-style PK/PD: list(pk_vec, pd_vec) + one proportion -> one protocol.
  if ( nProp == 1L && nProt != 1L ) {
    flat_lens = vapply( protocols, function( p ) length( .pfimFlattenFedorovProtocol( p ) ), integer( 1L ) )
    if ( sum( flat_lens ) == as.integer( nTimes ) && all( flat_lens > 0L ) ) {
      protocols = list( protocols )
      nProt = 1L
    }
  }

  if ( nProp != nProt )
    stop(
      "Fedorov-Wynn: proportionsOfSubjects length (", nProp,
      ") must match elementaryProtocols length (", nProt,
      "). For PK/PD, pass one protocol as list(list(pk, pd)) or a single flat ",
      "vector of length ", nTimes, " (outcomes concatenated).",
      call. = FALSE
    )

  if ( any( !is.finite( proportions ) ) || any( proportions < 0 ) )
    stop( "Fedorov-Wynn: proportionsOfSubjects must be finite and non-negative.", call. = FALSE )
  .validateUnitProportions( proportions, "proportionsOfSubjects" )

  list( protocols = protocols, proportions = as.numeric( proportions ) )
}

#' Match one normalized protocol vector to a unique sampling-grid row.
#' @noRd
#' @keywords internal
.pfimMatchFedorovProtocolRow = function( protocol, samplingsForFedorovWynn, index ) {
  nProtocols = nrow( samplingsForFedorovWynn )
  nTimes     = ncol( samplingsForFedorovWynn )
  flat = .pfimFlattenFedorovProtocol( protocol )
  if ( length( flat ) != nTimes )
    stop(
      "Fedorov-Wynn: elementary protocol ", index, " has length ", length( flat ),
      " but the sampling grid has ", nTimes, " times ",
      "(concatenate all outcomes in arm order, e.g. PK then PD).",
      call. = FALSE
    )
  # Near-equality: grid times are numeric; exact == is brittle for user input.
  scale = pmax( 1, abs( flat ) )
  tol   = sqrt( .Machine$double.eps ) * scale
  hits = which( vapply( seq_len( nProtocols ), function( r ) {
    all( abs( samplingsForFedorovWynn[ r, ] - flat ) <= tol )
  }, logical( 1L ) ) )
  if ( !length( hits ) )
    stop(
      "Fedorov-Wynn: elementary protocol ", index,
      " was not found in the sampling grid: (",
      paste( flat, collapse = ", " ), ").",
      call. = FALSE
    )
  if ( length( hits ) > 1L )
    warning(
      "Fedorov-Wynn: elementary protocol ", index,
      " matches ", length( hits ), " grid rows; using the first.",
      call. = FALSE
    )
  as.integer( hits[[ 1L ]] )
}

#' Fedorov-Wynn algorithm in Rcpp.
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
#' @param show_process Logical; print cycle progress.
#' @param delta     Relative optimality gap: stop when
#'   \eqn{\max_i \mathrm{Tr}(F^{-1} F_i) \le p(1+\delta)} (default \code{1e-4}).
#' @return A list with the results of the Fedorov-Wynn algorithm.
#' @name FedorovWynnAlgorithm_Rcpp
#' @keywords internal
NULL

#' Bar chart for algorithm weights / frequencies.
#' @noRd
#' @keywords internal
.algoBars = function( data, x, y, xlab = x, ylab = y, ylim = NULL, breaks = NULL ) {
  y_max = if ( is.null( ylim ) ) max( data[[ y ]], 1 ) else ylim[ 2L ]
  y_min = if ( is.null( ylim ) ) 0 else ylim[ 1L ]
  # NULL breaks = ggplot2 defaults; do not pass breaks = NULL explicitly.
  scale_y = if ( is.null( breaks ) ) {
    scale_y_continuous( limits = c( y_min, y_max ), expand = expansion( mult = c( 0, 0.06 ) ) )
  } else {
    scale_y_continuous(
      limits = c( y_min, y_max ), breaks = breaks,
      expand = expansion( mult = c( 0, 0.06 ) )
    )
  }
  ggplot( data, aes( x = stats::reorder( .data[[ x ]], .data[[ y ]] ), y = .data[[ y ]] ) ) +
    geom_col( fill = "gray50" ) +
    scale_y +
    scale_x_discrete( expand = c( 0, 0 ) ) +
    labs( x = xlab, y = ylab ) +
    coord_flip( clip = "off" ) +
    .pfimBaseTheme() +
    theme(
      # Flipped bars: keep category labels horizontal (override base theme angle = 90).
      axis.text.x = element_text( angle = 0, vjust = 0.5 ),
      axis.text.y = element_text( angle = 0, hjust = 1 ),
      panel.grid.major.x = element_line( color = "gray90", linewidth = 0.5 ),
      panel.grid.major.y = element_blank()
    )
}

#' FedorovWynnAlgorithm: S7 class
#' @title FedorovWynnAlgorithm
#' @description
#' Fedorov-Wynn exchange algorithm for D-optimal sampling-time design (Rcpp kernel).
#' Pass \code{elementaryProtocols}, \code{numberOfSubjects}, and
#' \code{proportionsOfSubjects} via \code{optimizerParameters} on \code{\link{Optimization}}.
#' Optional \code{delta} (default \code{1e-4}) is the relative optimality
#' gap \eqn{\max_i \phi_i \le p(1+\delta)}, same role as in the multiplicative algorithm.
#' Each elementary protocol may be a flat numeric vector (full grid row) or a list of
#' per-outcome vectors that concatenate to that row (typical PK/PD). A vignette-style
#' \code{list(pk, pd)} with a single proportion is accepted as one multi-outcome protocol.
#' @param FedorovWynnAlgorithmOutputs Output list from the Rcpp routine.
#' @return A \code{FedorovWynnAlgorithm} specification object for \code{run(Optimization)}.
#' @examples
#' \dontrun{
#' vignette("Example01")
#' }
#' @include Optimization.R
#' @export

FedorovWynnAlgorithm = new_class( "FedorovWynnAlgorithm",
                                  package    = "PFIM",
                                  properties = list(
                                    FedorovWynnAlgorithmOutputs = new_property( class_list, default = list() )
                                  )
)
S4_register( FedorovWynnAlgorithm )

plotFrequenciesFedorovWynnAlgorithm = new_generic(
  "plotFrequenciesFedorovWynnAlgorithm",
  c( "optimization", "optimizationAlgorithm" )
)

#' Fedorov-Wynn discrete D-optimal design (multi-design driver).
#'
#' Pipeline: enumerate constraint-grid FIMs -> match each
#' \code{elementaryProtocols} seed to a full grid row -> call C++ exchange ->
#' map frequencies to arms -> re-evaluate initial vs optimal designs.
#'
#' @param optimizationObject An \code{\link{Optimization}} project.
#' @param optimizationAlgorithm A \code{FedorovWynnAlgorithm} instance.
#' @return The same \code{Optimization} with \code{optimisationDesign} and
#'   \code{optimisationAlgorithmOutputs} filled.
#' @name optimizeDesign
#' @keywords internal

method( optimizeDesign, list( Optimization, FedorovWynnAlgorithm ) ) =
  function( optimizationObject, optimizationAlgorithm ) {
    .pfimDiscreteOptimizeDesigns(
      optimizationObject, optimizationAlgorithm, .optimizeFedorovWynnOneDesign
    )
  }

#' @noRd
#' @keywords internal
.optimizeFedorovWynnOneDesign = function( optimizationObject, optimizationAlgorithm ) {

    p                          = projectProp( optimizationObject, "optimizerParameters" )
    showProcess                = isTRUE( p$showProcess )
    # Relative optimality gap (same role as MultiplicativeAlgorithm delta).
    delta                      = if ( !is.null( p$delta ) ) as.numeric( p$delta ) else 1e-4
    if ( length( delta ) != 1L || !is.finite( delta ) || delta < 0 )
      stop( "Fedorov-Wynn: optimizerParameters$delta must be a finite >= 0 scalar.",
            call. = FALSE )
    initialSamplings           = p$elementaryProtocols
    totalNumberOfIndividuals   = p$numberOfSubjects
    proportionsOfSubjects      = p$proportionsOfSubjects
    totalCost                  = as.numeric( totalNumberOfIndividuals )

    design          = projectProp( optimizationObject, "designs" )[[ 1L ]]
    initialDesign   = .pfimCloneS7( design )
    optimalDesign   = .pfimCloneS7( design )
    designName      = prop( design, "name" )

    # One FIM (+ sampling row) per dose x sampling combination under constraints.
    fimsFromConstraints = generateFimsFromConstraints( optimizationObject )
    listArms            = fimsFromConstraints$listArms[[ designName ]]
    fim                 = projectProp( optimizationObject, "fim" )

    # Individual / Bayesian: optimum is a single protocol (covariance-mixture vertex).
    if ( .pfimIsSubjectLevelFim( fim ) ) {
      denseFims = fimsFromConstraints$listFimsAlgoMult[[ designName ]]
      best      = .pfimBestSubjectProtocol( denseFims )
      idx       = best$index
      prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" ) = list(
        listArms         = listArms[ idx ],
        optimalWeights   = 1,
        numberOfSubjects = as.double( totalNumberOfIndividuals ),
        algorithmOutput  = list(
          status    = 0L,
          converged = TRUE,
          method    = "bestSubjectProtocol"
        )
      )
      optimalArms = setOptimalArms( fim, optimizationAlgorithm )
      prop( optimalDesign, "arms" ) = optimalArms
      return( .pfimStoreOptimization(
        optimizationObject, initialDesign, optimalDesign,
        algorithmOutputs = list(
          optimizationAlgorithm = optimizationAlgorithm,
          optimalArms           = optimalArms,
          optimalWeights        = 1,
          frequencies           = 1,
          mixtureDcriterion     = best$Dcriterion,
          realisedDcriterion    = best$Dcriterion
        )
      ) )
    }

    # Stack list-of-rows into matrices expected by FedorovWynnAlgorithm_Rcpp.
    # fisherMatrices: one packed-triangle row per candidate (not dense pxp -
    # those go to Multiplicative via listFimsAlgoMult).
    samplingsForFedorovWynn = reduce( fimsFromConstraints$samplingsForFedorovWynnAlgo[[ designName ]], rbind )
    fisherMatrices          = reduce( fimsFromConstraints$listFimsAlgoFW[[            designName ]], rbind )

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
    buf     = .fedorovWynnBufferSizes( ndimFim, nProtocols, nTimes )
    nMaxPop = buf$nMaxPop
    nFisher = buf$nFisher
    nBuf    = buf$nBuf

    # Pre-allocate C++ output buffers (filled in place by the kernel).
    numprot = rep( 0, nBuf )
    freq    = rep( 0, nBuf )
    nbdata  = rep( 0, nBuf )
    vectps  = rep( 0, buf$nVectps )
    fisher  = rep( 0, nFisher )
    nok     = 0L

    # Match each user-supplied elementary protocol to one full grid row.
    # Accept flat vectors, nested per-outcome lists, and vignette-style
    # list(pk, pd) with a single proportion (normalized above).
    # Order of indexElemProt must match proportionsOfSubjects for the C++ init.
    normalized = .pfimNormalizeFedorovElementaryProtocols(
      initialSamplings, proportionsOfSubjects, nTimes
    )
    initialSamplings      = normalized$protocols
    proportionsOfSubjects = normalized$proportions

    indexElemProt = vapply(
      seq_along( initialSamplings ),
      function( i ) .pfimMatchFedorovProtocolRow(
        initialSamplings[[ i ]], samplingsForFedorovWynn, i
      ),
      integer( 1L )
    )
    nElemProt = length( indexElemProt )

    # C++ (0-based buffers): zeprot[0] = support size, zeprot[1:] = grid
    # indices in R's 1-based numbering; zefreq[0:(n-1)] = initial weights.
    zeprot = c( nElemProt, indexElemProt )
    zefreq = rep( 0, nBuf )
    zefreq[ seq_len( length( proportionsOfSubjects ) ) ] = proportionsOfSubjects

    if ( showProcess )
      message( "Fedorov-Wynn exchange algorithm started" )

    output = FedorovWynnAlgorithm_Rcpp(
      elementaryProtocolsFW, ndimen, npInit,
      numprot, freq, nbdata, vectps, fisher, nok,
      zeprot, zefreq, showProcess, delta
    )

    fw_converged = .pfimFedorovWynnCheckStatus( output$status )

    if ( showProcess )
      message( "Fedorov-Wynn exchange algorithm finished" )

    # C++ writes 1-based grid indices into numprot[1:active]; unused slots stay 0.
    # freq is renormalised so sum(active) == 1 (see kernel write-back).
    # Common mask keeps freq and numprot aligned (separate filters desync when
    # equivalent candidates leave tiny / zero weights in mismatched slots).
    nSlot  = min( length( output$freq ), length( output$numprot ) )
    active = seq_len( nSlot )[
      output$freq[ seq_len( nSlot ) ] > 0 &
        output$numprot[ seq_len( nSlot ) ] > 0
    ]
    activeFreq       = output$freq[ active ]
    indexOptimalArms = output$numprot[ active ]
    freq_sum         = sum( activeFreq )

    if ( length( activeFreq ) == 0L || abs( freq_sum - 1 ) > 1e-4 )
      stop(
        "Fedorov-Wynn algorithm did not converge (sum(freq) = ", round( freq_sum, 6 ), ").",
        call. = FALSE
      )

    optimalFrequencies = activeFreq
    listArms           = listArms[ indexOptimalArms ]
    # Subject counts that sum exactly to numberOfSubjects (Hamilton rounding).
    numberOfSubjects = .allocProportionalSubjects(
      p$numberOfSubjects, optimalFrequencies
    )

    prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" ) = list(
      listArms           = listArms,
      optimalWeights     = optimalFrequencies,
      numberOfSubjects   = numberOfSubjects,
      algorithmOutput    = list(
        status    = output$status,
        converged = fw_converged
      )
    )

    # Population FIM allocates sizes from numberOfSubjects; individual keeps size 1.
    fim         = projectProp( optimizationObject, "fim" )
    optimalArms = setOptimalArms( fim, optimizationAlgorithm )
    prop( optimalDesign, "arms" ) = optimalArms

    .pfimStoreOptimization(
      optimizationObject, initialDesign, optimalDesign,
      algorithmOutputs = list(
        optimizationAlgorithm = optimizationAlgorithm,
        optimalArms           = optimalArms,
        optimalWeights        = optimalFrequencies,
        frequencies           = optimalFrequencies
      )
    )
  }

#' Frequency trajectories of the Fedorov-Wynn algorithm
#' @name plotFrequenciesFedorovWynnAlgorithm
#' @return A \code{ggplot2} plot object.
#' @examples
#' \dontrun{
#' # help(plotFrequenciesFedorovWynnAlgorithm); vignette("Example01")
#' }
#' @keywords internal

method( plotFrequenciesFedorovWynnAlgorithm,
        list( Optimization, FedorovWynnAlgorithm ) ) =
  function( optimization, optimizationAlgorithm ) {

    plotDat = .pfimDiscreteMixturePlotData( optimization )
    .algoBars(
      data.frame(
        arm       = plotDat$data$label,
        frequency = plotDat$data$value
      ),
      x = "arm", y = "frequency",
      xlab = plotDat$xlab, ylab = "Frequency"
    )
  }

#' Constraint tables for optimization reports
#' @name constraintsTableForReport
#' @keywords internal

method( constraintsTableForReport, FedorovWynnAlgorithm ) =
  function( optimizationAlgorithm, arms ) {
    .pfimConstraintsTableDiscrete( optimizationAlgorithm, arms )
  }
