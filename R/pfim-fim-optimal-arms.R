# Discrete optimizer helpers: Hamilton allocation and setOptimalArms builders.

#' Proportional subject counts per arm that sum exactly to \code{total}.
#'
#' Uses largest-remainder (Hamilton) rounding at \code{digits} decimal places so
#' \code{sum(n) == round(total, digits)} (within floating tolerance). Default
#' \code{digits = 0} yields integer subject counts (implementable designs).
#' Naive \code{round(total * w)} can miss the study size by ±1 after rounding.
#' @include Fim.R
#' @noRd
#' @keywords internal
.allocProportionalSubjects = function( total, weights, digits = 0L ) {
  w = as.numeric( weights )
  if ( !length( w ) )
    return( numeric( 0 ) )
  sw = sum( w )
  if ( !is.finite( sw ) || sw <= 0 )
    .pfimInternalStop( "protocol weights must sum to a positive finite value." )
  w = w / sw
  total = as.numeric( total )
  if ( !is.finite( total ) || total < 0 )
    .pfimInternalStop( "total subject count must be a non-negative finite value." )

  scale  = 10^as.integer( digits )
  target = as.integer( round( total * scale ) )
  exact  = total * w * scale
  base   = floor( exact + 1e-10 )
  rem    = target - as.integer( sum( base ) )
  if ( rem != 0L ) {
    frac = exact - base
    # Largest remainders first when rem > 0; smallest first when rem < 0.
    ord  = if ( rem > 0L )
      order( frac, decreasing = TRUE, method = "radix" )
    else
      order( frac, decreasing = FALSE, method = "radix" )
    nAdj = abs( rem )
    step = if ( rem > 0L ) 1L else -1L
    base[ ord[ seq_len( nAdj ) ] ] = base[ ord[ seq_len( nAdj ) ] ] + step
  }
  as.numeric( base / scale )
}

#' Total N for population multiplicative allocation.
#'
#' Prefers \code{numberOfSubjects} (aligned with Fedorov-Wynn), then legacy
#' \code{numberOfArms} when it stores the study size, else the first cell arm size.
#' @noRd
#' @keywords internal
.multiplicativePopulationNTotal = function( optimizationAlgorithm, cells, weightsIndex,
                                            out = NULL ) {
  if ( is.null( out ) )
    out = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
  nSubjects = out$numberOfSubjects
  if ( is.null( nSubjects ) )
    nSubjects = out$numberOfArms
  if ( length( nSubjects ) == 1L && is.finite( nSubjects ) && nSubjects > 0 )
    return( as.double( nSubjects ) )
  first = .multCellAt( cells, weightsIndex, 1L )
  prop( .constraintCellArms( first )[[ 1L ]], "size" )
}


#' Resolve retained Mult cells from outputs (\code{listArms} preferred).
#' @noRd
#' @keywords internal
.multListArmsFromOut = function( out ) {
  cells = out$listArms
  if ( is.null( cells ) )
    cells = out$armFims
  if ( is.null( cells ) || !length( cells ) )
    stop(
      "MultiplicativeAlgorithm: no listArms/armFims available for setOptimalArms.",
      call. = FALSE
    )
  cells
}

#' Map Mult weight slots to cells (retained list or full grid).
#' @noRd
#' @keywords internal
.multCellAt = function( cells, weightsIndex, j ) {
  # Retained listArms: length matches weights; full grid: index by weightsIndex.
  if ( length( cells ) == length( weightsIndex ) )
    cells[[ j ]]
  else
    cells[[ weightsIndex[[ j ]] ]]
}

#' @noRd
#' @keywords internal
.setOptimalArmsMultiplicative = function( optimizationAlgorithm, out = NULL ) {
  # Same Hamilton allocation for all FIM types (weights -> implementable sizes).
  .setOptimalArmsMultiplicativePopulation( optimizationAlgorithm, out = out )
}

#' Build population arms from Multiplicative weights.
#'
#' Keeps candidates with weight \code{> weightThreshold}, then allocates subjects
#' with \code{.allocProportionalSubjects()} (exact sum to N). Joint (multi-outcome)
#' cells keep only the single best-weighted protocol.
#' @noRd
#' @keywords internal
.setOptimalArmsMultiplicativePopulation = function( optimizationAlgorithm, out = NULL ) {
  if ( is.null( out ) )
    out = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
  weightsIndex = out$weightsIndex
  weights      = out$optimalWeights
  cells        = .multListArmsFromOut( out )
  if ( !length( weightsIndex ) )
    weightsIndex = seq_along( weights )
  if ( length( weights ) != length( weightsIndex ) )
    stop( "MultiplicativeAlgorithm: optimalWeights/weightsIndex length mismatch.",
          call. = FALSE )

  # Joint multi-outcome cell: one protocol only (max weight among retained).
  # Split study N across that protocol's arms (same exact-sum rule as FW).
  firstCell = .multCellAt( cells, weightsIndex, 1L )
  if ( length( .constraintCellArms( firstCell ) ) > 1L ) {
    best  = which.max( weights )
    arms  = .cloneConstraintCellArms( .multCellAt( cells, weightsIndex, best ) )
    N_total = .multiplicativePopulationNTotal(
      optimizationAlgorithm, cells, weightsIndex, out = out
    )
    fracs = map_dbl( arms, ~ prop( .x, "size" ) )
    sizes = .allocProportionalSubjects( N_total, fracs )
    arms  = map( seq_along( arms ), function( j ) {
      prop( arms[[ j ]], "size" ) = sizes[[ j ]]
      arms[[ j ]]
    } )
    return( keep( arms, ~ prop( .x, "size" ) > 0 ) )
  }

  # Convert continuous weights to group sizes that sum exactly to N_total.
  N_total   = .multiplicativePopulationNTotal(
    optimizationAlgorithm, cells, weightsIndex, out = out
  )
  nPerGroup = .allocProportionalSubjects( N_total, weights )
  armList   = map2( seq_along( weightsIndex ), nPerGroup, function( j, size ) {
    arm = .cloneConstraintArm( .multCellAt( cells, weightsIndex, j ) )
    prop( arm, "size" ) = size
    prop( arm, "name" ) = paste0( "Arm", weightsIndex[[ j ]] )
    arm
  })
  # Drop empty arms (Hamilton can assign 0 when many tiny weights).
  armList = keep( armList, ~ prop( .x, "size" ) > 0 )
  armList[ order( map_dbl( armList, ~ prop( .x, "size" ) ), decreasing = TRUE ) ]
}

#' @noRd
#' @keywords internal
.setOptimalArmsFedorovWynn = function( optimizationAlgorithm, ... ) {
  # Same Hamilton allocation for all FIM types.
  .setOptimalArmsFedorovWynnPopulation( optimizationAlgorithm, ... )
}

#' Build population arms from Fedorov-Wynn frequencies.
#'
#' Each optimal protocol gets \code{N} from \code{numberOfSubjects}
#' (alias \code{numberOfIndividuals}). Counts sum to the study size via
#' \code{.allocProportionalSubjects()}. When a cell holds
#' several arm entries, sizes are split by relative entry sizes with the same
#' exact-sum allocator.
#' @noRd
#' @keywords internal
.setOptimalArmsFedorovWynnPopulation = function( optimizationAlgorithm, ... ) {
  outputs             = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
  numberOfIndividuals = outputs$numberOfSubjects
  if ( is.null( numberOfIndividuals ) )
    numberOfIndividuals = outputs$numberOfIndividuals
  arms = unlist( imap( outputs$listArms, function( cell, i ) {
    entries = if ( !is.null( cell$entries ) ) cell$entries else list( cell )
    arms    = map( entries, ~ .pfimCloneS7( .pfimAsArm( .x ) ) )
    N       = as.double( numberOfIndividuals[[ i ]] )
    if ( !is.finite( N ) || N <= 0 )
      return( list() )
    if ( length( arms ) == 1L ) {
      prop( arms[[ 1L ]], "size" ) = N
      prop( arms[[ 1L ]], "name" ) = paste0( "Arm", i )
      return( list( arms[[ 1L ]] ) )
    }
    # Split N across joint entries by relative original sizes (exact sum).
    fracs = map_dbl( arms, ~ prop( .x, "size" ) )
    sizes = .allocProportionalSubjects( N, fracs )
    map( seq_along( arms ), function( j ) {
      prop( arms[[ j ]], "size" ) = sizes[[ j ]]
      arms[[ j ]]
    })
  } ), recursive = FALSE )
  # Drop empty supports (equivalent candidates / tiny frequencies).
  keep( arms, ~ prop( .x, "size" ) > 0 )
}

