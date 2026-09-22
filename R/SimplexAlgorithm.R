#' @title SimplexAlgorithm
#' @description
#' Nelder-Mead (amoeba) optimization for sampling times (Rcpp kernel).
#' Pass \code{pctInitialSimplexBuilding}, \code{tolerance}, and \code{maxIteration}
#' via \code{optimizerParameters} on \code{\link{Optimization}}.
#' Search starts from each outcome's \code{initialSamplings} (same flat layout as PSO/PGBO).
#' Extra vertices move one coordinate toward the far end of its sampling window
#' (not toward 0), so unequal \code{numberOfTimesByWindows} stays feasible.
#' Infeasible \code{initialSamplings} are replaced by a feasible window start
#' with a warning.
#' \code{tolerance <= 0} disables the relative stop (\code{converged = NA}, no warning),
#' matching PSO/PGBO when stall stop is off.
#' @param optimizerOutputs List filled by \code{optimizeDesign()} (optimal arms, etc.).
#' @return A \code{SimplexAlgorithm} specification object.
#' @examples
#' \dontrun{
#' vignette("Example02")
#' }
#' @include Optimization.R
#' @export

SimplexAlgorithm = new_class( "SimplexAlgorithm", package = "PFIM",
                              properties = list(
                                optimizerOutputs = new_property( class_list, default = list() )
                              ) )
S4_register( SimplexAlgorithm )

#' Nelder-Mead simplex minimiser (Rcpp).
#'
#' Compiled implementation in \code{src/SimplexAlgorithm.cpp}.
#' The objective \code{funk} remains an R callback because FIM evaluation stays in R.
#' Vertices are flat sampling vectors (same layout as PSO); C++ never sees arms.
#'
#' @param p Numeric matrix of simplex vertices (rows).
#' @param y Numeric vector of objective values at vertices.
#' @param ftol Relative tolerance on the spread of \code{y}.
#' @param itmax Maximum number of iterations.
#' @param funk R function \code{funk(data, pr, outcomes)}.
#' @param outcomes Outcome structure passed to \code{funk} (unused on flat path).
#' @param data Optimization object passed to \code{funk}.
#' @param show_process Logical; print progress.
#' @return List with \code{p}, \code{y}, \code{iterations}, \code{converged}, \code{results}.
#' @name fun_amoeba_Rcpp
#' @keywords internal
NULL

#' Flat sampling vector -> \code{1/D} for Nelder-Mead (invalid -> large penalty).
#'
#' Used by \code{optimizeDesign(SimplexAlgorithm)}. Same flat layout as PSO;
#' fitness is cost \code{1/D} (PGBO evaluates raw D in its kernel).
#' @param simplex Numeric flat sampling vector (one vertex).
#' @param layout Flat sampling layout from \code{.buildFlatSamplingLayout()}.
#' @param evalCtx Reusable context from \code{.pfimMetaheuristicEvalContext()}.
#' @return Numeric cost (length 1).
#' @noRd
#' @keywords internal
.simplexFlatFitness = function( simplex, layout, evalCtx ) {
  flat = unname( as.numeric( simplex ) )
  if ( length( flat ) != length( layout$initial_flat ) )
    .pfimStop(
      paste0(
        "SimplexAlgorithm: vertex length (", length( flat ),
        ") does not match flat sampling layout (", length( layout$initial_flat ), ")."
      )
    )
  flat = .pfimSortFlatByGroups( flat, layout$sorting_groups )
  if ( !.isFlatValid( layout, flat ) )
    return( .metaheuristicFitnessPenalty )
  d = 1 / Dcriterion( evalCtx( flat ) )
  if ( !is.finite( d ) ) .metaheuristicFitnessPenalty else d
}

#' Nelder-Mead continuous D-optimal design (delegates to multi-design driver).
#'
#' @param optimizationObject An \code{\link{Optimization}} project.
#' @param optimizationAlgorithm A \code{SimplexAlgorithm} instance.
#' @return Updated \code{Optimization} after optimizing each design in turn.
#' @name optimizeDesign
#' @keywords internal

method( optimizeDesign, list( Optimization, SimplexAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm )
{
  .pfimContinuousOptimizeDesigns(
    optimizationObject,
    optimizationAlgorithm,
    .optimizeSimplexOneDesign
  )
}

#' @noRd
#' @keywords internal
.optimizeSimplexOneDesign = function( optimizationObject, optimizationAlgorithm )
{
  optimizerParameters = projectProp( optimizationObject, "optimizerParameters" )
  showProcess = optimizerParameters$showProcess
  pctInitialSimplexBuilding = optimizerParameters$pctInitialSimplexBuilding
  tolerance = optimizerParameters$tolerance
  maxIteration = optimizerParameters$maxIteration

  prep = .pfimPrepareContinuousDesign( optimizationObject )
  design        = prep$design
  arms          = prep$arms
  initialDesign = .pfimCloneS7( design )
  layout        = .buildFlatSamplingLayout( design )
  evalTemplate  = .evaluationFromOptimization( optimizationObject, initialDesign, name = "" )
  evalCtx       = .pfimMetaheuristicEvalContext( evalTemplate, design, arms )

  samplingsSimplex = .pfimSimplexStartMatrix( layout, pctInitialSimplexBuilding )
  .pfimRequireFeasibleSimplex( layout, samplingsSimplex, "initial simplex" )

  # Flat fitness with layout / eval context closed over once (cost = 1/D).
  # Sort-before-eval: spacing constraints assume ordered times within each group.
  funk = function( optimizationObject, simplex, outcomes ) {
    .simplexFlatFitness( simplex, layout, evalCtx )
  }

  y = map_dbl(
    seq_len( nrow( samplingsSimplex ) ),
    ~ funk( optimizationObject, samplingsSimplex[ .x, ], NULL )
  )

  opti = fun_amoeba_Rcpp(
    samplingsSimplex, y, tolerance, maxIteration,
    funk, optimizationObject, NULL, showProcess
  )

  feas = .pfimRequireFeasibleSimplex( layout, opti$p, "search result" )
  y_feas = as.numeric( opti$y )
  y_feas[ !feas ] = Inf

  .pfimWarnIfNotConverged(
    "SimplexAlgorithm", tolerance, opti$converged,
    extra = paste0(
      " within maxIteration = ", maxIteration,
      " (tolerance = ", tolerance, ")"
    )
  )

  best_row  = which.min( y_feas )[ 1L ]
  flat_best = .pfimSortFlatByGroups(
    unname( as.numeric( opti$p[ best_row, ] ) ), layout$sorting_groups
  )
  finalArms     = .applyFlatToArms( flat_best, arms )
  optimalDesign = .pfimOptimalDesignFrom( initialDesign, finalArms )

  algoOut = .pfimContinuousAlgoOutputs(
    list(
      converged  = opti$converged,
      iterations = as.integer( opti$iterations ),
      bestCost   = min( y_feas )
    ),
    tolerance = tolerance
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

method( constraintsTableForReport, SimplexAlgorithm ) = function( optimizationAlgorithm, arms ) {
  .pfimConstraintsTableContinuous( optimizationAlgorithm, arms )
}

#' Sort flat sampling coordinates within each outcome group.
#' @noRd
#' @keywords internal
.pfimSortFlatByGroups = function( flat, sorting_groups ) {
  reduce( sorting_groups, function( acc, idx ) {
    acc[ idx ] = sort( acc[ idx ] )
    acc
  }, .init = flat )
}

#' Sequential window bounds for each flat coordinate (after group-wise sort).
#' First \code{numberOfTimesByWindows[1]} times of a group belong to window 1.
#' @noRd
#' @keywords internal
.pfimFlatCoordinateWindows = function( layout ) {
  n  = length( layout$initial_flat )
  lo = rep( NA_real_, n )
  hi = rep( NA_real_, n )
  for ( g in seq_along( layout$group_specs ) ) {
    spec = layout$group_specs[[ g ]]
    idx  = layout$sorting_groups[[ g ]]
    tab  = .pfimSamplingWindowsTable( spec$constraint )
    n_by = as.integer( prop( spec$constraint, "numberOfTimesByWindows" ) )
    if ( length( n_by ) == 1L && nrow( tab ) > 1L )
      n_by = rep( n_by, nrow( tab ) )
    pos = 0L
    for ( w in seq_len( nrow( tab ) ) ) {
      k = n_by[[ w ]]
      if ( k <= 0L ) next
      take = idx[ pos + seq_len( k ) ]
      lo[ take ] = tab$min[[ w ]]
      hi[ take ] = tab$max[[ w ]]
      pos = pos + k
    }
  }
  list( lo = lo, hi = hi )
}

#' Feasible Nelder-Mead start: user times if they satisfy the windows, else
#' equally spaced points inside each window.
#' @noRd
#' @keywords internal
.pfimFeasibleSimplexStart = function( layout ) {
  start = .pfimSortFlatByGroups( layout$initial_flat, layout$sorting_groups )
  if ( .isFlatValid( layout, start ) )
    return( list( flat = start, reseeded = FALSE ) )
  list(
    flat = .sampleFlatFromConstraints( layout, .evenSamplingsFromSamplingConstraints ),
    reseeded = TRUE
  )
}

#' Nelder-Mead start matrix inside sampling windows.
#'
#' Vertex 1 is a feasible start. Vertex \eqn{i+1} moves coordinate \eqn{i}
#' toward the far end of its sequential window (not toward 0). Vertices that
#' still fail counts or \code{minSampling} are replaced by a feasible sample.
#' @noRd
#' @keywords internal
.pfimSimplexStartMatrix = function( layout, pct ) {
  seeded = .pfimFeasibleSimplexStart( layout )
  if ( seeded$reseeded )
    .pfimWarn(
      "SimplexAlgorithm: initialSamplings do not satisfy the sampling windows ",
      "(counts or minSampling). Starting from a feasible point inside the windows."
    )
  start = seeded$flat
  n = length( start )
  mat = matrix( start, nrow = n + 1L, ncol = n, byrow = TRUE )
  bounds = .pfimFlatCoordinateWindows( layout )
  frac = as.numeric( pct ) / 100
  if ( !is.finite( frac ) || frac < 0 ) frac = 0
  for ( i in seq_len( n ) ) {
    lo = bounds$lo[[ i ]]
    hi = bounds$hi[[ i ]]
    cur = start[[ i ]]
    if ( !is.finite( lo ) || !is.finite( hi ) ) next
    toward = if ( abs( cur - lo ) <= abs( cur - hi ) ) hi else lo
    moved = cur + frac * ( toward - cur )
    mat[ i + 1L, i ] = min( hi, max( lo, moved ) )
  }
  mat = t( apply(
    mat, 1L,
    function( row ) .pfimSortFlatByGroups( row, layout$sorting_groups )
  ) )
  for ( r in seq_len( nrow( mat ) ) ) {
    if ( .isFlatValid( layout, mat[ r, ] ) ) next
    mat[ r, ] = .sampleFlatFromConstraints( layout )
  }
  mat
}

#' Which simplex rows satisfy window counts and \code{minSampling}.
#' @noRd
#' @keywords internal
.pfimSimplexFeasibleRows = function( layout, mat ) {
  mat = as.matrix( mat )
  if ( !nrow( mat ) )
    return( logical( 0L ) )
  vapply(
    seq_len( nrow( mat ) ),
    function( i ) {
      flat = .pfimSortFlatByGroups(
        unname( as.numeric( mat[ i, ] ) ), layout$sorting_groups
      )
      .isFlatValid( layout, flat )
    },
    logical( 1L )
  )
}

#' Abort rather than report an all-infeasible amoeba as a converged optimum.
#' @return Logical feasibility mask when at least one row is valid.
#' @noRd
#' @keywords internal
.pfimRequireFeasibleSimplex = function( layout, mat, when ) {
  feas = .pfimSimplexFeasibleRows( layout, mat )
  if ( any( feas ) )
    return( feas )
  .pfimStop(
    "SimplexAlgorithm: ", when,
    " has no feasible sampling times. Check initialSamplings against ",
    "samplingsWindows, numberOfTimesByWindows, and minSampling."
  )
}
