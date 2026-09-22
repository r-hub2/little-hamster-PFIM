#' @title SamplingTimeConstraints
#' @description
#' Constraints on sampling times for design optimization (windows, fixed times,
#' number of points to optimize, minimum spacing).
#'
#' Metaheuristics (Simplex / PSO / PGBO) draw candidates with
#' \code{generateSamplingsFromSamplingConstraints()} and reject infeasible ones
#' via \code{checkSamplingTimeConstraintsForMetaheuristic()}. Discrete algorithms
#' instead enumerate combinations with \code{generateSamplingTimesCombination()}.
#' @param outcome Character string: outcome name.
#' @param initialSamplings Initial sampling-time vector. Continuous
#'   optimizers start from this vector (not from \code{Arm} \code{samplingTimes});
#'   a warning is issued when the two differ. When \code{samplingsWindows} is
#'   set, \code{initialSamplings} must already satisfy the window counts and
#'   \code{minSampling}.
#' @param fixedTimes Times that must remain fixed during optimization.
#' @param numberOfsamplingsOptimisable Total protocol size including
#'   \code{fixedTimes}, not the number of free times
#'   (\code{k = numberOfsamplingsOptimisable - length(fixedTimes)}).
#'   Historical property name (\code{s} lowercase);
#'   \code{numberOfSamplingsOptimisable} is accepted as an alias.
#' @param numberOfSamplingsOptimisable Alias of \code{numberOfsamplingsOptimisable}.
#' @param samplingsWindows List of time windows, in increasing time order.
#' @param numberOfTimesByWindows Number of samples per window.
#' @param minSampling Minimum spacing between consecutive samples.
#' @return An S7 object of class \code{SamplingTimeConstraints}.
#' @export

SamplingTimeConstraints = new_class( "SamplingTimeConstraints", package = "PFIM",

                                     properties = list( outcome = new_property(class_character, default = character(0)),
                                                        initialSamplings = new_property(class_vector, default = numeric(0)),
                                                        fixedTimes = new_property(class_vector, default = numeric(0)),
                                                        numberOfsamplingsOptimisable = new_property(class_double, default = 0.0),
                                                        samplingsWindows = new_property(class_list, default = list()),
                                                        numberOfTimesByWindows = new_property(class_vector, default = numeric(0)),
                                                        minSampling = new_property(class_vector, default = numeric(0))),
                                     validator = function( self ) {
                                       o = prop( self, "outcome" )
                                       if ( length( o ) > 1L )
                                         return( "SamplingTimeConstraints: outcome must be a single string." )
                                       if ( .pfimIsBlankScalar( o ) )
                                         return( "SamplingTimeConstraints: outcome must be a non-empty string." )
                                       if ( any( prop( self, "minSampling" ) < 0, na.rm = TRUE ) ||
                                            any( is.na( prop( self, "minSampling" ) ) ) )
                                         return( "SamplingTimeConstraints: minSampling must be non-negative." )
                                       n = prop( self, "numberOfsamplingsOptimisable" )
                                       if ( length( n ) && ( is.na( n ) || n < 0 ) )
                                         return( "SamplingTimeConstraints: numberOfsamplingsOptimisable must be non-negative." )
                                       init = prop( self, "initialSamplings" )
                                       if ( length( n ) && length( init ) && n > length( init ) )
                                         return( paste0(
                                           "SamplingTimeConstraints: numberOfsamplingsOptimisable (",
                                           n, ") cannot exceed length(initialSamplings) (", length( init ), ")."
                                         ) )
                                       nw = length( prop( self, "samplingsWindows" ) )
                                       nt = length( prop( self, "numberOfTimesByWindows" ) )
                                       ms = length( prop( self, "minSampling" ) )
                                       if ( nw && nt > 1L && nt != nw )
                                         return( "SamplingTimeConstraints: numberOfTimesByWindows length must be 1 or match samplingsWindows." )
                                       if ( nw && ms > 1L && ms != nw )
                                         return( "SamplingTimeConstraints: minSampling length must be 1 or match samplingsWindows." )
                                       orderMsg = .pfimWindowsOrderError( self )
                                       if ( !is.null( orderMsg ) )
                                         return( orderMsg )
                                       winMsg = .pfimInitialSamplingsWindowError( self )
                                       if ( !is.null( winMsg ) )
                                         return( winMsg )
                                       NULL
                                     },
                                     constructor = function( outcome = character(0),
                                                             initialSamplings = numeric(0),
                                                             fixedTimes = numeric(0),
                                                             numberOfsamplingsOptimisable = 0.0,
                                                             numberOfSamplingsOptimisable = NULL,
                                                             samplingsWindows = list(),
                                                             numberOfTimesByWindows = numeric(0),
                                                             minSampling = numeric(0) ) {
                                       if ( !is.null( numberOfSamplingsOptimisable ) ) {
                                         if ( !missing( numberOfsamplingsOptimisable ) &&
                                              numberOfsamplingsOptimisable != 0 &&
                                              !identical(
                                                as.numeric( numberOfsamplingsOptimisable ),
                                                as.numeric( numberOfSamplingsOptimisable )
                                              ) )
                                           stop(
                                             "SamplingTimeConstraints: numberOfSamplingsOptimisable and ",
                                             "numberOfsamplingsOptimisable disagree.",
                                             call. = FALSE
                                           )
                                         numberOfsamplingsOptimisable = numberOfSamplingsOptimisable
                                       }
                                       new_object(
                                         S7_object(),
                                         outcome                      = outcome,
                                         initialSamplings             = initialSamplings,
                                         fixedTimes                   = fixedTimes,
                                         numberOfsamplingsOptimisable = numberOfsamplingsOptimisable,
                                         samplingsWindows             = samplingsWindows,
                                         numberOfTimesByWindows       = numberOfTimesByWindows,
                                         minSampling                  = minSampling
                                       )
                                     } )

#' Draw random sampling times inside each window (respecting min spacing).
#' @param samplingTimeConstraints A \code{SamplingTimeConstraints} object.
#' @param ... Optional method arguments (unused by current methods).
#' @name generateSamplingsFromSamplingConstraints
#' @keywords internal
generateSamplingsFromSamplingConstraints = new_generic( "generateSamplingsFromSamplingConstraints", c( "samplingTimeConstraints" ) )

#' Check that proposed times meet window counts and min spacing.
#' @param samplingTimesConstraints A \code{SamplingTimeConstraints} object.
#' @param arm An \code{Arm} providing the candidate sampling schedule.
#' @param ... Method arguments: \code{newSamplings} (numeric times) and
#'   \code{outcome} (outcome name) for the metaheuristic check method.
#' @name checkSamplingTimeConstraintsForMetaheuristic
#' @keywords internal
checkSamplingTimeConstraintsForMetaheuristic = new_generic( "checkSamplingTimeConstraintsForMetaheuristic", c( "samplingTimesConstraints", "arm"  ) )

#' Build window table: min, max, delta (= minSampling), n (= times per window).
#' @noRd
#' @keywords internal
.pfimSamplingWindowsTable = function( samplingTimeConstraints ) {
  windows = prop( samplingTimeConstraints, "samplingsWindows" )
  tab = as.data.frame( do.call(
    "cbind",
    list(
      t( as.data.frame( windows ) ),
      as.data.frame( list(
        prop( samplingTimeConstraints, "minSampling" ),
        prop( samplingTimeConstraints, "numberOfTimesByWindows" )
      ) )
    )
  ) )
  colnames( tab ) = c( "min", "max", "delta", "n" )
  rownames( tab ) = NULL
  tab
}

#' Error if sampling windows are not declared in increasing time order.
#' @noRd
#' @keywords internal
.pfimWindowsOrderError = function( samplingTimeConstraints ) {
  win = prop( samplingTimeConstraints, "samplingsWindows" )
  if ( length( win ) < 2L )
    return( NULL )
  mins = map_dbl( win, function( w ) min( unlist( w ) ) )
  if ( any( !is.finite( mins ) ) )
    return( "SamplingTimeConstraints: samplingsWindows must have finite bounds." )
  if ( is.unsorted( mins, strictly = TRUE ) )
    return( "SamplingTimeConstraints: samplingsWindows must be declared in increasing order." )
  NULL
}

#' Split sorted times into the declared per-window counts, then check that
#' each block lies in its window.
#'
#' Contiguous windows such as \code{[0,4]}+\code{[4,12]} are a partition of the
#' protocol: the first \code{n[1]} times belong to window 1, the next
#' \code{n[2]} to window 2, and a shared endpoint can sit on either side
#' according to those counts. Times outside their assigned window fail.
#' @return List of numeric blocks (one per window), or \code{NULL} if the
#'   split is impossible (length mismatch).
#' @noRd
#' @keywords internal
.pfimTimesByDeclaredWindows = function( times, n_by, nw ) {
  times = as.numeric( times )
  if ( length( n_by ) == 1L && nw > 1L )
    n_by = rep( n_by, nw )
  if ( length( n_by ) != nw )
    return( NULL )
  n_by = as.integer( n_by )
  if ( any( n_by < 0L ) || length( times ) != sum( n_by ) )
    return( NULL )
  st = sort( times )
  groups = vector( "list", nw )
  pos = 0L
  for ( i in seq_len( nw ) ) {
    k = n_by[[ i ]]
    if ( k == 0L ) {
      groups[[ i ]] = numeric( 0L )
      next
    }
    groups[[ i ]] = st[ pos + seq_len( k ) ]
    pos = pos + k
  }
  groups
}

#' Whether each declared-count block lies inside its sampling window.
#' @noRd
#' @keywords internal
.pfimWindowBlocksInside = function( groups, mins, maxs ) {
  if ( is.null( groups ) )
    return( FALSE )
  all( vapply( seq_along( groups ), function( i ) {
    x = groups[[ i ]]
    length( x ) == 0L || all( x >= mins[[ i ]] & x <= maxs[[ i ]] )
  }, logical( 1L ) ) )
}

#' Error string when sampling times violate declared windows; else \code{NULL}.
#'
#' Skipped when windows or \code{times} are empty (discrete constraints, partial
#' objects). Occupancy splits the sorted times by \code{numberOfTimesByWindows}
#' (so a shared endpoint is not double-counted). \code{minSampling} uses the
#' same blocks.
#' @noRd
#' @keywords internal
.pfimTimesWindowError = function( samplingTimeConstraints, times,
                                  label = "initialSamplings" ) {
  win  = prop( samplingTimeConstraints, "samplingsWindows" )
  init = as.numeric( times )
  if ( !length( win ) || !length( init ) )
    return( NULL )
  n_by = prop( samplingTimeConstraints, "numberOfTimesByWindows" )
  if ( !length( n_by ) )
    return( paste0(
      "SamplingTimeConstraints: numberOfTimesByWindows is required when ",
      label, " and samplingsWindows are set."
    ) )
  nw = length( win )
  if ( length( n_by ) == 1L && nw > 1L )
    n_by = rep( n_by, nw )
  if ( length( n_by ) != nw )
    return( NULL )
  expected = sum( n_by )
  if ( expected != length( init ) )
    return( paste0(
      "SamplingTimeConstraints: length(", label, ") (", length( init ),
      ") must equal sum(numberOfTimesByWindows) (", expected, ")."
    ) )
  bounds = map( win, unlist )
  mins = map_dbl( bounds, min )
  maxs = map_dbl( bounds, max )
  groups = .pfimTimesByDeclaredWindows( init, n_by, nw )
  if ( is.null( groups ) || !.pfimWindowBlocksInside( groups, mins, maxs ) ) {
    outside = is.null( groups ) || !all( vapply( init, function( t ) {
      any( t >= mins & t <= maxs )
    }, logical( 1L ) ) )
    if ( outside )
      return( paste0(
        "SamplingTimeConstraints: ", label,
        " fall outside every sampling window."
      ) )
    got = if ( is.null( groups ) ) {
      rep( 0L, nw )
    } else {
      vapply( seq_len( nw ), function( i ) {
        x = groups[[ i ]]
        sum( x >= mins[[ i ]] & x <= maxs[[ i ]] )
      }, integer( 1L ) )
    }
    return( paste0(
      "SamplingTimeConstraints: ", label, " do not place ",
      paste( n_by, collapse = ", " ),
      " time(s) in each sampling window (got ",
      paste( got, collapse = ", " ), ")."
    ) )
  }
  delta = prop( samplingTimeConstraints, "minSampling" )
  if ( !length( delta ) )
    return( NULL )
  if ( length( delta ) == 1L )
    delta = rep( delta, nw )
  for ( i in seq_len( nw ) ) {
    x = groups[[ i ]]
    if ( length( x ) < 2L ) next
    if ( any( diff( x ) < delta[[ i ]] ) )
      return( paste0(
        "SamplingTimeConstraints: ", label, " violate minSampling inside a window."
      ) )
  }
  NULL
}

#' @noRd
#' @keywords internal
.pfimInitialSamplingsWindowError = function( samplingTimeConstraints ) {
  .pfimTimesWindowError(
    samplingTimeConstraints,
    prop( samplingTimeConstraints, "initialSamplings" )
  )
}

#' Draw feasible sampling times under window and spacing constraints.
#'
#' For each window \code{[min, max]} with \code{n} points and minimum gap
#' \code{delta}, places ordered uniforms so consecutive samples are at least
#' \code{delta} apart: remaining slack is \code{max - min - (n-1)*delta}.
#' @name generateSamplingsFromSamplingConstraints
#' @return Numeric vector of generated sampling times (all windows concatenated).
#' @usage NULL
#' @keywords internal

method( generateSamplingsFromSamplingConstraints, SamplingTimeConstraints ) = function( samplingTimeConstraints ) {
  .samplingsFromWindowsTable( samplingTimeConstraints, random = TRUE )
}

#' Place \code{n} times in each window using the same slack formula as the
#' random generator. Mid-bin quantiles avoid shared window endpoints.
#' @noRd
#' @keywords internal
.evenSamplingsFromSamplingConstraints = function( samplingTimeConstraints ) {
  .samplingsFromWindowsTable( samplingTimeConstraints, random = FALSE )
}

#' @noRd
#' @keywords internal
.samplingsFromWindowsTable = function( samplingTimeConstraints, random ) {
  tab = .pfimSamplingWindowsTable( samplingTimeConstraints )
  unlist( pmap( tab, function( min, max, delta, n ) {
    n = as.integer( n )
    if ( n <= 0L ) return( numeric( 0L ) )
    distance = max - min - ( n - 1 ) * delta
    u = if ( isTRUE( random ) )
      sort( stats::runif( n, 0, 1 ) )
    else
      ( seq_len( n ) - 0.5 ) / n
    min + distance * u + delta * seq( 0, n - 1 )
  } ), use.names = FALSE )
}

#' Check sampling feasibility for PSO, PGBO or simplex.
#'
#' Returns two logicals: whether each window contains exactly
#' \code{numberOfTimesByWindows} points, and whether consecutive points inside
#' each window respect \code{minSampling}. Shared endpoints are partitioned by
#' \code{numberOfTimesByWindows} (same split as \code{.pfimTimesByDeclaredWindows()}).
#' @name checkSamplingTimeConstraintsForMetaheuristic
#' @return List with \code{constraintWindowsLength} and \code{constraintMinimalSampling}.
#' @usage NULL
#' @keywords internal

method( checkSamplingTimeConstraintsForMetaheuristic, list( SamplingTimeConstraints, Arm ) ) = function( samplingTimesConstraints, arm, newSamplings, outcome ) {
  tab = .pfimSamplingWindowsTable( samplingTimesConstraints )
  groups = .pfimTimesByDeclaredWindows( newSamplings, tab$n, nrow( tab ) )
  countsOk = .pfimWindowBlocksInside( groups, tab$min, tab$max )
  gapsOk = FALSE
  if ( countsOk ) {
    gapsOk = all( vapply( seq_along( groups ), function( i ) {
      x = groups[[ i ]]
      if ( length( x ) < 2L ) return( TRUE )
      all( diff( x ) >= tab$delta[[ i ]] )
    }, logical( 1L ) ) )
  }
  list(
    constraintWindowsLength   = countsOk,
    constraintMinimalSampling = gapsOk
  )
}
