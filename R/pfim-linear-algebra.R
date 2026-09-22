# Matrix inverses and FIM packing helpers (C++ kernels in pfim-linalg.hpp).
#
# Thin R wrappers around chol_inv / safe_solve, plus lower-triangle packing used
# by the Fedorov-Wynn PackedFim layout, and correlation from the FIM inverse.

#' Cholesky-based inverse of a symmetric positive-definite matrix.
#' @noRd
#' @keywords internal
.safeCholInv = function( V ) chol_inv_Rcpp( V )

#' Pack the lower triangle of a symmetric FIM row-wise (Fedorov-Wynn \code{PackedFim} order).
#'
#' Returns a 1-row matrix of length \eqn{p(p+1)/2} with elements
#' \eqn{M_{11}, M_{21}, M_{22}, M_{31}, \ldots}.
#' @noRd
#' @keywords internal
.packFisherLowerTriangle = function( M ) {
  # Row-wise lower triangle of a symmetric M equals column-major upper.tri.
  matrix( M[ upper.tri( M, diag = TRUE ) ], nrow = 1L )
}

#' Estimator correlation matrix from a Fisher information matrix.
#'
#' Computes \code{cov2cor(solve(M))} via a safe Cholesky inverse. Singular FIMs
#' return an NA matrix instead of aborting (same cases as \code{.fimBuildSeAndRse}).
#' @noRd
#' @keywords internal
.fimCorrelationMatrix = function( M ) {
  tryCatch(
    stats::cov2cor( .safeCholInv( M ) ),
    error = function( e ) {
      p = nrow( M )
      matrix( NA_real_, p, p, dimnames = dimnames( M ) )
    }
  )
}

#' General square-matrix inverse (singular -> error from C++).
#' @noRd
#' @keywords internal
.safeSolve = function( M ) safe_solve_Rcpp( M )

#' Relative eigenvalue cutoff on the correlation scale (unit-invariant rank).
#'
#' Machine-\eqn{\varepsilon\cdot p\cdot\lambda_{\max}} on the raw FIM treats
#' roundoff of a large block as information (finite SE instead of Inf).
#' @noRd
#' @keywords internal
.pfimPsdRelTol = 1e-10

#' Correlation-scale eigen-decomposition of a symmetric matrix.
#'
#' Forms \eqn{R=D^{-1/2}MD^{-1/2}} with \eqn{D=\mathrm{diag}(M)_+} so the rank
#' decision does not depend on parameter units. Eigenvalues of \eqn{R} below
#' \code{relTol * max(lambda)} are a numerical null. Non-positive diagonal
#' entries get scale 0 (no information on that parameter).
#' @noRd
#' @keywords internal
.pfimPsdCorrelationEigen = function( M, relTol = .pfimPsdRelTol ) {
  M = as.matrix( 0.5 * ( M + t( M ) ) )
  p = nrow( M )
  if ( !p )
    return( list(
      p = 0L, M = M, s = numeric( 0L ),
      values = numeric( 0L ), vectors = matrix( 0, 0L, 0L ),
      pos = logical( 0L ), pos_diag = logical( 0L )
    ) )
  if ( any( !is.finite( M ) ) )
    return( list(
      p = p, M = M, s = rep( 0, p ),
      values = rep( 0, p ), vectors = diag( p ),
      pos = rep( FALSE, p ), pos_diag = rep( FALSE, p ),
      nonfinite = TRUE
    ) )

  d = diag( M )
  pos_diag = is.finite( d ) & d > 0
  s = numeric( p )
  s[ pos_diag ] = 1 / sqrt( d[ pos_diag ] )
  R = M * tcrossprod( s )

  if ( !any( pos_diag ) ) {
    return( list(
      p = p, M = M, s = s, R = R,
      values = rep( 0, p ), vectors = diag( p ),
      pos = rep( FALSE, p ), pos_diag = pos_diag
    ) )
  }

  ev    = eigen( R, symmetric = TRUE )
  scale = max( ev$values, 0 )
  tol   = scale * relTol
  list(
    p = p, M = M, s = s, R = R,
    values = ev$values, vectors = ev$vectors,
    pos = ev$values > tol, pos_diag = pos_diag, tol = tol
  )
}

#' Zero off-diagonal and set Inf variance on selected coordinates.
#' @noRd
#' @keywords internal
.pfimZeroCrossInfDiag = function( C, idx ) {
  if ( !length( idx ) )
    return( C )
  C[ idx, ] = 0
  C[ , idx ] = 0
  C[ cbind( idx, idx ) ] = Inf
  C
}

#' Moore-Penrose inverse with Inf variance on the correlation-scale null space.
#'
#' Identifiable directions: \eqn{C=D^{-1/2}R^+D^{-1/2}}. Parameters that
#' participate in a null eigenvector, or that have a non-positive diagonal,
#' get \code{Inf} variance and zero cross-covariances.
#' @noRd
#' @keywords internal
.pfimPsdPseudoInverse = function( M, relTol = .pfimPsdRelTol ) {
  eg = .pfimPsdCorrelationEigen( M, relTol )
  p  = eg$p
  if ( !p )
    return( matrix( numeric( 0L ), 0L, 0L ) )
  if ( isTRUE( eg$nonfinite ) )
    return( matrix( Inf, p, p ) )

  C = matrix( 0, p, p )
  if ( any( eg$pos ) ) {
    Vpos = eg$vectors[ , eg$pos, drop = FALSE ]
    Rinv = Vpos %*% ( t( Vpos ) / eg$values[ eg$pos ] )
    C = Rinv * tcrossprod( eg$s )
  }
  if ( any( !eg$pos ) ) {
    V0 = eg$vectors[ , !eg$pos, drop = FALSE ]
    null_mass = as.numeric( rowSums( V0^2 ) )
    C = .pfimZeroCrossInfDiag( C, which( null_mass > 1e-8 | !eg$pos_diag ) )
  }
  0.5 * ( C + t( C ) )
}
