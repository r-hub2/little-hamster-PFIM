# =============================================================================
# Residual variance diagonals + dV/dsigma (Combined1 / Combined2)
# =============================================================================
#
# Primary path: Rcpp Armadillo (`residualErrorDerivatives_Rcpp`) - vectorized,
# no R allocations in the inner loop. R fallback keeps source / uncompiled
# installs working (CRAN checks, load_all without compile).
#
# Forms (string or mapped to C++ form_id 1/2):
#   combined1 : V = (a + b f^c)^2 ,  dV/da = 2 base , dV/db = 2 base f^c
#   combined2 : V = a^2 + (b f^c)^2, dV/da = 2 a    , dV/db = 2 (b f^c) f^c
#
# Estimable derivatives are omitted when !\.pfimSigmaIsEstimable(value, fixed)
# so FIM column layout matches .sigmaNames() / .sigmaValues().
#
# Return shape: Matrix::Diagonal for V and each dV/dsigma - multi-output
# bdiag stays sparse until FIM kernels densify with as.matrix().
# =============================================================================

#' R fallback for residual variance / sigma derivatives (same math as C++).
#' @noRd
#' @keywords internal
.pfimErrorModelDerivativesR = function( sigmaInter, sigmaSlope, sigmaInterFixed,
                                        sigmaSlopeFixed, cError, f,
                                        form = "combined1" ) {
  f = as.numeric( f )
  n = length( f )
  if ( !n )
    stop( "residual error: prediction vector f has length 0.", call. = FALSE )
  if ( !( cError > 0 ) )
    stop( "residual error: cError must be positive.", call. = FALSE )

  # f^c once; cError == 1 is the common case (skip pow).
  f_term = if ( cError == 1 ) f else f^cError

  if ( identical( form, "combined2" ) ) {
    prop_sd  = sigmaSlope * f_term
    variance = sigmaInter^2 + prop_sd^2
    d_inter  = rep( 2 * sigmaInter, n )
    d_slope  = 2 * prop_sd * f_term
  } else {
    base     = sigmaInter + sigmaSlope * f_term
    variance = base^2
    d_inter  = 2 * base
    d_slope  = 2 * base * f_term
  }

  list(
    variance = variance,
    d_inter  = if ( .pfimSigmaIsEstimable( sigmaInter, sigmaInterFixed ) )
      d_inter else numeric( 0 ),
    d_slope  = if ( .pfimSigmaIsEstimable( sigmaSlope, sigmaSlopeFixed ) )
      d_slope else numeric( 0 )
  )
}

#' Wrap a numeric vector as \code{Matrix::Diagonal} (empty -> 0x0).
#' @noRd
#' @keywords internal
.pfimDiagonalFromVec = function( x ) {
  x = as.numeric( x )
  if ( !length( x ) ) return( Matrix::Diagonal( n = 0L ) )
  Matrix::Diagonal( x = x )
}

#' Embed a diagonal block into a sparse \code{total x total} matrix.
#'
#' Used when stacking multi-output observation orderings: each outcome's
#' dV/dsigma lives on its sampling-time slice without densifying zeros.
#' @param d Numeric diagonal, or square diagonal matrix.
#' @param offset 0-based start index in the full observation vector.
#' @param total Full observation dimension.
#' @noRd
#' @keywords internal
.pfimEmbedDiagonalBlock = function( d, offset, total ) {
  if ( is.matrix( d ) || inherits( d, "Matrix" ) )
    d = diag( as.matrix( d ) )
  d = as.numeric( d )
  n = length( d )
  if ( !n || !total )
    return( Matrix::Diagonal( n = as.integer( total ) ) )
  i = as.integer( offset ) + seq_len( n )
  Matrix::sparseMatrix(
    i = i, j = i, x = d,
    dims = c( as.integer( total ), as.integer( total ) )
  )
}

#' Residual V and estimable dV/dsigma diagonals (Rcpp with R fallback).
#' @noRd
#' @keywords internal
.pfimErrorModelDerivatives = function( sigmaInter, sigmaSlope, sigmaInterFixed,
                                       sigmaSlopeFixed, cError, f,
                                       form = "combined1" ) {
  form_id = if ( identical( form, "combined2" ) ) 2L else 1L
  # tryCatch: missing DLL / older builds still evaluate correctly in R.
  out = tryCatch(
    residualErrorDerivatives_Rcpp(
      sigmaInter, sigmaSlope, sigmaInterFixed, sigmaSlopeFixed,
      cError, as.numeric( f ), form_id
    ),
    error = function( e ) {
      .pfimErrorModelDerivativesR(
        sigmaInter, sigmaSlope, sigmaInterFixed, sigmaSlopeFixed,
        cError, f, form
      )
    }
  )

  # Drop empty derivative vectors; wrap survivors as Diagonal for sparse bdiag.
  sigmaDerivatives = list(
    sigmaInter = out$d_inter,
    sigmaSlope = out$d_slope
  ) |>
    keep( ~ length( .x ) > 0L ) |>
    map( .pfimDiagonalFromVec )

  list(
    sigmaDerivatives = sigmaDerivatives,
    errorVariance    = .pfimDiagonalFromVec( out$variance )
  )
}
