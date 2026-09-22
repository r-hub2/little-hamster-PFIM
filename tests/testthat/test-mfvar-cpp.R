# Variance-FIM trace block: C++ kernel vs R reference.

test_that( "computeMFVar_Rcpp matches R on random SPD system", {
  set.seed( 42L )
  p = 6L
  V = crossprod( matrix( rnorm( p * p ), p, p ) ) + diag( p )
  V_inv = PFIM:::.safeCholInv( V )
  dV_list = replicate(
    4L,
    {
      M = matrix( rnorm( p * p ), p, p )
      M + t( M )
    },
    simplify = FALSE
  )

  ref = PFIM:::.computeMFVar_R( V_inv, dV_list )
  out = computeMFVar_Rcpp( V_inv, dV_list )
  expect_equal( out, ref, tolerance = 1e-10 )
  expect_true( isSymmetric( out, tol = 1e-9 ) )
} )

test_that( "computeMFVar_mixed_Rcpp matches dense kernel on rank-1 columns", {
  set.seed( 43L )
  p = 20L
  n = 5L
  V = crossprod( matrix( rnorm( p * p ), p, p ) ) + diag( p )
  V_inv = PFIM:::.safeCholInv( V )
  W = matrix( rnorm( p * n ), p, n )
  dV_list = lapply( seq_len( n ), function( i ) tcrossprod( W[ , i ] ) )

  ref = computeMFVar_Rcpp( V_inv, dV_list )
  out = computeMFVar_mixed_Rcpp( V_inv, W, list() )
  expect_equal( out, ref, tolerance = 1e-10 )
} )

test_that( "computeMFVar_Rcpp rejects wrong matrix dimensions", {
  V_inv = diag( 3 )
  bad = list( matrix( 1, 2, 2 ) )
  expect_error( computeMFVar_Rcpp( V_inv, bad ), "not 3x3" )
} )

test_that( "chol_inv rejects non-finite input", {
  M = matrix( c( 1, NA, NA, 1 ), 2, 2 )
  expect_no_warning(
    expect_error( chol_inv_Rcpp( M ), "finite|not positive definite" )
  )
} )

test_that( "computeMFVar_mixed_Rcpp matches dense kernel on mixed rank-1 and full", {
  set.seed( 44L )
  p = 12L
  V = crossprod( matrix( rnorm( p * p ), p, p ) ) + diag( p )
  V_inv = PFIM:::.safeCholInv( V )
  W = matrix( rnorm( p * 3L ), p, 3L )
  full = replicate( 2L, {
    M = matrix( rnorm( p * p ), p, p )
    M + t( M )
  }, simplify = FALSE )
  dV_list = c(
    lapply( seq_len( ncol( W ) ), function( i ) tcrossprod( W[ , i ] ) ),
    full
  )

  ref = computeMFVar_Rcpp( V_inv, dV_list )
  out = computeMFVar_mixed_Rcpp( V_inv, W, full )
  expect_equal( out, ref, tolerance = 1e-10 )
} )
