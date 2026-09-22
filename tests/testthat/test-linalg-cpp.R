# C++ linalg kernels vs R references.

test_that( "chol_inv_Rcpp matches chol2inv reference", {
  cases = list(
    .random_spd( 4L, 1L ),
    .random_spd( 12L, 2L ),
    diag( 8 ) + 1e-12
  )
  purrr::walk( cases, function( V ) {
    expect_equal( chol_inv_Rcpp( V ), .safeCholInv_ref( V ), tolerance = 1e-10 )
  } )
} )

test_that( "chol_inv on cas1 population V block", {
  cas = Filter( function( x ) x$id == "cas1_noCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  V   = .pop_fim_simple_V( run( build_cov_iov_evaluation( cas ) ) )
  expect_equal( chol_inv_Rcpp( V ), .safeCholInv_ref( V ), tolerance = 1e-9 )
} )

test_that( "chol_inv fails on indefinite matrix", {
  M = matrix( c( 1, 2, 2, 4 ), 2, 2 )
  expect_error( chol_inv_Rcpp( M ), "not positive definite" )
} )

test_that( "chol_inv fails on the all-zero matrix", {
  V = matrix( 0, 2, 2 )
  expect_error( chol_inv_Rcpp( V ), "zero or not positive definite|not positive definite" )
} )

test_that( ".packFisherLowerTriangle matches Fedorov PackedFim order", {
  M = matrix( c( 1, 2, 3, 2, 4, 5, 3, 5, 6 ), 3, 3 )
  packed = as.numeric( PFIM:::.packFisherLowerTriangle( M ) )
  expect_equal( packed, c( 1, 2, 4, 3, 5, 6 ) )
  legacy = M[ rev( lower.tri( t( M ), diag = TRUE ) ) ]
  expect_equal( packed, legacy )
} )

test_that( "safe_solve_Rcpp matches solve reference", {
  set.seed( 7L )
  cases = list(
    matrix( c( 2, 1, 1, 3 ), 2, 2 ),
    crossprod( matrix( rnorm( 25 ), 5, 5 ) ),
    diag( 4 ) * c( 1, 1e-14, 1, 1e-14 )
  )
  purrr::walk( cases, function( M ) {
    expect_equal( safe_solve_Rcpp( M ), .safeSolve_ref( M ), tolerance = 1e-9 )
  } )
} )

test_that( "safe_solve fails on singular matrix", {
  M = matrix( c( 1, 2, 2, 4 ), 2, 2 )
  expect_error( safe_solve_Rcpp( M ), "singular or ill-conditioned" )
} )

test_that( ".safeCholInv and .safeSolve call Rcpp kernels", {
  V = .random_spd( 5L, 4L )
  M = matrix( c( 3, 1, 0.5, 2 ), 2, 2 )
  expect_equal( PFIM:::.safeCholInv( V ), chol_inv_Rcpp( V ) )
  expect_equal( PFIM:::.safeSolve( M ), safe_solve_Rcpp( M ) )
} )
