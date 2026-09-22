# Residual error variance forms (Combined1 / Combined2).

test_that("Combined2 uses quadrature variance (PopED-style)", {
  err = Combined2( output = "RespPK", sigmaInter = 0.2, sigmaSlope = 0.3, cError = 1 )
  f   = c( 1, 2 )
  res = PFIM:::evaluateErrorModelDerivatives( err, f )
  expect_equal( diag( as.matrix( res$errorVariance ) ), 0.2^2 + ( 0.3 * f )^2 )
  expect_equal( as.matrix( res$sigmaDerivatives$sigmaInter )[ 1L, 1L ], 2 * 0.2 )
  expect_equal( as.matrix( res$sigmaDerivatives$sigmaSlope )[ 2L, 2L ], 2 * 0.3 * 4 )
})

test_that("Combined2 respects cError and n = 1", {
  err = Combined2( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.1, cError = 2 )
  res = PFIM:::evaluateErrorModelDerivatives( err, 3 )
  expect_equal( as.matrix( res$errorVariance ), matrix( 0.5^2 + ( 0.1 * 9 )^2, 1L, 1L ) )
  expect_equal( as.matrix( res$sigmaDerivatives$sigmaSlope )[ 1L, 1L ], 2 * 0.1 * 81 )
})

test_that("Combined2 drops fixed or zero sigma derivatives", {
  err = Combined2(
    output = "RespPK", sigmaInter = 0, sigmaSlope = 0.2,
    sigmaInterFixed = TRUE
  )
  res = PFIM:::evaluateErrorModelDerivatives( err, c( 1, 2 ) )
  expect_null( res$sigmaDerivatives$sigmaInter )
  expect_named( res$sigmaDerivatives, "sigmaSlope" )
})

test_that("cError applies power to Combined1 proportional term", {
  err = Combined1( output = "RespPK", sigmaInter = 0.2, sigmaSlope = 0.3, cError = 2 )
  f   = c( 1, 2 )
  res = PFIM:::evaluateErrorModelDerivatives( err, f )
  expect_equal( diag( as.matrix( res$errorVariance ) ), ( 0.2 + 0.3 * f^2 )^2 )
  expect_equal( as.matrix( res$sigmaDerivatives$sigmaInter )[ 1L, 1L ], 2 * ( 0.2 + 0.3 ) )
  expect_equal( as.matrix( res$sigmaDerivatives$sigmaSlope )[ 2L, 2L ], 2 * ( 0.2 + 0.3 * 4 ) * 4 )
})

test_that("cError = 1 matches linear Combined1 proportional term", {
  f = c( 1.5, 2 )
  err_c = Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.2, cError = 1 )
  err_l = Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.2 )
  expect_equal(
    PFIM:::evaluateErrorModelDerivatives( err_c, f ),
    PFIM:::evaluateErrorModelDerivatives( err_l, f )
  )
})

test_that("Constant rejects nonzero sigmaSlope", {
  expect_error(
    Constant( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.05 ),
    "sigmaSlope must be 0",
    fixed = TRUE
  )
})

test_that("Proportional rejects nonzero sigmaInter", {
  expect_error(
    Proportional( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.2 ),
    "sigmaInter must be 0",
    fixed = TRUE
  )
})

test_that("getModelErrorData reports cError", {
  err = Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.2, cError = 1.5 )
  row = getModelErrorData( err )
  expect_equal( row$cError, "1.5" )
})
