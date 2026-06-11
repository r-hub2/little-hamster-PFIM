# End-to-end optimization smoke tests for all algorithms.

test_that("FedorovWynnAlgorithm: run() completes with valid optimal design", {
  opt = .minimal_discrete_opt( "FedorovWynnAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("PSOAlgorithm: run() completes with valid optimal design", {
  opt = .minimal_continuous_opt( "PSOAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("PGBOAlgorithm: run() completes with valid optimal design", {
  opt = .minimal_continuous_opt( "PGBOAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("SimplexAlgorithm: run() completes with valid optimal design", {
  opt = .minimal_continuous_opt( "SimplexAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("MultiplicativeAlgorithm: run() completes with valid optimal design", {
  opt = .minimal_discrete_opt( "MultiplicativeAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("generateFimsFromConstraints returns constraint-grid FIMs", {
  opt = .minimal_discrete_opt( "MultiplicativeAlgorithm" )
  res = generateFimsFromConstraints( opt )
  expect_type( res, "list" )
  expect_true( length( res$listFimsAlgoMult ) >= 1L )
  expect_gte( nrow( res$listFimsAlgoMult[[ 1L ]][[ 1L ]] ), 1L )
})
