# Helpers: rebuildEvalModel, beta dict, cas10-style individual/Bayesian (light)

test_that("rebuildEvalModel sets outputNames", {
  evaluation = cas10_evaluation( "test_rebuild_ind", "individual" )
  m = rebuildEvalModel( evaluation, finiteDifference = FALSE )
  expect_gt( length( prop( m, "outputNames" ) ), 0L )
  expect_equal( prop( m, "outputNames" ), "RespPK" )

  covs = cas10_covariates()
  beta_names = .betaInternalNamesFromCovariates( covs )
  expect_length( beta_names, 2L )
  beta_vals = .betaValuesFromCovariates( covs, beta_names )
  expect_length( beta_vals, 2L )
  expect_equal( beta_vals[ 1L ], log( 1.2 ), tolerance = 1e-8 )
  expect_equal( beta_vals[ 2L ], log( 1.1 ), tolerance = 1e-8 )
})

test_that("cas10 individual: run, FIM 6x6, covariateTest significance on beta only", {
  evaluation = cas10_evaluation( "cas10_ind", "individual" )
  evaluation = run( evaluation )

  fim = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( ncol( prop( fim, "fisherMatrix" ) ), 6L )

  ct  = covariateTest( evaluation )
  sig = prop( ct, "significance" )
  expect_false( any( startsWith( sig$Parameter, "\u03bc_" ) ) )
  expect_equal( nrow( sig ), 2L )
})

test_that("cas10 Bayesian: run, FIM 5x5, covariateTest", {
  evaluation = cas10_evaluation( "cas10_bayes", "Bayesian" )
  evaluation = run( evaluation )

  fim_b = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( ncol( prop( fim_b, "fisherMatrix" ) ), 5L )
  se = prop( fim_b, "SEAndRSE" )$SEAndRSE
  expect_equal( nrow( se ), 5L )
  expect_s3_class( covariateTest( evaluation ), "PFIM::CovariateTest" )
})
