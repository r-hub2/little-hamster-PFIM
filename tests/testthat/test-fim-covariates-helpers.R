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

test_that("numberOfOccasions: NA infers; explicit value must match covariates/IOV", {
  fx = cas10_design()
  base = list(
    name                    = "occ_test",
    modelFromLibrary        = cas10_model_from_library(),
    modelParameters         = cas10_model_parameters(),
    modelCovariates         = cas10_covariates(),
    modelCovariatesEquation = "exponential",
    modelError              = fx$modelError,
    designs                 = list( fx$design1 ),
    fimType                 = "individual",
    outputs                 = list( "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
  ev_infer = do.call( Evaluation, c( base, list( numberOfOccasions = NA_real_ ) ) )
  m_infer  = defineModelType( ev_infer )
  expect_equal( prop( m_infer, "numberOfOccasions" ), 4L )

  ev_ok = do.call( Evaluation, c( base, list( numberOfOccasions = 4L ) ) )
  m_ok  = defineModelType( ev_ok )
  expect_equal( prop( m_ok, "numberOfOccasions" ), 4L )

  ev_bad = do.call( Evaluation, c( base, list( numberOfOccasions = 2L ) ) )
  expect_error(
    defineModelType( ev_bad ),
    "numberOfOccasions \\(2\\) is inconsistent with covariates/IOV \\(expected 4\\)"
  )
})

test_that("gamma-only IOV: numberOfOccasions >= 2 accepted", {
  args = list(
    name                    = "gamma_iov_occ",
    modelParameters         = cov_iov_parameters( with_iov = TRUE ),
    modelCovariates         = list(),
    modelCovariatesEquation = "exponential",
    modelEquations          = cov_iov_model_equations(),
    modelError              = cov_iov_model_error(),
    designs                 = list( cov_iov_evaluation_design() ),
    fimType                 = "population",
    outputs                 = list( RespPK = "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 ),
    numberOfOccasions       = 3L
  )
  ev = do.call( Evaluation, args )
  m  = defineModelType( ev )
  expect_equal( prop( m, "numberOfOccasions" ), 3L )

  args$numberOfOccasions = 1L
  expect_error(
    do.call( Evaluation, args ) |> defineModelType(),
    "numberOfOccasions must be >= 2 when gamma"
  )
})

test_that("cas10 individual: run, FIM 4x4 (mu+sigma; no beta), covariateTest empty betas", {
  evaluation = cas10_evaluation( "cas10_ind", "individual" )
  evaluation = run( evaluation )

  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( ncol( prop( fim, "fisherMatrix" ) ), 4L )

  ct  = covariateTest( evaluation )
  sig = prop( ct, "covariate" )
  # Individual FIM omits beta; significance table has no covariate rows.
  expect_equal( nrow( sig ), 0L )
})

test_that("cas10 Bayesian: run, FIM 3x3 (mu only), covariateTest empty betas", {
  evaluation = cas10_evaluation( "cas10_bayes", "Bayesian" )
  evaluation = run( evaluation )

  fim_b = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( ncol( prop( fim_b, "fisherMatrix" ) ), 3L )
  se = prop( fim_b, "SEAndRSE" )$SEAndRSE
  expect_equal( nrow( se ), 3L )
  ct = covariateTest( evaluation )
  expect_equal( nrow( prop( ct, "covariate" ) ), 0L )
})

test_that( "Bayesian FIM with covariates omits beta (mu block only)", {
  evaluation = run( cas10_evaluation( "cas10_bayes_cross", "Bayesian" ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M   = prop( fim, "fisherMatrix" )
  fe  = PFIM:::.fimFixedEffectLabels( evaluation )
  n_mu = length( fe$columnNamesMu )
  expect_equal( ncol( M ), n_mu )
  expect_false( any( grepl( "beta|β", colnames( M ) ) ) )
})

test_that( "rebuildEvalModel returns independent copies from cache", {
  ev = cas10_evaluation( "rebuild_clone", "population" )
  PFIM:::.pfimClearFimCaches()
  m1 = rebuildEvalModel( ev )
  m2 = rebuildEvalModel( ev )
  prop( m1, "modelCovariates" ) = list()
  m3 = rebuildEvalModel( ev )
  expect_gt( length( prop( m3, "modelCovariates" ) ), 0L )
  expect_true( usesCovariateOccasionStructure( m3 ) )
} )

test_that( ".fdModelEvaluations clones before stripping covariates", {
  ev  = cas10_evaluation( "fd_strip", "population" )
  arm = prop( prop( ev, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  PFIM:::.pfimClearFimCaches()
  invisible( rebuildEvalModel( ev, finiteDifference = TRUE ) )
  cacheId = PFIM:::.pfimModelCacheId( ev )
  cached  = PFIM:::.pfimEvalModelCache[[ cacheId ]]$fd
  invisible( PFIM:::.fdModelEvaluations( cached, arm ) )
  expect_true( usesCovariateOccasionStructure( cached ) )
  expect_gt( length( prop( cached, "modelCovariates" ) ), 0L )
} )
