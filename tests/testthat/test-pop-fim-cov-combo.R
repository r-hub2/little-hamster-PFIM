# Population covariate / IOV FIM — C++ vs R reference.

.popFimComboInputs = function( evaluation, iter = 1L, arm = NULL ) {
  model = rebuildEvalModel( evaluation, finiteDifference = FALSE )
  if ( is.null( arm ) )
    arm = prop( prop( evaluation, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]

  parameters = prop( model, "modelParameters" )
  distrib    = map( parameters, ~ prop( .x, "distribution" ) )
  allGrad    = prop( arm, "evaluationGradients" )
  varRes     = prop( arm, "evaluationVariance" )
  outNames   = prop( model, "outputNames" )
  nOcc       = length( varRes[[ iter ]]$variances )

  grads = map( seq_len( nOcc ), function( occ ) {
    gpo = map( outNames, ~ t( allGrad[[ iter ]]$gradients[[ occ ]]$gradient[[ .x ]] ) )
    if ( length( gpo ) == 1L ) gpo[[ 1L ]] else do.call( cbind, gpo )
  } )

  sigma_n = length( varRes[[ iter ]]$variances[[ 1L ]]$variance$sigmaDerivatives )

  list(
    gradients = if ( length( grads ) == 1L ) grads[[ 1L ]] else reduce( grads, cbind ),
    mu_values = PFIM:::.pfimPopMuChainFactors( parameters ),
    omega_iiv = unname( map_dbl( distrib, ~ prop( .x, "omega" )^2 ) ),
    gamma_values = vapply( parameters, function( x ) pluck( x, "gamma", .default = 0 ), numeric( 1L ) ),
    error_variance = if ( length( grads ) == 1L ) {
      as.matrix( varRes[[ iter ]]$variances[[ 1L ]]$variance$errorVariance )
    } else {
      as.matrix( bdiag( map( seq_len( nOcc ), ~ as.matrix(
        varRes[[ iter ]]$variances[[ .x ]]$variance$errorVariance
      ) ) ) )
    },
    occ_col_widths = vapply( grads, ncol, integer( 1L ) ),
    sigma_derivatives = map( seq_len( sigma_n ), function( i ) {
      sdo = map( seq_len( nOcc ), ~ varRes[[ iter ]]$variances[[ .x ]]$variance$sigmaDerivatives[[ i ]] )
      if ( length( sdo ) == 1L ) as.matrix( sdo[[ 1L ]] ) else as.matrix( bdiag( sdo ) )
    } ),
    has_iov = any( vapply( parameters, function( x ) pluck( x, "gamma", .default = 0 ), numeric( 1L ) ) > 0 )
  )
}

.callPopFimComboCpp = function( inp ) {
  do.call(
    computePopFimCombo_Rcpp,
    c(
      inp[ c( "gradients", "mu_values", "omega_iiv", "error_variance",
              "occ_col_widths", "sigma_derivatives" ) ],
      list( gamma = inp$gamma_values, has_iov = inp$has_iov )
    )
  )
}

expect_combo_blocks_match = function( evaluation, case_id, tol = 1e-8 ) {
  arm = prop( prop( evaluation, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]
  purrr::walk( seq_along( prop( arm, "evaluationGradients" ) ), function( iter ) {
    inp = .popFimComboInputs( evaluation, iter = iter, arm = arm )
    expect_equal(
      .callPopFimComboCpp( inp ),
      do.call( .computePopFimCombo_R, inp ),
      tolerance = tol,
      info = paste( case_id, "iter", iter )
    )
  } )
}

purrr::walk( cov_iov_evaluation_cases(), function( case ) {
  if ( case$id == "cas1_noCov_noIOV" ) return()
  local( {
    c = case
    test_that( paste0( "C++ combo matches R: ", c$id ), {
      evaluation = run( build_cov_iov_evaluation( c ) )
      expect_combo_blocks_match( evaluation, c$id )
    } )
  } )
})

test_that( "cas6 population FIM D matches reference (C++ path)", {
  cas = Filter(
    function( x ) x$id == "cas6_noCovFixed_withCovOccasion_withIOV",
    cov_iov_evaluation_cases()
  )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas6_occasion_iov$D, tolerance = 1e-6 )
} )

test_that( "IOV gamma block splits per occasion when n_occasions > 2", {
  n_omega = 3L
  n_per   = 2L
  n_occ   = 3L
  p       = n_per * n_occ
  set.seed( 42L )
  G = matrix( rnorm( n_omega * p ), nrow = n_omega, ncol = p )
  mu    = c( 1, 2, 3 )
  omega = rep( 0.09, n_omega )
  gamma = c( 0, sqrt( 0.0225 ), sqrt( 0.0225 ) )
  R     = diag( p ) * 0.01

  inp = list(
    gradients         = G,
    mu_values         = mu,
    omega_iiv         = omega,
    gamma_values      = gamma,
    error_variance    = R,
    occ_col_widths    = rep( n_per, n_occ ),
    sigma_derivatives = list( diag( p ) * 0.001 ),
    has_iov           = TRUE
  )

  expect_equal(
    .callPopFimComboCpp( inp ),
    do.call( .computePopFimCombo_R, inp ),
    tolerance = 1e-10
  )

  wrong = inp
  wrong$occ_col_widths = c( n_per, n_per * 2L )
  expect_failure(
    expect_equal(
      .callPopFimComboCpp( wrong ),
      do.call( .computePopFimCombo_R, inp ),
      tolerance = 1e-10
    )
  )
} )

test_that( "population covariate path on cas10", {
  evaluation = run( cas10_evaluation( "pop_cov_combo", "population" ) )
  M = prop( PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation ), "fisherMatrix" )
  expect_gt( det( M ), 0 )
  expect_gt( ncol( M ), 6L )
} )

test_that( ".pfimCovOccasionCacheGet increments hits on reuse", {
  local_pfim_opts( list(
    covariate.occasion.cache = TRUE,
    covariate.occasion.cache.hits = 0L,
    fim.cache.scope = NULL
  ) )
  PFIM:::.pfimCovOccasionCacheClear()

  evaluation = cas10_evaluation( "cov_occ_cache", "population" )
  pfim_set_option( fim.cache.scope = PFIM:::.pfimProjectScopeId( evaluation ) )
  PFIM:::.pfimSetCacheEvaluation( evaluation )
  arm = prop( prop( evaluation, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  occ = list( evaluation = list( RespPK = data.frame( time = 1, RespPK = 1 ) ) )
  params = c( ka = 1, V = 3.5, Cl = 2 )
  PFIM:::.pfimCovOccasionCacheSet(
    arm, "combo1", "occ1", params, occ, wantModel = TRUE, wantGrad = FALSE
  )
  PFIM:::.pfimCovOccasionCacheGet(
    arm, "combo1", "occ1", params, wantModel = TRUE, wantGrad = FALSE
  )
  PFIM:::.pfimCovOccasionCacheGet(
    arm, "combo1", "occ1", params, wantModel = TRUE, wantGrad = FALSE
  )
  expect_gte( pfim_get_option( "covariate.occasion.cache.hits" ), 1L )
} )

test_that( ".pfimCovOccasionCacheKey stays within env name limit after ODE grid expansion", {
  modelEquations = list(
    "Deriv_RespPK" = "dose_RespPK/V * ka  * exp( -ka * t ) - Cl/V * RespPK",
    "Deriv_RespPD" = "Rin*(1-Imax*RespPK/(RespPK+C50))-kout*RespPD"
  )
  modelParameters = list(
    ModelParameter( name = "V",  distribution = LogNormal( mu = 8, omega = sqrt( 0.02 ) ), gamma = sqrt( 0.0225 ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 0.13, omega = sqrt( 0.06 ) ), gamma = sqrt( 0.0225 ) ),
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1.6, omega = sqrt( 0.70 ) ), gamma = sqrt( 0.0225 ) ),
    ModelParameter( name = "Rin",  distribution = LogNormal( mu = 5.4, omega = sqrt( 0.2 ) ) ),
    ModelParameter( name = "kout", distribution = LogNormal( mu = 0.06, omega = sqrt( 0.02 ) ) ),
    ModelParameter( name = "Imax", distribution = LogNormal( mu = 1, omega = 0 ), fixedMu = TRUE, fixedOmega = TRUE ),
    ModelParameter( name = "C50",  distribution = LogNormal( mu = 1.2, omega = sqrt( 0.01 ) ) )
  )
  modelError = list(
    Combined1( output = "RespPK", sigmaInter = 0.6, sigmaSlope = 0.07 ),
    Constant( output = "RespPD", sigmaInter = 4 )
  )
  cov = Covariate(
    name = "CYP2C9",
    categories = c( "Wild", "Others" ),
    categoriesProportions = c( 0.6, 0.4 ),
    effects = list( "Others" = c( "Cl" = log( 0.5 ) ) )
  )
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  arm = Arm(
    name = "arm1", size = 32, administrations = list( admin ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 2, 6, 12, 24, 48, 120 ) ),
      SamplingTimes( outcome = "RespPD", samplings = c( 0, 24, 48, 72, 120, 144 ) )
    ),
    initialCondition = list( RespPK = 0, RespPD = 0 )
  )
  design = Design( name = "design1", arms = list( arm ) )
  ev = run( Evaluation(
    name = "ode_cov_iov_cache_key",
    modelEquations = modelEquations, modelParameters = modelParameters,
    modelError = modelError, modelCovariates = list( cov ),
    modelCovariatesEquation = "exponential", designs = list( design ),
    outputs = list( RespPK = "RespPK", RespPD = "RespPD" ),
    numberOfOccasions = 2L, fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  ) )

  PFIM:::.pfimClearGradientPerfCaches()
  PFIM:::.invalidateEvalModelCache( ev )
  model = rebuildEvalModel( ev, finiteDifference = TRUE )
  arm = updateSamplingTimes( arm, getSamplingData( arm ) )
  model = PFIM:::.pfimPrepareModelForEvaluation( PFIM:::.pfimCloneS7( model ), arm )
  combo = prop( model, "covariatesCombination" )$combinations$name[[ 1L ]]
  params = prop( model, "modelParametersWithCovariates" )[[ combo ]][[ "occ1" ]]
  key = PFIM:::.pfimCovOccasionCacheKey( arm, combo, "occ1", params, model )
  expect_lte( nchar( key ), 10000L )

  plotOptions = list( unitTime = "h", unitOutcomes = c( RespPK = "mg/L", RespPD = "" ) )
  expect_error( plotEvaluation( ev, plotOptions ), NA )
  expect_error( plotSensitivityIndices( ev, plotOptions ), NA )
} )
