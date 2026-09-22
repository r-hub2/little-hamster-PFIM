test_that("gradient performance options preserve FIM results (ODE 1cpt)", {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) )
  arm   = Arm(
    name = "arm", size = 40,
    administrations = list( admin ),
    samplingTimes   = list( st )
  )
  design = Design( name = "d", arms = list( arm ) )
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  err = Constant( output = "RespPK", sigmaInter = 0.1 )
  ev = Evaluation(
    name = "perf_test",
    modelFromLibrary = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" ),
    modelParameters  = mp,
    modelError       = list( err ),
    outputs          = list( "RespPK" ),
    designs          = list( design ),
    fimType          = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )

  opts_off = list(
    perf.adminCache    = FALSE,
    perf.fdCache       = FALSE,
    perf.odeTimesCache = FALSE
  )
  opts_on = list(
    perf.adminCache    = TRUE,
    perf.fdCache       = TRUE,
    perf.odeTimesCache = TRUE
  )

  run_baseline = function( opts ) {
    old = setNames( lapply( names( opts ), pfim_get_option ), names( opts ) )
    do.call( pfim_set_option, opts )
    on.exit( do.call( pfim_set_option, old ), add = TRUE )
    out = run( ev )
    list(
      det = getDeterminant( out ),
      rse = getRSE( out )$RSE
    )
  }

  base = run_baseline( opts_off )
  fast = run_baseline( opts_on )

  expect_equal( fast$det, base$det, tolerance = 1e-6 )
  expect_equal( fast$rse, base$rse, tolerance = 1e-6 )
})

test_that( "finiteDifferenceHessian grid includes fixedMu when omega > 0", {
  pfim_reset_session()
  local_pfim_opts( list( perf.fdLinearOnly = FALSE ) )
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ), fixedMu = TRUE ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  ev = Evaluation(
    name = "fd_fixed_mu",
    modelFromLibrary = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" ),
    modelParameters  = mp,
    modelError       = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs          = list( "RespPK" ),
    designs          = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) ) )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )
  model = rebuildEvalModel( ev, finiteDifference = TRUE )
  grid  = prop( model, "parametersForComputingGradient" )

  # Bayesian MAP needs ∂f/∂μ_V even when fixedMu; population later drops the column.
  expect_equal( grid$freeIdx, c( 1L, 2L, 3L ) )
  expect_equal( ncol( grid$shifted ), 10L )

  ref = run( ev )
  expect_gt( getDeterminant( ref ), 0 )
})

test_that( "fixedMu population FIM excludes mu_V from the Fisher matrix", {
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ), fixedMu = TRUE ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  ev = Evaluation(
    name = "fd_fixed_fim",
    modelFromLibrary = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" ),
    modelParameters  = mp,
    modelError       = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs          = list( "RespPK" ),
    designs          = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) ) )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )

  pfim_reset_session()
  on.exit( pfim_reset_session(), add = TRUE )

  out = run( ev )
  M   = prop( prop( out, "fim" ), "fisherMatrix" )
  rn  = rownames( M )

  expect_false( any( grepl( "^mu_V$|^\u03bc_V$", rn ) ) )
  expect_true( any( grepl( "omega^2_V", rn, fixed = TRUE ) | grepl( "\u03c9\u00b2_V", rn, fixed = TRUE ) ) )
  expect_gt( getDeterminant( out ), 0 )
  expect_true( all( is.finite( M ) ) )
})

test_that( "fixedMu individual FIM with covariates excludes mu_V (aligns with population)", {
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ), fixedMu = TRUE ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  sex = Covariate(
    name = "Sex", categories = c( "M", "F" ), categoriesProportions = c( 0.5, 0.5 ),
    effects = list( "F" = c( "Cl" = log( 1.2 ) ) )
  )
  ev = Evaluation(
    name = "ind_fixed_mu_cov",
    modelFromLibrary        = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" ),
    modelParameters         = mp,
    modelCovariates         = list( sex ),
    modelCovariatesEquation = "exponential",
    modelError              = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs                 = list( "RespPK" ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) ) )
      ) )
    ) ),
    fimType = "individual",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )

  pfim_reset_session()
  on.exit( pfim_reset_session(), add = TRUE )

  out = run( ev )
  fim = PFIM:::setEvaluationFim( prop( out, "fim" ), out )
  M   = prop( fim, "fisherMatrix" )
  rn  = rownames( M )

  expect_false( any( grepl( "mu_V|μ_V|_V$", rn ) & grepl( "mu|μ", rn ) ) )
  expect_false( any( grepl( "^mu_V$|^μ_V$", rn ) ) )
  expect_true( any( grepl( "mu_ka|\u03bc_ka|<U\\+03BC>_ka", rn ) ) )
  expect_true( any( grepl( "mu_Cl|\u03bc_Cl|<U\\+03BC>_Cl", rn ) ) )
  # Individual FIM omits covariate beta columns.
  expect_false( any( grepl( "beta_|β_", rn ) ) )
  expect_true( all( is.finite( M ) ) )
  expect_true( is.finite( det( M ) ) && abs( det( M ) ) > 0 )
})

test_that( "FD grid includes fixedMu when gamma > 0 and omega = 0", {
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter(
      name = "V", distribution = LogNormal( mu = 3.5, omega = 0 ),
      fixedMu = TRUE, gamma = 0.2
    ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  expect_equal( PFIM:::.pfimFdFreeIndices( mp ), c( 1L, 2L, 3L ) )
})

test_that( ".paramOmegaEstimable ignores fixedMu", {
  p_free = ModelParameter(
    name = "V", distribution = LogNormal( mu = 3.5, omega = 0.3 )
  )
  p_fix_mu = ModelParameter(
    name = "V", distribution = LogNormal( mu = 3.5, omega = 0.3 ),
    fixedMu = TRUE
  )
  p_fix_om = ModelParameter(
    name = "V", distribution = LogNormal( mu = 3.5, omega = 0.3 ),
    fixedOmega = TRUE
  )
  p_zero = ModelParameter(
    name = "V", distribution = LogNormal( mu = 3.5, omega = 0 ),
    fixedMu = TRUE
  )
  expect_true( PFIM:::.paramOmegaEstimable( p_free ) )
  expect_true( PFIM:::.paramOmegaEstimable( p_fix_mu ) )
  expect_false( PFIM:::.paramOmegaEstimable( p_fix_om ) )
  expect_false( PFIM:::.paramOmegaEstimable( p_zero ) )
})

test_that( "FD grid defaults to linear central stencil (1 + 2*n columns)", {
  grid = PFIM:::.pfimFiniteDifferenceGrid( c( 1, 3.5, 2 ), seq_len( 3L ) )
  expect_equal( ncol( grid$shifted ), 7L )
} )

test_that( "FD quadratic stencil has 1 + 2*n + n(n-1)/2 columns", {
  local_pfim_opts( list( perf.fdLinearOnly = FALSE ) )
  grid = PFIM:::.pfimFiniteDifferenceGrid( c( 1, 3.5, 2 ), seq_len( 3L ) )
  expect_equal( ncol( grid$shifted ), 10L )
} )

test_that( "scaleToOdeTol does not inflate classic FD steps when already >> tol", {
  relStep = .Machine$double.eps^( 1 / 3 )
  pars = c( 1, 3.5, 2 )
  base = PFIM:::.pfimFiniteDifferenceGrid( pars, seq_len( 3L ) )
  ode  = PFIM:::.pfimFiniteDifferenceGrid(
    pars, seq_len( 3L ),
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ),
    scaleToOdeTol = TRUE
  )
  # Classic relative step must be kept (not replaced by odeTol^(1/3) ≈ 2e-3).
  expect_equal( ode$frac[ 2L ], pars[ 1L ] * relStep, tolerance = 1e-12 )
  expect_equal( ode$shifted, base$shifted )
} )

test_that( "FD XcolsInv depends only on nFree, not mu", {
  PFIM:::.pfimEnvClear( PFIM:::.pfimFdSchemeCache )
  on.exit( PFIM:::.pfimEnvClear( PFIM:::.pfimFdSchemeCache ), add = TRUE )
  g1 = PFIM:::.pfimFiniteDifferenceGrid( c( 1, 2, 3 ), seq_len( 3L ) )
  g2 = PFIM:::.pfimFiniteDifferenceGrid( c( 10, 0.5, 7 ), seq_len( 3L ) )
  expect_equal( g1$XcolsInv, g2$XcolsInv )
  expect_false( isTRUE( all.equal( g1$frac, g2$frac ) ) )
  expect_equal( length( PFIM:::.pfimCacheKeys( PFIM:::.pfimFdSchemeCache ) ), 1L )
})

test_that( "tiny |mu| warns once per mu signature about the FD step floor", {
  pfim_reset_session()
  on.exit( pfim_reset_session(), add = TRUE )
  expect_warning(
    PFIM:::.pfimFiniteDifferenceGrid( c( 1e-6, 1 ), seq_len( 2L ) ),
    "absolute step floor"
  )
  expect_no_warning(
    PFIM:::.pfimFiniteDifferenceGrid( c( 1e-6, 1 ), seq_len( 2L ) )
  )
  # A second badly scaled design (distinct tiny mus) must still warn.
  expect_warning(
    PFIM:::.pfimFiniteDifferenceGrid( c( 1e-8, 1 ), seq_len( 2L ) ),
    "absolute step floor"
  )
})

test_that( "scaleToOdeTol raises tiny absolute steps with a warning", {
  pfim_reset_session()
  on.exit( pfim_reset_session(), add = TRUE )
  grid = NULL
  expect_warning(
    {
      grid <- PFIM:::.pfimFiniteDifferenceGrid(
        c( 1e-12, 1e-12 ), seq_len( 2L ),
        odeSolverParameters = list( atol = 1e-6, rtol = 1e-6 ),
        scaleToOdeTol = TRUE
      )
    },
    "absolute step floor"
  )
  # Floor is 10 * odeTol = 1e-5 (pmax(|θ|,1e-4)*eps^(1/3) is ~6e-10 without floor).
  expect_true( all( grid$frac[ 2:3 ] >= 1e-5 - 1e-15 ) )

  expect_no_warning(
    PFIM:::.pfimFiniteDifferenceGrid(
      c( 1e-12, 1e-12 ), seq_len( 2L ),
      odeSolverParameters = list( atol = 1e-6, rtol = 1e-6 ),
      scaleToOdeTol = TRUE
    )
  )

  local_pfim_opts( list( verbose = TRUE ) )
  expect_warning(
    PFIM:::.pfimFiniteDifferenceGrid(
      c( 1e-12, 1e-12 ), seq_len( 2L ),
      odeSolverParameters = list( atol = 1e-6, rtol = 1e-6 ),
      scaleToOdeTol = TRUE
    ),
    "Finite-difference step raised"
  )
})

test_that( "perf.fdLinearOnly reduces FD grid to 1 + 2*n columns", {
  local_pfim_opts( list( perf.fdLinearOnly = TRUE ) )
  grid = PFIM:::.pfimFiniteDifferenceGrid( c( 1, 3.5, 2 ), seq_len( 3L ) )
  expect_equal( ncol( grid$shifted ), 7L )
} )

test_that( "linear FD stencil is close to quadratic on 1cpt ODE", {
  local_pfim_opts( list( perf.fdLinearOnly = FALSE, perf.fdCache = FALSE ) )
  ev = Evaluation(
    name = "fd_linear_cmp",
    modelFromLibrary = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" ),
    modelParameters  = list(
      ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
      ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ) ),
      ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
    ),
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs    = list( "RespPK" ),
    designs    = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) ) )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list()
  )

  run_fd = function( linearOnly ) {
    local_pfim_opts( list( perf.fdLinearOnly = linearOnly, perf.fdCache = FALSE ) )
    run( ev )
  }

  ev_quad = run_fd( FALSE )
  ev_lin  = run_fd( TRUE )
  expect_gt( getDeterminant( ev_quad ), 0 )
  expect_gt( getDeterminant( ev_lin ), 0 )
  expect_equal( getDeterminant( ev_lin ), getDeterminant( ev_quad ), tolerance = 1e-4 )
  rse_lin  = getRSE( ev_lin )
  rse_quad = getRSE( ev_quad )
  expect_equal( as.numeric( rse_lin$RSE ), as.numeric( rse_quad$RSE ), tolerance = 1e-3 )
} )

test_that( "FD hessian session cache key includes ODE tolerances", {
  local_pfim_opts( list( perf.fdCache = TRUE ) )
  PFIM:::.pfimEnvClear( PFIM:::.pfimFdSchemeCache )
  on.exit( PFIM:::.pfimEnvClear( PFIM:::.pfimFdSchemeCache ), add = TRUE )

  params = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 1e-12, omega = 0.3 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 1e-12, omega = 0.3 ) )
  )
  mk = function( atol ) {
    ModelODE(
      name = "fd_cache_tol",
      modelParameters = params,
      odeSolverParameters = list( atol = atol, rtol = atol )
    )
  }
  g_loose = prop(
    suppressWarnings( PFIM:::.pfimFiniteDifferenceHessianCached( mk( 1e-4 ) ) ),
    "parametersForComputingGradient"
  )
  g_tight = prop(
    suppressWarnings( PFIM:::.pfimFiniteDifferenceHessianCached( mk( 1e-8 ) ) ),
    "parametersForComputingGradient"
  )
  expect_false( isTRUE( all.equal( g_loose$frac, g_tight$frac ) ) )
} )
