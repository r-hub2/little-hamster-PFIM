test_that( ".modelOmegaIIVVariance reads omegaWithIOV cache", {
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.4 ) )
  )
  model = ModelAnalytic( modelParameters = mp )
  model = evaluateOmegaMatrixFromCovariates( model )
  expect_equal(
    PFIM:::.modelOmegaIIVVariance( model ),
    c( ka = 0.3^2, V = 0.4^2 )
  )
} )

test_that( "infusion ODE admin cache matches full build after mu shift", {
  pfim_reset_session()
  ev = Evaluation(
    name = "infusion_admin_cache",
    modelFromLibrary = list( PKModel = "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2" ),
    modelParameters = list(
      ModelParameter( name = "Vm",  distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
      ModelParameter( name = "Km",  distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) ),
      ModelParameter( name = "V1",  distribution = LogNormal( mu = 10,   omega = sqrt( 0.20 ) ) ),
      ModelParameter( name = "V2",  distribution = LogNormal( mu = 20,   omega = sqrt( 0.20 ) ) ),
      ModelParameter( name = "k12", distribution = LogNormal( mu = 0.5,  omega = sqrt( 0.20 ) ) ),
      ModelParameter( name = "k21", distribution = LogNormal( mu = 0.3,  omega = sqrt( 0.20 ) ) )
    ),
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs = list( RespPK = "C1" ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list(
          Administration( outcome = "RespPK", Tinf = 2, timeDose = 0, dose = 50 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 6, 12 ) )
        ),
        initialConditions = list( C1 = 0, C2 = 0 )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )
  arm = prop( prop( ev, "designs" )[[ 1L ]], "arms" )[[ 1L ]]

  local_pfim_opts( list( perf.adminCache = TRUE ) )
  pfim_reset_session()
  PFIM:::.pfimClearGradientPerfCaches()
  model = rebuildEvalModel( ev, finiteDifference = TRUE )
  invisible( PFIM:::.pfimPrepareModelForEvaluation( PFIM:::.pfimCloneS7( model ), arm ) )

  params = prop( model, "modelParameters" )
  prop( prop( params[[ 1L ]], "distribution" ), "mu" ) = 0.09
  shifted = PFIM:::.pfimCloneS7( model )

  cached = PFIM:::.pfimPrepareModelForEvaluation( shifted, arm )
  full   = defineModelAdministration( PFIM:::.pfimCloneS7( shifted ), arm )

  expect_equal( evaluateModel( cached, arm ), evaluateModel( full, arm ), tolerance = 1e-6 )
} )

test_that( "infusion ODE admin cache preserves population FIM with covariates", {
  ev = Evaluation(
    name = "infusion_ode_admin_fim",
    modelFromLibrary = list( PKModel = "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2" ),
    modelParameters = list(
      ModelParameter( name = "Vm",  distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
      ModelParameter( name = "Km",  distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) ),
      ModelParameter( name = "V1",  distribution = LogNormal( mu = 10,   omega = sqrt( 0.20 ) ) ),
      ModelParameter( name = "V2",  distribution = LogNormal( mu = 20,   omega = sqrt( 0.20 ) ) ),
      ModelParameter( name = "k12", distribution = LogNormal( mu = 0.5,  omega = sqrt( 0.20 ) ) ),
      ModelParameter( name = "k21", distribution = LogNormal( mu = 0.3,  omega = sqrt( 0.20 ) ) )
    ),
    modelCovariates = list(
      Covariate(
        name = "Sex", categories = c( "M", "F" ),
        categoriesProportions = c( 0.5, 0.5 ),
        effects = list( "F" = c( "V1" = log( 1.2 ) ) )
      )
    ),
    modelCovariatesEquation = "exponential",
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs = list( RespPK = "C1" ),
    designs = list( Design(
      name = "d1",
      arms = list( Arm(
        name = "a1", size = 40,
        administrations = list(
          Administration( outcome = "RespPK", Tinf = 2, timeDose = 0, dose = 50 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 6, 12 ) )
        ),
        initialConditions = list( C1 = 0, C2 = 0 )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )

  run_once = function( admin_cache ) {
    pfim_reset_session()
    local_pfim_opts( list(
      perf.adminCache = admin_cache,
      perf.fdCache = FALSE,
      perf.odeTimesCache = FALSE
    ) )
    getDeterminant( run( ev ) )
  }

  expect_equal( run_once( TRUE ), run_once( FALSE ), tolerance = 1e-6 )
} )
