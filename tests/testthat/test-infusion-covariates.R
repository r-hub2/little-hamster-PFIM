# Covariate dispatch and FIM on infusion / steady-state / ODE-infusion models.

.prep_eval_model_type = function( ev ) {
  if ( length( projectProp( ev, "modelFromLibrary" ) ) > 0L )
    prop( ev, "modelEquations" ) = defineModelEquationsFromLibraryOfModel( ev )
  defineModelType( ev )
}

.infusion_cov_sex = function() {
  Covariate(
    name                    = "Sex",
    categories              = c( "M", "F" ),
    categoriesProportions   = c( 0.5, 0.5 ),
    effects                 = list( "F" = c( "V" = log( 1.2 ) ) )
  )
}

.infusion_analytic_cov_ev = function( fimType = "population" ) {
  admin = Administration( outcome = "RespPK", Tinf = 2, timeDose = 0, dose = 30 )
  arm   = Arm(
    name              = "a1",
    size              = 40,
    administrations   = list( admin ),
    samplingTimes     = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 8 ) )
    )
  )
  Evaluation(
    name                    = "inf_analytic_cov",
    modelFromLibrary        = list( PKModel = "Linear1InfusionSingleDose_ClV" ),
    modelParameters         = list(
      ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
      ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) )
    ),
    modelCovariates         = list( .infusion_cov_sex() ),
    modelCovariatesEquation = "exponential",
    modelError              = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    designs                 = list( Design( name = "d1", arms = list( arm ) ) ),
    fimType                 = fimType,
    outputs                 = list( "RespPK" )
  )
}

.infusion_steady_cov_ev = function( fimType = "population" ) {
  admin = Administration( outcome = "RespPK", Tinf = 5, tau = 5, dose = 20 )
  arm   = Arm(
    name              = "a1",
    size              = 40,
    administrations   = list( admin ),
    samplingTimes     = list(
      SamplingTimes(
        outcome   = "RespPK",
        samplings = c( 0, 1, 2, 5, 7, 8, 10, 12, 14, 15, 16, 20 )
      )
    )
  )
  Evaluation(
    name                    = "inf_ss_cov",
    modelFromLibrary        = list( PKModel = "Linear1InfusionSteadyState_ClVtau" ),
    modelParameters         = list(
      ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
      ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) )
    ),
    modelCovariates         = list( .infusion_cov_sex() ),
    modelCovariatesEquation = "exponential",
    modelError              = list( Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 ) ),
    designs                 = list( Design( name = "d1", arms = list( arm ) ) ),
    fimType                 = fimType,
    outputs                 = list( "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
}

.bolus_steady_cov_ev = function( fimType = "population" ) {
  admin = Administration( outcome = "RespPK", tau = 12, timeDose = 0, dose = 100 )
  arm   = Arm(
    name              = "a1",
    size              = 40,
    administrations   = list( admin ),
    samplingTimes     = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 1, 4, 8, 12, 16, 20 ) )
    )
  )
  Evaluation(
    name                    = "bolus_ss_cov",
    modelFromLibrary        = list( PKModel = "Linear1BolusSteadyState_ClVtau" ),
    modelParameters         = list(
      ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
      ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) )
    ),
    modelCovariates         = list( .infusion_cov_sex() ),
    modelCovariatesEquation = "exponential",
    modelError              = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    designs                 = list( Design( name = "d1", arms = list( arm ) ) ),
    fimType                 = fimType,
    outputs                 = list( "RespPK" )
  )
}

test_that( "covariate evaluation core dispatches infusion and steady-state models", {
  ev_inf  = .infusion_analytic_cov_ev()
  ev_iss  = .infusion_steady_cov_ev()
  ev_bss  = .bolus_steady_cov_ev()

  m_inf  = rebuildEvalModel( ev_inf )
  m_iss  = rebuildEvalModel( ev_iss )
  m_bss  = rebuildEvalModel( ev_bss )

  expect_identical( PFIM:::.covariateEvaluationCore( m_inf ), PFIM:::evaluateAnalyticInfusionCore )
  expect_identical( PFIM:::.covariateEvaluationCore( m_iss ), PFIM:::evaluateAnalyticInfusionSteadyStateCore )
  expect_identical( PFIM:::.covariateEvaluationCore( m_bss ), PFIM:::evaluateAnalyticSteadyStateCore )
  expect_true( PFIM:::.hasCovariateEvaluationCore( m_inf ) )
  expect_true( PFIM:::.hasCovariateEvaluationCore( m_iss ) )
  expect_true( PFIM:::.hasCovariateEvaluationCore( m_bss ) )
} )

test_that( "infusion analytic population FIM with covariate expands combinations", {
  ev  = run( .infusion_analytic_cov_ev( "population" ) )
  arm = prop( prop( ev, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]
  expect_gt( length( prop( arm, "evaluationGradients" ) ), 1L )
  .expect_valid_fim( ev, min_dim = 3L )
} )

test_that( "infusion steady-state population FIM with covariate runs", {
  ev = run( .infusion_steady_cov_ev( "population" ) )
  expect_s7_class( .prep_eval_model_type( .infusion_steady_cov_ev() ), ModelAnalyticInfusionSteadyState )
  arm = prop( prop( ev, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]
  expect_gt( length( prop( arm, "evaluationGradients" ) ), 1L )
  .expect_valid_fim( ev, min_dim = 3L )
} )

test_that( "bolus steady-state population FIM with covariate uses steady-state core", {
  ev = run( .bolus_steady_cov_ev( "population" ) )
  expect_s7_class( .prep_eval_model_type( .bolus_steady_cov_ev() ), ModelAnalyticSteadyState )
  arm = prop( prop( ev, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]
  expect_gt( length( prop( arm, "evaluationGradients" ) ), 1L )
  .expect_valid_fim( ev, min_dim = 3L )
} )

test_that( "ODE infusion population FIM with covariate runs", {
  pk_model = "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2"
  params   = list(
    ModelParameter( name = "Vm",  distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km",  distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) ),
    ModelParameter( name = "V1",  distribution = LogNormal( mu = 10,   omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "V2",  distribution = LogNormal( mu = 20,   omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "k12", distribution = LogNormal( mu = 0.5,  omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "k21", distribution = LogNormal( mu = 0.3,  omega = sqrt( 0.20 ) ) )
  )
  ev = Evaluation(
    name                    = "ode_inf_cov",
    modelFromLibrary        = list( PKModel = pk_model ),
    modelParameters         = params,
    modelCovariates         = list( Covariate(
      name = "Sex", categories = c( "M", "F" ),
      categoriesProportions = c( 0.5, 0.5 ),
      effects = list( "F" = c( "V1" = log( 1.2 ) ) )
    ) ),
    modelCovariatesEquation = "exponential",
    modelError              = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs                 = list( RespPK = "C1" ),
    designs = list( Design(
      name = "d1",
      arms = list( Arm(
        name              = "a1",
        size              = 40,
        administrations   = list(
          Administration( outcome = "RespPK", Tinf = 2, timeDose = 0, dose = 50 )
        ),
        samplingTimes     = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 6, 12 ) )
        ),
        initialCondition  = list( C1 = 0, C2 = 0 )
      ) )
    ) ),
    fimType                 = "population",
    odeSolverParameters     = .ode_solver_params()
  )
  model = .prep_eval_model_type( ev )
  expect_s7_class( model, ModelODEInfusionDoseInEquation )
  expect_identical(
    PFIM:::.covariateEvaluationCore( rebuildEvalModel( ev ) ),
    PFIM:::.odeInfusionEvaluateModelCore
  )
  ev = run( ev )
  arm = prop( prop( ev, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]
  expect_gt( length( prop( arm, "evaluationGradients" ) ), 1L )
  .expect_valid_fim( ev, min_dim = 7L )
} )
