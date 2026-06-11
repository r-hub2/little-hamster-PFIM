# PK+PD ODE with categorical covariate and unequal sampling times per output.

test_that("PKPD population FIM with covariate runs when outputs have different n samplings", {
  modelEquationsPKPD = list(
    "Deriv_RespPK" = "dose_RespPK/V * ka  * exp( -ka * t ) - Cl/V * RespPK",
    "Deriv_RespPD" = "Rin-kout*(1-Imax*RespPK/(RespPK+C50))*RespPD"
  )
  modelParametersPKPD = list(
    ModelParameter( name = "V",    distribution = LogNormal( mu = 8, omega = sqrt( 0.02 ) ) ),
    ModelParameter( name = "Cl",   distribution = LogNormal( mu = 0.13, omega = sqrt( 0.06 ) ) ),
    ModelParameter( name = "ka",   distribution = LogNormal( mu = 1.6, omega = sqrt( 0.70 ) ) ),
    ModelParameter( name = "Rin",  distribution = LogNormal( mu = 5.4, omega = sqrt( 0.2 ) ) ),
    ModelParameter( name = "kout", distribution = LogNormal( mu = 0.06, omega = sqrt( 0.02 ) ) ),
    ModelParameter(
      name = "Imax", distribution = LogNormal( mu = 1, omega = 0 ),
      fixedMu = TRUE, fixedOmega = TRUE
    ),
    ModelParameter( name = "C50", distribution = LogNormal( mu = 1.2, omega = sqrt( 0.01 ) ) )
  )
  modelErrorPKPD = list(
    Combined1( output = "RespPK", sigmaInter = 0.6, sigmaSlope = 0.07 ),
    Constant( output = "RespPD", sigmaInter = 4 )
  )
  administrationRespPK = Administration(
    outcome = "RespPK", timeDose = 0, dose = 100
  )
  samplingTimesRespPK = SamplingTimes(
    outcome = "RespPK",
    samplings = c( 0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120 )
  )
  samplingTimesRespPD = SamplingTimes(
    outcome = "RespPD",
    samplings = c( 0, 24, 36, 48, 72, 96, 120, 144 )
  )
  empiricalArmPKPD = Arm(
    name = "empiricalArm",
    size = 32,
    administrations = list( administrationRespPK ),
    samplingTimes = list( samplingTimesRespPK, samplingTimesRespPD ),
    initialCondition = list( RespPK = 0, RespPD = 0 )
  )
  empiricalDesignPKPD = Design(
    name = "empiricalDesign", arms = list( empiricalArmPKPD )
  )
  CYP2C9 = Covariate(
    name = "CYP2C9",
    categories = c( "Wild", "Others" ),
    categoriesProportions = c( 0.6, 0.4 ),
    effects = list( "Others" = c( "Cl" = log( 0.5 ) ) )
  )

  evaluation = Evaluation(
    name = "",
    modelEquations = modelEquationsPKPD,
    modelParameters = modelParametersPKPD,
    modelError = modelErrorPKPD,
    modelCovariates = list( CYP2C9 ),
    modelCovariatesEquation = "exponential",
    designs = list( empiricalDesignPKPD ),
    outputs = list( RespPK = "RespPK", RespPD = "RespPD" ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-12, rtol = 1e-12 )
  )

  evaluation = run( evaluation )
  fim = prop( evaluation, "fim" )
  M   = prop( fim, "fisherMatrix" )
  expect_true( is.matrix( M ) )
  expect_true( all( is.finite( M ) ) )
  expect_gt( ncol( M ), 6L )

  se_tbl = prop( fim, "SEAndRSE" )[[ "SEAndRSE" ]]
  expect_true( is.data.frame( se_tbl ) )
  expect_gt( nrow( se_tbl ), 6L )

  show_lines = utils::capture.output( show( evaluation ) )
  expect_true( any( grepl( "RSE\\(%\\)", show_lines, fixed = FALSE ) ) )
  expect_true( any( grepl( "mu_Cl|μ_Cl", show_lines, fixed = FALSE ) ) )
})
