# PK analytique + PD ODE : mapping outputs (RespPK/RespPD -> compartiments).

test_that("Linear1BolusSingleDose_ClV + PD ODE, outputs Cc/E", {
  modelFromLibrary = list(
    PKModel = "Linear1BolusSingleDose_ClV",
    PDModel = "TurnoverkoutEmax_RinEmaxCC50koutE"
  )
  modelParameters = list(
    ModelParameter( name = "V",    distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
    ModelParameter( name = "Cl",   distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) ),
    ModelParameter( name = "Emax", distribution = LogNormal( mu = 10,  omega = sqrt( 0.5 ) ) ),
    ModelParameter( name = "Rin",  distribution = LogNormal( mu = 5,   omega = sqrt( 0.1 ) ) ),
    ModelParameter( name = "kout", distribution = LogNormal( mu = 2,   omega = sqrt( 0.5 ) ) ),
    ModelParameter( name = "C50",  distribution = LogNormal( mu = 2,   omega = sqrt( 0.2 ) ) )
  )
  modelError = list(
    Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 ),
    Constant( output = "RespPD", sigmaInter = 0.8 )
  )
  administration = Administration(
    outcome = "Cc", timeDose = c( 0, 10 ), dose = c( 10, 20 )
  )
  arm1 = Arm(
    name = "BrasTest",
    size = 40,
    administrations = list( administration ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0, 1, 2, 5, 10, 12, 20 ) ),
      SamplingTimes( outcome = "RespPD", samplings = c( 0, 2, 10, 20, 30 ) )
    ),
    initialConditions = list( Cc = 0, E = 0 )
  )
  design1 = Design( name = "design1", arms = list( arm1 ) )

  evaluation = Evaluation(
    name = "",
    modelFromLibrary = modelFromLibrary,
    modelParameters = modelParameters,
    modelError = modelError,
    outputs = list( "RespPK" = "Cc", "RespPD" = "E" ),
    designs = list( design1 ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-4, rtol = 1e-4 )
  )

  evaluation = run( evaluation )
  .expect_valid_fim( evaluation, min_dim = 6L )
})

test_that("outputs Cc/E prime sur les CI (C1/C2)", {
  modelFromLibrary = list(
    PKModel = "Linear1BolusSingleDose_ClV",
    PDModel = "TurnoverkoutEmax_RinEmaxCC50koutE"
  )
  modelParameters = list(
    ModelParameter( name = "V",    distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
    ModelParameter( name = "Cl",   distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) ),
    ModelParameter( name = "Emax", distribution = LogNormal( mu = 10,  omega = sqrt( 0.5 ) ) ),
    ModelParameter( name = "Rin",  distribution = LogNormal( mu = 5,   omega = sqrt( 0.1 ) ) ),
    ModelParameter( name = "kout", distribution = LogNormal( mu = 2,   omega = sqrt( 0.5 ) ) ),
    ModelParameter( name = "C50",  distribution = LogNormal( mu = 2,   omega = sqrt( 0.2 ) ) )
  )
  modelError = list(
    Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 ),
    Constant( output = "RespPD", sigmaInter = 0.8 )
  )
  administration = Administration(
    outcome = "Cc", timeDose = c( 0, 10 ), dose = c( 10, 20 )
  )
  arm1 = Arm(
    name = "BrasTest",
    size = 40,
    administrations = list( administration ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0, 1, 2, 5, 10, 12, 20 ) ),
      SamplingTimes( outcome = "RespPD", samplings = c( 0, 2, 10, 20, 30 ) )
    ),
    initialConditions = list( C1 = 0, C2 = 0 )
  )
  design1 = Design( name = "design1", arms = list( arm1 ) )

  evaluation = Evaluation(
    name = "",
    modelFromLibrary = modelFromLibrary,
    modelParameters = modelParameters,
    modelError = modelError,
    outputs = list( "RespPK" = "Cc", "RespPD" = "E" ),
    designs = list( design1 ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-4, rtol = 1e-4 )
  )

  eq = defineModelEquationsFromLibraryOfModel( evaluation )
  derivNames = names( unlist( eq, use.names = TRUE ) )
  expect_true( all( c( "Deriv_Cc", "Deriv_E" ) %in% derivNames ) )
  expect_false( any( c( "Deriv_C1", "Deriv_C2" ) %in% derivNames ) )
})

test_that("Linear1InfusionSingleDose_ClV + PD ODE, outputs Cc/E", {
  modelFromLibrary = list(
    PKModel = "Linear1InfusionSingleDose_ClV",
    PDModel = "TurnoverkoutEmax_RinEmaxCC50koutE"
  )
  modelParameters = list(
    ModelParameter( name = "V",    distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
    ModelParameter( name = "Cl",   distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) ),
    ModelParameter( name = "Emax", distribution = LogNormal( mu = 10,  omega = sqrt( 0.5 ) ) ),
    ModelParameter( name = "Rin",  distribution = LogNormal( mu = 5,   omega = sqrt( 0.1 ) ) ),
    ModelParameter( name = "kout", distribution = LogNormal( mu = 2,   omega = sqrt( 0.5 ) ) ),
    ModelParameter( name = "C50",  distribution = LogNormal( mu = 2,   omega = sqrt( 0.2 ) ) )
  )
  modelError = list(
    Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 ),
    Constant( output = "RespPD", sigmaInter = 0.8 )
  )
  administration = Administration(
    outcome = "Cc",
    Tinf = c( 2, 5, 10 ),
    timeDose = c( 0, 10, 20 ),
    dose = c( 10, 20, 30 )
  )
  arm1 = Arm(
    name = "BrasTest",
    size = 40,
    administrations = list( administration ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0, 1, 2, 5, 10, 12, 20 ) ),
      SamplingTimes( outcome = "RespPD", samplings = c( 0, 2, 10, 20, 30 ) )
    ),
    initialConditions = list( Cc = 0, E = 0 )
  )
  design1 = Design( name = "design1", arms = list( arm1 ) )

  evaluation = Evaluation(
    name = "",
    modelFromLibrary = modelFromLibrary,
    modelParameters = modelParameters,
    modelError = modelError,
    outputs = list( "RespPK" = "Cc", "RespPD" = "E" ),
    designs = list( design1 ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-4, rtol = 1e-4 )
  )

  evaluation = run( evaluation )
  .expect_valid_fim( evaluation, min_dim = 6L )
})
