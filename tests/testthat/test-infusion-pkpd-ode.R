# PK infusion analytique + PD ODE (tests_PFIM pk_analytic_infusion_pd_ode).

test_that("Linear1InfusionSingleDose_ClV + Turnover PD ODE, outputs C1/C2", {
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
    Combined1(output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1),
    Constant(output = "RespPD", sigmaInter = 0.8)
  )

  administration = Administration(
    outcome = "C1",
    Tinf = c(2, 5, 10),
    timeDose = c(0, 10, 20),
    dose = c(10, 20, 30)
  )
  arm1 = Arm(
    name = "BrasTest",
    size = 40,
    administrations = list(administration),
    samplingTimes = list(
      SamplingTimes(outcome = "RespPK", samplings = c(0, 1, 2, 5, 10, 12, 20)),
      SamplingTimes(outcome = "RespPD", samplings = c(0, 2, 10, 20, 30))
    ),
    initialConditions = list( C1 = 0, C2 = 0 )
  )
  design1 = Design(name = "design1", arms = list(arm1))

  evaluation = Evaluation(
    name = "",
    modelFromLibrary = modelFromLibrary,
    modelParameters = modelParameters,
    modelError = modelError,
    outputs = list( "RespPK" = "C1", "RespPD" = "C2" ),
    designs = list(design1),
    fimType = "population",
    odeSolverParameters = list(atol = 1e-4, rtol = 1e-4)
  )

  evaluation = run( evaluation )
  .expect_valid_fim( evaluation, min_dim = 6L )
})
