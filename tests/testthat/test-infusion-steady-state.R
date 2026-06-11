# Integration test: library model Linear1InfusionSteadyState_ClVtau
# (mirrors tests_PFIM/Evaluation_tests/library_of_models/pk_analytic_infusion_steady_state)

test_that("Linear1InfusionSteadyState_ClVtau population FIM runs", {
  modelFromLibrary = list(PKModel = "Linear1InfusionSteadyState_ClVtau")
  modelParameters = list(
    ModelParameter(name = "V",  distribution = LogNormal(mu = 3.5, omega = 0.09)),
    ModelParameter(name = "Cl", distribution = LogNormal(mu = 2,   omega = 0.09))
  )
  modelError = list(Combined1(output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1))

  administration = Administration(
    outcome = "RespPK",
    Tinf = 5,
    tau = 5,
    dose = 20
  )
  samplingTimes = SamplingTimes(
    outcome = "RespPK",
    samplings = c(0, 1, 2, 5, 7, 8, 10, 12, 14, 15, 16, 20)
  )
  arm1 = Arm(
    name = "BrasTest",
    size = 40,
    administrations = list(administration),
    samplingTimes = list(samplingTimes)
  )
  design1 = Design(name = "design1", arms = list(arm1))

  evaluation = Evaluation(
    name = "",
    modelFromLibrary = modelFromLibrary,
    modelParameters = modelParameters,
    modelError = modelError,
    outputs = list("RespPK"),
    designs = list(design1),
    fimType = "population",
    odeSolverParameters = list(atol = 1e-8, rtol = 1e-8)
  )

  prop( evaluation, "modelEquations" ) = defineModelEquationsFromLibraryOfModel( evaluation )
  expect_s7_class( defineModelType( evaluation ), ModelAnalyticInfusionSteadyState )
  evaluation = run(evaluation)
  expect_true(length(prop(evaluation, "evaluationDesign")) > 0L)
})
