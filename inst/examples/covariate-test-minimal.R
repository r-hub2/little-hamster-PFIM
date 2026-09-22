library(PFIM)
ev = Evaluation(
  name = "covariate_test",
  modelEquations = list(
    RespPK = "dose_RespPK/V * ka/(ka - Cl/V) * (exp(-Cl/V * t) - exp(-ka * t))"
  ),
  modelParameters = list(
    ModelParameter(name = "ka", distribution = LogNormal(mu = 1, omega = sqrt(0.09))),
    ModelParameter(name = "V",  distribution = LogNormal(mu = 3.5, omega = sqrt(0.09))),
    ModelParameter(name = "Cl", distribution = LogNormal(mu = 2, omega = sqrt(0.09)))
  ),
  modelCovariates = list(
    Covariate(
      name = "Sex", categories = c("M", "F"), categoriesProportions = c(0.5, 0.5),
      effects = list(F = c(V = log(1.2)))
    )
  ),
  modelCovariatesEquation = "exponential",
  modelError = list(Constant(output = "RespPK", sigmaInter = 0.1)),
  outputs = list(RespPK = "RespPK"),
  designs = list(Design(
    name = "design1",
    arms = list(Arm(
      name = "arm1", size = 40,
      administrations = list(Administration(outcome = "RespPK", timeDose = 0, dose = 30)),
      samplingTimes = list(SamplingTimes(outcome = "RespPK", samplings = c(0.5, 2, 4, 6, 8)))
    ))
  )),
  fimType = "population"
)
ev = run(ev)
ct = covariateTest(ev)
