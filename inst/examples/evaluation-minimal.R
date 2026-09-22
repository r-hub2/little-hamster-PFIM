library(PFIM)
ev = Evaluation(
  name = "Linear1BolusSingleDose_kV",
  modelFromLibrary = list(PKModel = "Linear1BolusSingleDose_kV"),
  modelParameters = list(
    ModelParameter(name = "k", distribution = LogNormal(mu = 0.25, omega = sqrt(0.25))),
    ModelParameter(name = "V", distribution = LogNormal(mu = 15, omega = sqrt(0.1)))
  ),
  modelError = list(Combined1(output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15)),
  outputs = list("RespPK"),
  designs = list(Design(
    name = "design1",
    arms = list(Arm(
      name = "arm1", size = 200,
      administrations = list(Administration(outcome = "RespPK", timeDose = 0, dose = 100)),
      samplingTimes = list(SamplingTimes(outcome = "RespPK", samplings = c(0.33, 1.5, 5, 12)))
    ))
  )),
  fimType = "population"
)
ev = run(ev)
dim(getFisherMatrix(ev)$fisherMatrix)
