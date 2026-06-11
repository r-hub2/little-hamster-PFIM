# PK MODELS

test_that("Linear1BolusSingleDose_kV", {

  modelFromLibrary = list("PKModel" = "Linear1BolusSingleDose_kV")

  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = sqrt(0.25) ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 15, omega = sqrt(0.1) ) ) )

  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
  modelError = list( errorModelRespPK )

  # administration
  administration = Administration( outcome = "RespPK", timeDose = c( 0 ), dose = c( 100 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.33, 1.5, 5, 12 ) )

  # arm
  arm1 = Arm( name = "BrasTest", size = 200, administrations  = list( administration ) , samplingTimes    = list( samplingTimes ) )

  # design
  design1 = Design( name = "design1", arms = list( arm1 ) )

  # Evaluation
  evaluationFIM = Evaluation( name = "Linear1BolusSingleDose_kV",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK" ),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )
  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  valueDetPopulationFim = 1456484950098966784
  tol = 1e-6
  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol)
})

############################################################################################################################

test_that("Linear1BolusSingleDose_ClV", {

  modelFromLibrary = list("PKModel" = "Linear1BolusSingleDose_ClV")

  # model modelParameters
  modelParameters = list( ModelParameter( name = "Cl", distribution = LogNormal( mu = 3.75, omega = sqrt(0.25) ) ),
                          ModelParameter( name = "V",  distribution = LogNormal( mu = 15, omega = sqrt(0.1) ) ) )

  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
  modelError = list( errorModelRespPK )

  # administration
  administration = Administration( outcome = "RespPK", timeDose = c( 0 ), dose = c( 100 ) )
  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.33, 1.5, 5, 12 ) )
  # arm
  arm1 = Arm( name = "BrasTest",
              size = 200,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes ) )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear1BolusSingleDose_ClV",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK" ),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )

  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  valueDetPopulationFim = 7292815928272215
  tol = 1e-6
  expect_equal(detPopulationFim, valueDetPopulationFim, tolerance = tol)
})

test_that("Linear1InfusionSingleDose_kV", {

  # --------------------------------------
  # model definition


  modelFromLibrary = list("PKModel" = "Linear1InfusionSingleDose_kV")


  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt(0.09) ) ),
    ModelParameter( name = "k",  distribution = LogNormal( mu = 0.6, omega = sqrt(0.09) ) )
  )

  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 )
  modelError = list( errorModelRespPK )

  # administration
  administration = Administration( outcome = "RespPK",  Tinf = c(2), timeDose = c( 0 ), dose = c( 30 ) )
  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 4, 8 ) )
  # arm
  arm1 = Arm( name = "BrasTest",
              size = 40,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes )
  )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear1InfusionSingleDose_kV",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK" ),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )

  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  detPopulationFim

  valueDetPopulationFim = 1976729482139640576.0
  tol = 1e-6
  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol )

})


test_that("Model PK 1cpt : Linear1InfusionSingleDose_VCl", {

  # --------------------------------------
  # model definition

  # model equations

  modelFromLibrary = list("PKModel" = "Linear1InfusionSingleDose_ClV") #Linear1InfusionSingleDose_VCl


  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "V", distribution = LogNormal( mu = 3.5, omega = sqrt(0.09) ) ),
    ModelParameter( name = "Cl",  distribution = LogNormal( mu = 2, omega = sqrt(0.09) ) )
  )

  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 )
  modelError = list( errorModelRespPK )


  # administration
  administration = Administration( outcome = "RespPK",
                                   Tinf = c(2),
                                   timeDose = c( 0 ),
                                   dose = c( 30 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 4, 8 ) )

  # arm
  arm1 = Arm( name = "BrasTest",
              size = 40,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes )
  )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear1InfusionSingleDose_ClV",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK" ),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )

  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  detPopulationFim


  valueDetPopulationFim = 155857928552697888.0

  tol = 1e-6

  expect_equal( detPopulationFim, valueDetPopulationFim, tolerance = tol )

})


