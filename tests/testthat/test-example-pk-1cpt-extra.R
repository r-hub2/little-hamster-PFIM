test_that("Model PK 1cpt : Linear1FirstOrderSingleDose_kaClV", {

  modelFromLibrary = list("PKModel" = "Linear1FirstOrderSingleDose_kaClV")

  modelParameters = list(
    ModelParameter( name = "V", distribution = LogNormal( mu = 63.000, omega = 0 ) ),
    ModelParameter( name = "Cl",  distribution = LogNormal( mu = 0.513, omega = 0 ) ),
    ModelParameter( name = "ka",  distribution = LogNormal( mu = 1.050, omega = sqrt(0.1) ) )
  )

  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0, sigmaSlope = 0.0676 )
  modelError = list( errorModelRespPK )

  administration = Administration( outcome = "RespPK", timeDose = c( 0 ), dose = c( 5500 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.01, 1, 3, 5, 7, 10, 13, 17, 24 ) )

  # arm
  arm1 = Arm( name = "BrasTest",
              size = 25,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes )
  )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear1FirstOrderSingleDose_kaClV",
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


  valueDetPopulationFim = 88354126194397.8

  tol = 1e-6

  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol )

})

###################################################################################################################################


test_that("Model PK 1cpt : Linear1FirstOrderSingleDose_kaClV", {

  # --------------------------------------
  # model definition

  # model equations

  modelFromLibrary = list("PKModel" = "Linear1FirstOrderSingleDose_kaClV")


  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "V", distribution = LogNormal( mu = 3.5, omega =  sqrt( 0.09 ) ) ),
    ModelParameter( name = "Cl",  distribution = LogNormal( mu = 2.0, omega =  sqrt( 0.09 ) ) ),
    ModelParameter( name = "ka",  distribution = LogNormal( mu = 1.0, omega = sqrt(0.09) ) )
  )


  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 )
  modelError = list( errorModelRespPK )


  # administration
  administration = Administration( outcome = "RespPK",
                                   timeDose = c( 0,12,24,36,48 ),
                                   dose = c( 30,30,30,30,30 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 4, 8, 12.5, 13, 16, 20, 24.5, 25, 28, 32, 36.5, 37, 40, 44, 48.5, 49, 52, 56 ) )

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
  evaluationFIM = Evaluation( name = "Linear1FirstOrderSingleDose_kaClV",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK" ),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )

  # get the determinant of the Fisher matrix
  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  detPopulationFim

  valueDetPopulationFim = 2835909452801529708800884.0

  tol = 1e-6

  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol )

})

###################################################################################################################################


test_that("Model PK 1cpt : Linear1FirstOrderSteadyState_kaClVtau", {

  # --------------------------------------
  # model definition

  # model equations

  modelFromLibrary = list("PKModel" = "Linear1FirstOrderSteadyState_kaClVtau")


  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "ka",  distribution = LogNormal( mu = 1.050, omega = sqrt(0.1) ) ),
    ModelParameter( name = "Cl",  distribution = LogNormal( mu = 0.513, omega =  0 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 63.000, omega =  0 ) )
  )


  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0, sigmaSlope = 0.0676 )
  modelError = list( errorModelRespPK )


  # administration
  administration = Administration( outcome = "RespPK",
                                   tau = c(24),
                                   dose = c( 5500 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.01, 1, 3, 5, 7, 10, 13, 17, 24 ) )

  # arm
  arm1 = Arm( name = "BrasTest",
              size = 25,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes )
  )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear1FirstOrderSteadyState_kaClVtau",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK" ),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )

  # get the determinant of the Fisher matrix
  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  detPopulationFim


  valueDetPopulationFim = 119307107145

  tol = 1e-6

  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol )

})


############################################################################################################################
# END CODE
############################################################################################################################
