test_that("Model PK 2cpts : Linear2BolusSingleDose_ClQV1V2", {

  # --------------------------------------
  # model definition

  # model equations

  modelFromLibrary = list("PKModel" = "Linear2BolusSingleDose_ClQV1V2")


  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "Cl",  distribution = LogNormal( mu = 0.4, omega = sqrt(0.2) ) ),
    ModelParameter( name = "V1", distribution = LogNormal( mu = 10, omega = sqrt(0.1) ) ),
    ModelParameter( name = "Q", distribution = LogNormal( mu = 2, omega = sqrt(0.05) ) ),
    ModelParameter( name = "V2", distribution = LogNormal( mu = 50, omega = sqrt(0.4) ) )

  )


  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.6, sigmaSlope = 0.07 )
  modelError = list( errorModelRespPK )


  # administration
  administration = Administration( outcome = "RespPK",
                                   timeDose = c( 0 ),
                                   dose = c( 100 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120) )

  # arm
  arm1 = Arm( name = "BrasTest",
              size = 32,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes )
  )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear2BolusSingleDose_ClQV1V2",
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
  valueDetPopulationFim =  3587145686924841472
  tol = 1e-6
  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol )

})

###################################################################################################################################


test_that("Model PK 2cpts : Linear2BolusSingleDose_kk12k21V", {

  # --------------------------------------
  # model definition

  # model equations

  modelFromLibrary = list("PKModel" = "Linear2BolusSingleDose_kk12k21V")


  # model modelParameters
  modelParameters = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = sqrt(0.25) ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 15.00, omega = sqrt(0.10) ) ),
    ModelParameter( name = "k12",  distribution = LogNormal( mu = 1.00, omega = sqrt(0.40) ) ),
    ModelParameter( name = "k21",  distribution = LogNormal( mu = 0.80, omega = sqrt(0.30) ) )
  )


  # Error Model
  errorModelRespPK = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
  modelError = list( errorModelRespPK )


  # administration
  administration = Administration( outcome = "RespPK",
                                   timeDose = c( 0 ),
                                   dose = c( 100 ) )

  # sampling times
  samplingTimes = SamplingTimes( outcome = "RespPK",
                                 samplings = c(0.33, 1.5, 3, 5, 8, 12 ) )

  # arm
  arm1 = Arm( name = "BrasTest",
              size = 200,
              administrations  = list( administration ) ,
              samplingTimes    = list( samplingTimes )
  )

  # design
  design1 = Design( name = "design1",
                    arms = list( arm1 ) )

  # --------------------------------------
  # Evaluation

  # Evaluate the Fisher Information Matrix for the PopulationFIM
  evaluationFIM = Evaluation( name = "Linear2BolusSingleDose_kk12k21V",
                              modelFromLibrary = modelFromLibrary,
                              modelParameters = modelParameters,
                              modelError = modelError,
                              outputs = list( "RespPK"),
                              designs = list( design1 ),
                              fimType = "population",
                              odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 ) )

  evaluationFIM = run( evaluationFIM )



  FisherMatrix = getFisherMatrix(evaluationFIM )
  detPopulationFim =det(  FisherMatrix$fisherMatrix )
  detPopulationFim

  valueDetPopulationFim =  2137812605034096941482868

  tol = 1e-6

  expect_equal(detPopulationFim,valueDetPopulationFim, tolerance = tol )


})

###################################################################################################################################


