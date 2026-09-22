# Combined2 residual error across population / individual / Bayesian FIM.

.combined2_mini_evaluation = function( fimType ) {
  Evaluation(
    name = paste0( "c2_", fimType ),
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters = list(
      ModelParameter(
        name = "k",
        distribution = LogNormal( mu = 0.25, omega = sqrt( 0.25 ) )
      ),
      ModelParameter(
        name = "V",
        distribution = LogNormal( mu = 15, omega = sqrt( 0.1 ) )
      )
    ),
    modelError = list(
      Combined2( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
    ),
    outputs = list( "RespPK" ),
    designs = list( Design(
      name = "design1",
      arms = list( Arm(
        name = "arm1", size = 50,
        administrations = list(
          Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 2, 4, 8, 12 ) )
        )
      ) )
    ) ),
    fimType = fimType
  )
}

test_that( "Combined2 works for population, individual, and Bayesian FIM", {
  walk( c( "population", "individual", "Bayesian" ), function( fimType ) {
    ev = run( .combined2_mini_evaluation( fimType ) )
    fim = getFisherMatrix( ev )
    expect_true( is.matrix( fim$fisherMatrix ) )
    expect_true( all( is.finite( fim$fisherMatrix ) ) )
    rse = getRSE( ev )
    expect_true( all( is.finite( rse$RSE ) ) )
    expect_true( any( grepl( "slope|inter", rownames( rse ), ignore.case = TRUE ) ) ||
                   identical( fimType, "Bayesian" ) )
  } )
} )

test_that( "plotEvaluation densifies the time grid for display", {
  ev = run( Evaluation(
    name = "plot_dense",
    modelEquations = list( Deriv_Cc = "-k*Cc" ),
    modelParameters = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.2, omega = 0.1 ) )
    ),
    modelError = list(
      Combined2( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.05 )
    ),
    outputs = list( RespPK = "Cc" ),
    designs = list( Design(
      name = "d1",
      arms = list( Arm(
        name = "a1", size = 10,
        administrations = list(
          Administration( outcome = "Cc", timeDose = 0, dose = 10 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 2, 4, 8 ) )
        ),
        initialConditions = list( Cc = "dose_Cc" )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  ) )

  t0 = proc.time()[[ "elapsed" ]]
  plots = plotEvaluation( ev, list() )
  elapsed = proc.time()[[ "elapsed" ]] - t0
  expect_true( !is.null( plots$d1$a1$RespPK ) )
  expect_gt( nrow( plots$d1$a1$RespPK$data ), 5L )
  expect_true( all( c( 0.5, 1, 2, 4, 8 ) %in% plots$d1$a1$RespPK$data$time ) )
  expect_lt( elapsed, 5 )

  si = plotSensitivityIndices( ev, list() )
  p_si = si$d1$a1$RespPK[[ 1L ]]
  expect_s3_class( p_si, "ggplot" )
  expect_gt( nrow( p_si$data ), 5L )
} )
