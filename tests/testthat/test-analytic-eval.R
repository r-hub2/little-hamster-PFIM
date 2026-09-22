# Analytic evaluation helpers (shared by ModelAnalytic*).

test_that( ".analyticSamplingGrid is sorted unique", {
  arm = Arm(
    name = "a", size = 10,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes   = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 2, 0.5, 2 ) )
    )
  )
  expect_equal( PFIM:::.analyticSamplingGrid( arm ), c( 0.5, 2 ) )
} )

test_that( ".analyticOutcomeFrame drops empty passive columns", {
  expect_equal(
    names( PFIM:::.analyticOutcomeFrame( 1.5, NULL ) ),
    "Admin"
  )
  expect_equal(
    names( PFIM:::.analyticOutcomeFrame( 1.5, 0.2 ) ),
    c( "Admin", "NoAdmin" )
  )
} )

test_that( ".analyticFinishEvaluation restricts by outcome sampling", {
  tmp = data.frame(
    time = c( 0.5, 2, 6 ),
    RespPK = c( 1, 2, 3 ),
    RespPD = c( 10, 20, 30 )
  )
  arm = Arm(
    name = "a", size = 10,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 6 ) ),
      SamplingTimes( outcome = "RespPD", samplings = 2 )
    )
  )
  out = PFIM:::.analyticFinishEvaluation( tmp, c( "RespPK", "RespPD" ), arm )
  expect_equal( out$RespPK$time, c( 0.5, 6 ) )
  expect_equal( out$RespPD$time, 2 )
} )

test_that( ".pfimZeroCrossInfDiag isolates Inf coordinates", {
  C = diag( 3 )
  C = PFIM:::.pfimZeroCrossInfDiag( C, 2L )
  expect_equal( C[ 1L, 1L ], 1 )
  expect_true( is.infinite( C[ 2L, 2L ] ) )
  expect_equal( C[ 1L, 2L ], 0 )
  expect_equal( C[ 2L, 3L ], 0 )
} )
