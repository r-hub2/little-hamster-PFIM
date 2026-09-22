purrr::walk( cov_iov_evaluation_cases(), function( case ) {
  local( {
    c = case
    test_that( paste0( "cov/IOV: ", c$id ), {
      evaluation = run( build_cov_iov_evaluation( c ) )
      fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
      M   = prop( fim, "fisherMatrix" )
      expect_true( is.matrix( M ) )
      expect_equal( ncol( M ), c$ncol )
      expect_gt( det( M ), 0 )

      if ( isTRUE( c$cov_test ) ) {
        ct = covariateTest( evaluation, target_power = 0.90, alpha = 0.05 )
        expect_s3_class( ct, "PFIM::CovariateTest" )
        expect_gt( nrow( prop( ct, "covariate" ) ), 0L )
      }
    } )
  } )
})

test_that( "cas2 IOV without covariates matches population FIM reference", {
  cas = Filter( function( x ) x$id == "cas2_noCov_withIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M   = prop( fim, "fisherMatrix" )

  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas2_noCov_withIOV$D, tolerance = 1e-6 )
  expect_equal( M[ 1L, 1L ], .pfimGold$cas2_noCov_withIOV$M11, tolerance = 1e-4 )
  expect_equal( M[ 8L, 8L ], .pfimGold$cas2_noCov_withIOV$M88, tolerance = 1e-3 )
  expect_equal( M[ 10L, 10L ], .pfimGold$cas2_noCov_withIOV$M1010, tolerance = 1e-4 )
} )

test_that( "cas1 population FIM matches reference", {
  cas = Filter( function( x ) x$id == "cas1_noCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas1_noCov_noIOV$D, tolerance = 1e-6 )
} )

test_that( "cas3 population FIM matches Monolix reference", {
  cas = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M   = prop( fim, "fisherMatrix" )

  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas3_withCov_noIOV$D, tolerance = 1e-6 )
  expect_equal( M[ 1L, 1L ], .pfimGold$cas3_withCov_noIOV$M11, tolerance = 1e-4 )
  expect_equal( M[ 5L, 5L ], .pfimGold$cas3_withCov_noIOV$M55, tolerance = 1e-4 )
  expect_equal( M[ 8L, 8L ], .pfimGold$cas3_withCov_noIOV$M88, tolerance = 1e-4 )
} )

test_that( "cas4 population FIM matches reference", {
  cas = Filter( function( x ) x$id == "cas4_withCov_withIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas4_withCov_withIOV$D, tolerance = 1e-6 )
} )

test_that( "cas5 population FIM matches reference", {
  cas = Filter( function( x ) x$id == "cas5_withCov_mixed", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas5_withCov_mixed$D, tolerance = 1e-6 )
} )

test_that( "cas6 occasion covariate + IOV matches reference", {
  cas = Filter(
    function( x ) x$id == "cas6_noCovFixed_withCovOccasion_withIOV",
    cov_iov_evaluation_cases()
  )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M   = prop( fim, "fisherMatrix" )

  expect_equal( as.numeric( PFIM:::Dcriterion( fim ) ), .pfimGold$cas6_occasion_iov$D, tolerance = 1e-6 )
  expect_equal( M[ 1L, 1L ], .pfimGold$cas6_occasion_iov$M11, tolerance = 1e-4 )
  expect_equal( M[ 10L, 10L ], .pfimGold$cas6_occasion_iov$M1010, tolerance = 1e-3 )
} )
