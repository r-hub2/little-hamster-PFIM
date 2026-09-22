#' Mock Mult thinOut with joint multi-outcome cells.
#' @keywords internal
.mock_joint_mult_thinOut = function( n_weights = 2L, n_arms_per_cell = 2L ) {
  cell = function( tag ) {
    entries = map( seq_len( n_arms_per_cell ), function( i ) {
      list( arm = .minimal_discrete_arm( paste0( "arm", i, "_", tag ) ) )
    } )
    list( entries = entries )
  }
  w = c( 0.3, 0.7 )[ seq_len( n_weights ) ]
  list(
    listArms       = map( seq_len( n_weights ), cell ),
    weightsIndex   = seq_len( n_weights ),
    optimalWeights = w
  )
}

test_that( ".constraintCellJoint detects multi-outcome cells", {
  thin = .mock_joint_mult_thinOut()
  expect_true( PFIM:::.constraintCellJoint( thin$listArms ) )
  single = list( list( arm = .minimal_discrete_arm( "arm1" ) ) )
  expect_false( PFIM:::.constraintCellJoint( single ) )
} )

test_that( ".pfimShrinkJointMultMixture: joint length coincidence shrinks to one protocol", {
  thin = .mock_joint_mult_thinOut( 2L, 2L )
  expect_equal( length( thin$optimalWeights ), 2L )
  expect_equal( length( thin$listArms ), 2L )
  out = PFIM:::.pfimShrinkJointMultMixture( thin )
  expect_equal( length( out$optimalWeights ), 1L )
  expect_equal( out$optimalWeights, 1 )
  expect_equal( length( out$listArms ), 1L )
  expect_equal( out$weightsIndex, 2L )
} )

test_that( ".pfimShrinkJointMultMixture: no-op when not joint", {
  thin = list(
    listArms       = list( list( arm = .minimal_discrete_arm( "Arm1" ) ) ),
    weightsIndex   = 1L,
    optimalWeights = 1
  )
  out = PFIM:::.pfimShrinkJointMultMixture( thin )
  expect_identical( out, thin )
} )

test_that( ".pfimShrinkJointMultMixture: empty grid and zero-arm cell are no-ops", {
  empty = list(
    listArms       = list(),
    weightsIndex   = integer(),
    optimalWeights = numeric()
  )
  expect_false( PFIM:::.constraintCellJoint( list() ) )
  expect_identical( PFIM:::.pfimShrinkJointMultMixture( empty ), empty )

  zero_cell = list(
    listArms       = list( list( entries = list() ) ),
    weightsIndex   = 1L,
    optimalWeights = 1
  )
  expect_false( PFIM:::.constraintCellJoint( zero_cell$listArms ) )
  expect_identical( PFIM:::.pfimShrinkJointMultMixture( zero_cell ), zero_cell )

  null_cell = list(
    listArms       = list( list( arm = NULL ) ),
    weightsIndex   = 1L,
    optimalWeights = 1
  )
  expect_false( PFIM:::.constraintCellJoint( null_cell$listArms ) )
  expect_identical( PFIM:::.pfimShrinkJointMultMixture( null_cell ), null_cell )

  empty_w = .mock_joint_mult_thinOut()
  empty_w$optimalWeights = numeric()
  expect_identical( PFIM:::.pfimShrinkJointMultMixture( empty_w ), empty_w )
} )

test_that( ".alignOptimalWeightsToArms: needs Arm* names (joint arms are not Arm*)", {
  arms = list(
    .minimal_discrete_arm( "arm1" ),
    .minimal_discrete_arm( "arm2" )
  )
  expect_error(
    PFIM:::.alignOptimalWeightsToArms( arms, c( 1L, 2L ), c( 0.4, 0.6 ) )
  )
} )