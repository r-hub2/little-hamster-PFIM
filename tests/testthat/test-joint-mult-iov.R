test_that( "joint Mult grid + shrink works with IOV (gamma)", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_discrete_opt_multi_arm(
    n_arms = 2L, fimType = "population", name = "iov_joint_mult"
  )
  params = list(
    ModelParameter(
      name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ), gamma = 0.15
    ),
    ModelParameter(
      name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ), gamma = 0.10
    )
  )
  projectProp( opt, "modelParameters" ) = params
  projectProp( opt, "modelCovariatesEquation" ) = "exponential"
  projectProp( opt, "numberOfOccasions" ) = 2L

  res = PFIM:::generateFimsFromConstraints( opt )
  cells = res$listArms[[ 1L ]]
  expect_true( PFIM:::.constraintCellJoint( cells ) )
  expect_true( all( map_lgl( res$listFimsAlgoMult[[ 1L ]], ~ all( is.finite( .x ) ) ) ) )

  n_keep = min( 2L, length( cells ) )
  thin = list(
    listArms       = cells[ seq_len( n_keep ) ],
    weightsIndex   = seq_len( n_keep ),
    optimalWeights = if ( n_keep == 1L ) 1 else c( 0.35, 0.65 )
  )
  out = PFIM:::.pfimShrinkJointMultMixture( thin )
  expect_equal( length( out$optimalWeights ), 1L )
  expect_equal( out$optimalWeights, 1 )
  expect_length( PFIM:::.constraintCellArms( out$listArms[[ 1L ]] ), 2L )
} )
