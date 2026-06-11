# End-to-end MultiplicativeAlgorithm optimization on a small constraint grid.

test_that("MultiplicativeAlgorithm: run() fills optimal design and improves D-criterion", {
  opt = .minimal_mult_opt()

  combos = generateSamplingTimesCombination( projectProp( opt, "designs" )[[ 1L ]] )
  expect_gte( length( combos$opt_arm ), 2L )

  opt = run( opt )

  od = prop( opt, "optimisationDesign" )
  expect_type( od, "list" )
  expect_s7_class( od$evaluationInitialDesign, Evaluation )
  expect_s7_class( od$evaluationOptimalDesign, Evaluation )

  det_init = getDeterminant( od$evaluationInitialDesign )
  det_opt  = getDeterminant( od$evaluationOptimalDesign )
  expect_true( is.finite( det_init ) && det_init > 0 )
  expect_true( is.finite( det_opt ) && det_opt > 0 )

  weights = prop( opt, "optimisationAlgorithmOutputs" )$optimalWeights
  expect_true( length( weights ) >= 1L )
  expect_true( all( weights > 0 ) )

  expect_gt( getDcriterion( opt ), 0 )
})
