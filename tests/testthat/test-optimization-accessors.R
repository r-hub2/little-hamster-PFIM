# Optimization accessors, show(), and arm constraint tables.

test_that("show(Optimization) runs after multiplicative optimization", {
  opt = run( .minimal_mult_opt() )
  expect_output( show( opt ), "Initial design" )
  expect_output( show( opt ), "Optimal design" )
})

test_that("plotWeights and getArmConstraints work for MultiplicativeAlgorithm", {
  opt = run( .minimal_mult_opt() )
  algo = prop( opt, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  arm  = prop( projectProp( opt, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  expect_s3_class( plotWeights( opt ), "ggplot" )
  cons = PFIM:::getArmConstraints( arm, algo )
  expect_type( cons, "list" )
  expect_gt( length( cons ), 0L )
})

test_that("plotFrequencies works for FedorovWynnAlgorithm", {
  opt = run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_plot" ) )
  algo = prop( opt, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  arm  = prop( projectProp( opt, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  expect_type( PFIM:::getArmConstraints( arm, algo ), "list" )
  expect_s3_class( plotFrequencies( opt ), "ggplot" )
})

test_that("getArmConstraints works for PSOAlgorithm", {
  opt  = run( .minimal_continuous_opt( "PSOAlgorithm", name = "pso_cons" ) )
  algo = prop( opt, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  arm  = prop( projectProp( opt, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  cons = PFIM:::getArmConstraints( arm, algo )
  expect_type( cons, "list" )
  expect_gt( length( cons ), 0L )
})

test_that("arm administration and arm data helpers return tabular fields", {
  arm = prop( projectProp( .minimal_mult_opt(), "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  admin = PFIM:::armAdministration( arm )
  expect_type( admin, "list" )
  data  = PFIM:::getArmData( arm )
  expect_type( data, "list" )
  expect_true( "Sampling times" %in% names( data[[ 1L ]] ) )
})

test_that("Optimization getters return valid structures after run", {
  opt = run( .minimal_mult_opt() )
  se = getSE( opt )
  expect_true( is.data.frame( se ) || is.matrix( se ) )
  expect_true( nrow( se ) >= 1L )
  expect_true( is.matrix( getCorrelationMatrix( opt ) ) )
  fm = getFisherMatrix( opt )
  expect_true( is.matrix( fm$fisherMatrix ) )
})
