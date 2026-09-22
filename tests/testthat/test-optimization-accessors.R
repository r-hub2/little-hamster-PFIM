# Optimization accessors, show(), and arm constraint tables.

test_that("show(Optimization) runs after multiplicative optimization", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_mult_opt() )
  expect_output( show( opt ), "Optimal design" )
  expect_output( show( opt ), "Fisher Matrix" )
})

test_that("getMixtureDcriterion and getRealisedDcriterion require Multiplicative", {
  .skip_optimizer_run_on_cran()
  opt_mult = run( .minimal_mult_opt() )
  expect_true( is.finite( getMixtureDcriterion( opt_mult ) ) )
  expect_equal( getRealisedDcriterion( opt_mult ), getDcriterion( opt_mult ), tolerance = 1e-8 )
  opt_fw = suppressWarnings(
    run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_dc_gate" ) )
  )
  expect_error( getMixtureDcriterion( opt_fw ), "MultiplicativeAlgorithm" )
  expect_error( getRealisedDcriterion( opt_fw ), "MultiplicativeAlgorithm" )
})

test_that("plotWeights and getArmConstraints work for MultiplicativeAlgorithm", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_mult_opt() )
  algo = .opt_algo( opt )
  arm  = .first_arm( opt )
  expect_s3_class( plotWeights( opt ), "ggplot" )
  cons = PFIM:::getArmConstraints( arm, algo )
  expect_type( cons, "list" )
  expect_gt( length( cons ), 0L )
})

test_that("plotFrequencies works for FedorovWynnAlgorithm", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_plot" ) )
  algo = .opt_algo( opt )
  arm  = .first_arm( opt )
  expect_type( PFIM:::getArmConstraints( arm, algo ), "list" )
  expect_s3_class( plotFrequencies( opt ), "ggplot" )
})

test_that(".pfimAsArm unwraps FW entry-shaped optimalArms", {
  arm = Arm( name = "Arm1", size = 10 )
  entry = list( arm = arm, samplingsForFW = c( 0.25, 2, 4, 6 ) )
  expect_identical( PFIM:::.pfimAsArm( arm ), arm )
  expect_identical( PFIM:::.pfimAsArm( entry ), arm )
  expect_error( PFIM:::.pfimAsArm( list( samplingsForFW = 1 ) ), "Expected an Arm" )
})

test_that("plotWeights and plotFrequencies require matching optimizer class", {
  .skip_optimizer_run_on_cran()
  opt_mult = run( .minimal_mult_opt() )
  opt_fw   = suppressWarnings(
    run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_gate" ) )
  )
  expect_error( plotFrequencies( opt_mult ), "plotFrequencies" )
  expect_error( plotWeights( opt_fw ), "plotWeights" )
})

test_that(".pfimDiscreteMixturePlotData uses Protocol labels on length mismatch", {
  arm1 = Arm( name = "a", size = 10 )
  arm2 = Arm( name = "b", size = 10 )
  opt  = Optimization(
    name = "plot_mismatch",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
      ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
    ),
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 ) ),
    optimizer = "MultiplicativeAlgorithm",
    optimizerParameters = .mult_test_algo_params(),
    designs = list( Design( name = "d", arms = list( arm1 ) ) ),
    fimType = "individual",
    fim = IndividualFim(),
    outputs = list( "RespPK" )
  )
  prop( opt, "optimisationAlgorithmOutputs" ) = list(
    optimalArms    = list( arm1, arm2 ),
    optimalWeights = c( 0.7, 0.2, 0.1 )
  )
  dat = PFIM:::.pfimDiscreteMixturePlotData( opt )
  expect_equal( dat$xlab, "Protocol" )
  expect_equal( dat$data$label, paste0( "Protocol", 1:3 ) )
  expect_equal( dat$data$value, c( 0.7, 0.2, 0.1 ) )
})

test_that("getArmConstraints works for PSOAlgorithm", {
  .skip_optimizer_run_on_cran()
  opt  = run( .minimal_continuous_opt( "PSOAlgorithm", name = "pso_cons" ) )
  algo = .opt_algo( opt )
  arm  = .first_arm( opt )
  cons = PFIM:::getArmConstraints( arm, algo )
  expect_type( cons, "list" )
  expect_gt( length( cons ), 0L )
})

test_that("arm administration and arm data helpers return tabular fields", {
  arm = .first_arm( .minimal_mult_opt() )
  admin = PFIM:::armAdministration( arm, "design1" )
  expect_equal( admin[[ 1L ]][[ "Design name" ]], "design1" )
  expect_type( admin, "list" )
  data  = PFIM:::getArmData( arm )
  expect_type( data, "list" )
  expect_true( "Sampling times" %in% names( data[[ 1L ]] ) )
})

test_that("Optimization getters return valid structures after run", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_mult_opt() )
  se = getSE( opt )
  expect_true( is.data.frame( se ) || is.matrix( se ) )
  expect_true( nrow( se ) >= 1L )
  expect_true( is.matrix( getCorrelationMatrix( opt ) ) )
  fm = getFisherMatrix( opt )
  expect_true( is.matrix( fm$fisherMatrix ) )
  expect_true( is.matrix( prop( prop( opt, "fim" ), "fisherMatrix" ) ) )
  expect_equal(
    PFIM:::Dcriterion( prop( opt, "fim" ) ),
    getDcriterion( opt ),
    tolerance = 1e-10
  )
})

test_that("getCorrelationMatrix uses inverse FIM (not cor of FIM columns)", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_mult_opt() )
  M   = .fim_matrix( PFIM:::.getOptimalFim( opt ) )
  C   = getCorrelationMatrix( opt )
  expect_equal( C, PFIM:::.fimCorrelationMatrix( M ), tolerance = 1e-10 )
  expect_equal( diag( C ), rep( 1, nrow( C ) ), tolerance = 1e-10 )
})
