# covariateTest — significance (covariate/parameter), TOST, relevance.

test_that( "covariateTest significance helpers", {
  z = qnorm( 0.975 )
  expect_equal( .powerSignificance( 0.5, 0.1, z ), 1, tolerance = 0.002 )
  expect_lt( .powerSignificance( 0.05, 0.1, z ), 0.5 )
  expect_true( is.na( .nRequiredSignificance( 0, 1, z, 0.9 ) ) )
  expect_gt( .nRequiredSignificance( 0.5, 1, z, 0.9 ), 2 )
} )

test_that( "tost is an alias for covariateTest", {
  expect_identical( tost, covariateTest )
} )

test_that( "covariateTest exposes thetaL, thetaU, tests and three beta tests", {
  evaluation = run( cas10_evaluation( "ct_pop", "population" ) )
  ct         = covariateTest( evaluation, thetaL = log( 0.80 ), thetaU = log( 1.25 ) )

  expect_true( "thetaL" %in% names( formals( covariateTest ) ) )
  expect_true( "tests" %in% names( formals( covariateTest ) ) )
  expect_gt( nrow( prop( ct, "covariate" ) ), 0L )
  expect_gt( nrow( prop( ct, "nonRelevance" ) ) +
               nrow( prop( ct, "relevance" ) ), 0L )
  expect_false( any( startsWith( prop( ct, "covariate" )$Parameter, "\u03bc_" ) ) )
  expect_equal( names( prop( ct, "covariate" ) ),
                c( "Parameter", "Value", "SE", "RSE", "Power", "N_Required" ) )
  expect_false( "Status" %in% names( prop( ct, "nonRelevance" ) ) )
  expect_false( "Power_max" %in% names( prop( ct, "nonRelevance" ) ) )
} )

test_that( "covariateTest tests argument selects output slots", {
  evaluation = run( cas10_evaluation( "ct_pick", "population" ) )
  ct         = covariateTest( evaluation, tests = "nonRelevance" )
  expect_equal( nrow( prop( ct, "covariate" ) ), 0L )
  expect_gt( nrow( prop( ct, "nonRelevance" ) ), 0L )
  expect_equal( prop( ct, "settings" )$tests, "nonRelevance" )
} )

test_that( "covariateTest on cas3 population: beta significance only (no mu)", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[1L]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  ct         = covariateTest( evaluation )
  cov        = prop( ct, "covariate" )
  par        = prop( ct, "parameter" )

  expect_gt( nrow( cov ), 0L )
  expect_equal( nrow( par ), 0L )
  expect_false( any( startsWith( cov$Parameter, "\u03bc_" ) ) )
  expect_true( all( startsWith( cov$Parameter, "\u03b2_" ) | startsWith( cov$Parameter, "beta_" ) ) )
} )

test_that( "getCovariateTestTables returns kable slots", {
  evaluation = run( cas10_evaluation( "ct_tables", "population" ) )
  tables     = getCovariateTestTables( covariateTest( evaluation ) )

  expect_true( "significance" %in% names( tables ) )
  expect_true( "nonRelevance" %in% names( tables ) )
  expect_null( tables$parameter )
  expect_s3_class( tables$significance, "knitr_kable" )
} )

test_that( "covariateTest individual has no beta rows (N_scale = 1)", {
  evaluation = run( cas10_evaluation( "ct_nscale", "individual" ) )
  ct         = covariateTest( evaluation, tests = "significance" )
  st         = prop( ct, "settings" )
  expect_equal( st$N0, 40 )
  expect_equal( st$N_scale, 1 )
  expect_equal( st$fimType, "individual" )
  expect_equal( nrow( prop( ct, "covariate" ) ), 0L )
  expect_false( isTRUE( st$additive ) )
} )

test_that( "covariateTest population uses N_scale = N0", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  ct         = covariateTest( evaluation, tests = "significance" )
  st         = prop( ct, "settings" )
  expect_equal( st$N_scale, st$N0 )
  expect_equal( st$fimType, "population" )
} )

test_that( "covariateTest additive link uses (1+beta) as Ratio", {
  sc = PFIM:::.covariateTestEffectScale(
    values = 0.2, SE = 0.05, z_one = qnorm( 0.95 ),
    thetaL = log( 0.8 ), thetaU = log( 1.25 ), additive = TRUE
  )
  expect_equal( sc$effect, 1.2 )
  expect_equal( sc$Binf, exp( log( 0.8 ) ) - 1 )
  expect_equal( sc$Bsup, exp( log( 1.25 ) ) - 1 )
  sc_exp = PFIM:::.covariateTestEffectScale(
    values = log( 1.2 ), SE = 0.05, z_one = qnorm( 0.95 ),
    thetaL = log( 0.8 ), thetaU = log( 1.25 ), additive = FALSE
  )
  expect_equal( sc_exp$effect, 1.2, tolerance = 1e-12 )
} )

test_that( "categoriesProportions must sum to 1", {
  expect_error(
    CategoricalCovariate(
      name = "Sex",
      categories = c( "M", "F" ),
      categoriesProportions = c( 0.6, 0.3 )
    ),
    "categoriesProportions must sum to 1"
  )
} )

test_that( "sequencesProportions must sum to 1", {
  expect_error(
    CategoricalCovariateWithIOV(
      name = "Trt",
      categories = c( "A", "B" ),
      sequences = list( c( "A", "B" ), c( "B", "A" ) ),
      sequencesProportions = c( 0.7, 0.2 )
    ),
    "sequencesProportions must sum to 1"
  )
} )

test_that( "flat sort by groups keeps outcome groups separate", {
  row = c( 8, 2, 50, 0.5 )
  groups = list( 1:2, 3:4 )
  out = PFIM:::.pfimSortFlatByGroups( row, groups )
  expect_equal( out[ 1:2 ], c( 2, 8 ) )
  expect_equal( out[ 3:4 ], c( 0.5, 50 ) )
} )

test_that( "saveCovariateTest returns the output path invisibly", {
  evaluation = run( cas10_evaluation( "ct_save", "population" ) )
  ct = covariateTest( evaluation, tests = "significance" )
  f = tempfile( fileext = ".txt" )
  on.exit( unlink( f ), add = TRUE )
  path = suppressMessages( saveCovariateTest( ct, f ) )
  expect_equal( path, f )
  expect_true( file.exists( f ) )
} )
