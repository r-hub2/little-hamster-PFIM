# registries v3 + project access

defer_rm = function( key, env, envir = parent.frame() ) {
  withr::defer( rm( list = key, envir = env, inherits = FALSE ), envir = envir )
}

.pkg_root = function() {
  normalizePath( file.path( testthat::test_path(), "..", ".." ), mustWork = FALSE )
}

.pkg_source_r_lines = function( filename ) {
  path = file.path( .pkg_root(), "R", filename )
  skip_if_not( file.exists( path ), "R sources unavailable (installed-package check)" )
  readLines( path, warn = FALSE )
}

test_that( "fim type registry", {
  pfim_register_fim_type( "testfim", function() IndividualFim() )
  defer_rm( "testfim", PFIM:::.pfimFimTypeRegistry )
  out = defineFim( PFIMProject( fimType = "testfim" ) )
  expect_true( S7::S7_inherits( out, IndividualFim ) )
} )

test_that( "optimizer registry hook", {
  pfim_register_optimizer( "TestOptimizer", function() MultiplicativeAlgorithm() )
  defer_rm( "TestOptimizer", PFIM:::.pfimOptimizerRegistry )
  opt = Optimization(
    name = "x", optimizer = "TestOptimizer",
    fimType = "individual", fim = IndividualFim(), designs = list()
  )
  expect_s7_class( defineOptimizationAlgorithm( opt ), MultiplicativeAlgorithm )
})

test_that("getEvaluationDesign by name or index", {
  d1 = Design( name = "d1", arms = list(), fim = IndividualFim() )
  d2 = Design( name = "d2", arms = list(), fim = IndividualFim() )
  ev = Evaluation(
    name = "ev", fimType = "individual", fim = IndividualFim(),
    designs = list( d1, d2 ), evaluationDesign = list( d1, d2 )
  )
  expect_equal( prop( getEvaluationDesign( ev, 2 ), "name" ), "d2" )
  expect_equal( prop( getEvaluationDesign( ev, "d1" ), "name" ), "d1" )
})

test_that( "modelClass bypasses detect()", {
  proj = PFIMProject(
    modelClass = "ModelAnalytic",
    modelEquations = list( weird = "no Deriv_ here" ),
    modelParameters = list(), modelError = list()
  )
  expect_s7_class( defineModelType( proj ), ModelAnalytic )
} )

test_that("continuous setOptimalArms", {
  arm = Arm(
    name = "a1", size = 1,
    administrations = list( Administration( outcome = "y", timeDose = 0, dose = 1 ) ),
    samplingTimes = list( SamplingTimes( outcome = "y", samplings = c( 1, 2 ) ) )
  )
  algo = SimplexAlgorithm()
  prop( algo, "optimizerOutputs" ) = list( optimalArms = list( arm ) )
  expect_equal( setOptimalArms( IndividualFim(), algo ), list( arm ) )
})

test_that( "reset session clears design cache keys", {
  assign( "zz", 1L, envir = PFIM:::.pfimFimDesignCache )
  pfim_reset_session()
  expect_false( exists( "zz", envir = PFIM:::.pfimFimDesignCache, inherits = FALSE ) )
} )

test_that( "reset session clears eval-model rebuild cache", {
  assign( "zz", list( a = 1L ), envir = PFIM:::.pfimEvalModelCache )
  pfim_reset_session()
  expect_false( exists( "zz", envir = PFIM:::.pfimEvalModelCache, inherits = FALSE ) )
} )

test_that( "reset session clears ODE sim-time cache", {
  pfimOdeSimTimesCacheClear_Rcpp()
  pfimOdeSimTimesCached_Rcpp( "t", c( 1, 2 ), c( 1.5 ), enabled = TRUE )
  expect_equal( pfimOdeSimTimesCacheSize_Rcpp(), 1L )
  pfim_reset_session()
  expect_equal( pfimOdeSimTimesCacheSize_Rcpp(), 0L )
} )

test_that( "evaluateFimConstraintsCell is defined once", {
  skip_if( nzchar( Sys.getenv( "R_COVR" ) ), "source text rewritten under covr" )
  r_dir = file.path( .pkg_root(), "R" )
  skip_if_not( dir.exists( r_dir ), "R sources unavailable (installed-package check)" )
  n_def = sum( vapply(
    list.files( r_dir, "\\.R$", full.names = TRUE ),
    function( f ) sum( grepl( "^\\.evaluateFimConstraintsCell\\s*=", readLines( f, warn = FALSE ) ) ),
    integer( 1L )
  ) )
  expect_equal( n_def, 1L )
} )

test_that( "armAdministration does not hardcode design label", {
  arm_r = .pkg_source_r_lines( "Arm.R" )
  expect_false( any( grepl( "Design optimized", arm_r, fixed = TRUE ) ) )
  expect_true( any( grepl( "designName", arm_r, fixed = TRUE ) ) )
} )

test_that( "CovariateSamplingMethods stub removed", {
  expect_false( file.exists( file.path( .pkg_root(), "R", "CovariateSamplingMethods.R" ) ) )
  desc = readLines( system.file( "DESCRIPTION", package = "PFIM" ), warn = FALSE )
  expect_false( any( grepl( "CovariateSampling", desc, fixed = TRUE ) ) )
} )

test_that( "design FIM cache has LRU cap by default", {
  expect_equal( pfim_get_option( "fim.cache.maxEntries" ), 2048L )
  local_pfim_opts( list( fim.cache.maxEntries = 2L ) )
  env = PFIM:::.pfimFimDesignCache
  PFIM:::.pfimCacheStore( env, "a", 1L )
  PFIM:::.pfimCacheStore( env, "b", 2L )
  PFIM:::.pfimCacheStore( env, "c", 3L )
  expect_equal( length( PFIM:::.pfimCacheKeys( env ) ), 2L )
  rm( list = ls( env, all.names = TRUE ), envir = env )
} )
