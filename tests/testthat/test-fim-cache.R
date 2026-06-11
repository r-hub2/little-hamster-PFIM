# FIM design cache (Phase 2).

local_pfim_opts <- function( opts, env = parent.frame() ) {
  old = setNames( lapply( names( opts ), pfim_get_option ), names( opts ) )
  do.call( pfim_set_option, opts )
  withr::defer( do.call( pfim_set_option, old ), envir = env )
}

test_that("FIM design cache returns same result without re-running", {
  local_pfim_opts( list( fim.cache = TRUE, fim.cache.hits = 0L ) )

  opt = .minimal_mult_opt()
  PFIM:::.pfimFimCacheBegin( opt )
  ev  = PFIM:::.evaluationFromProject( opt )

  r1 = PFIM:::.pfimRunEvaluationCached( ev )
  r2 = PFIM:::.pfimRunEvaluationCached( ev )

  expect_equal( getDeterminant( r1 ), getDeterminant( r2 ) )
  expect_gte( pfim_get_option( "fim.cache.hits" ), 1L )
} )

test_that("optimization reuses FIM cache for final design evaluations", {
  local_pfim_opts( list( fim.cache = TRUE ) )

  opt = .minimal_mult_opt()
  opt = run( opt )

  stats = PFIM:::.pfimFimCacheStats()
  expect_gte( stats$hits, 1L )
  expect_gt( stats$size, 0L )
} )

test_that("constraint grid uses batch model rebuild", {
  local_pfim_opts( list( fim.cache = TRUE ) )

  opt = .minimal_mult_opt()
  invisible( generateFimsFromConstraints( opt ) )

  expect_true( isTRUE( pfim_get_option( "eval.batch", FALSE ) ) ||
    !is.null( pfim_get_option( "fim.cache.scope" ) ) )
  stats = PFIM:::.pfimFimCacheStats()
  expect_gt( stats$size, 0L )
  expect_false( is.null( pfim_get_option( "fim.cache.scope" ) ) )
} )

test_that("fim.cache = FALSE bypasses cache", {
  local_pfim_opts( list( fim.cache = FALSE, fim.cache.hits = 0L ) )

  opt = .minimal_mult_opt()
  PFIM:::.pfimFimCacheBegin( opt )
  ev = PFIM:::.evaluationFromProject( opt )

  PFIM:::.pfimRunEvaluationCached( ev )
  PFIM:::.pfimRunEvaluationCached( ev )

  expect_equal( pfim_get_option( "fim.cache.hits" ), 0L )
} )

test_that("constraint grid respects constraints.maxTasks", {
  local_pfim_opts( list( constraints.maxTasks = 3L ) )

  total = 20L
  idx = PFIM:::.pfimConstraintTaskIndices( total, numberOfDoses = 4L, nCombinations = 5L )
  expect_equal( length( idx ), 3L )
  expect_true( all( idx >= 1L & idx <= total ) )
} )

test_that(".evaluationFromOptimization matches project helper fields", {
  opt = .minimal_mult_opt()
  design = pluck( projectProp( opt, "designs" ), 1L )
  ev = PFIM:::.evaluationFromOptimization( opt, design, name = "test" )
  expect_s7_class( ev, Evaluation )
  expect_equal( prop( ev, "name" ), "test" )
  expect_equal( prop( ev, "fimType" ), projectProp( opt, "fimType" ) )
} )
