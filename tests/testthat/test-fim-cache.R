# FIM design cache.

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
  invisible( PFIM:::generateFimsFromConstraints( opt ) )

  stats = PFIM:::.pfimFimCacheStats()
  expect_gt( stats$size, 0L )
} )

test_that("constraint grid uses batch model rebuild", {
  local_pfim_opts( list( fim.cache = TRUE ) )

  opt = .minimal_mult_opt()
  invisible( PFIM:::generateFimsFromConstraints( opt ) )

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
  expect_equal( idx, c( 1L, 6L, 11L ) )
} )

test_that("constraint subsample is independent of the global RNG", {
  local_pfim_opts( list( constraints.maxTasks = 3L ) )
  idx1 = PFIM:::.pfimConstraintTaskIndices( 20L, 4L, 5L )
  set.seed( 999L )
  idx2 = PFIM:::.pfimConstraintTaskIndices( 20L, 4L, 5L )
  expect_equal( idx1, idx2 )
} )

test_that("fim.cache.maxEntries trims oldest entries (LRU)", {
  local_pfim_opts( list( fim.cache.maxEntries = 2L ) )
  env = PFIM:::.pfimFimDesignCache
  rm( list = ls( env, all.names = TRUE ), envir = env )
  PFIM:::.pfimCacheStore( env, "a", 1L )
  PFIM:::.pfimCacheStore( env, "b", 2L )
  expect_equal( PFIM:::.pfimCacheKeys( env ), c( "a", "b" ) )
  PFIM:::.pfimCacheStore( env, "c", 3L )
  expect_equal( PFIM:::.pfimCacheKeys( env ), c( "b", "c" ) )
  PFIM:::.pfimCacheLruTouch( env, "b" )
  PFIM:::.pfimCacheStore( env, "d", 4L )
  expect_equal( PFIM:::.pfimCacheKeys( env ), c( "b", "d" ) )
} )

test_that(".evaluationFromOptimization matches project helper fields", {
  opt = .minimal_mult_opt()
  design = pluck( projectProp( opt, "designs" ), 1L )
  ev = PFIM:::.evaluationFromOptimization( opt, design, name = "test" )
  expect_s7_class( ev, Evaluation )
  expect_equal( prop( ev, "name" ), "test" )
  expect_equal( prop( ev, "fimType" ), projectProp( opt, "fimType" ) )
} )

.infusion_fim_cache_ev = function( tinf = 2 ) {
  admin = Administration( outcome = "RespPK", Tinf = tinf, timeDose = 0, dose = 30 )
  arm = Arm(
    name = "a1",
    size = 40,
    administrations = list( admin ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 8 ) )
    )
  )
  Evaluation(
    name = "cache_tinf",
    modelFromLibrary = list( PKModel = "Linear1InfusionSingleDose_ClV" ),
    modelParameters = list(
      ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ) ),
      ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ) )
    ),
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    designs = list( Design( name = "d1", arms = list( arm ) ) ),
    fimType = "population",
    outputs = list( "RespPK" = "RespPK" )
  )
}

test_that("design signature includes Tinf and tau", {
  d1 = prop( .infusion_fim_cache_ev( 2 ), "designs" )[[ 1L ]]
  d2 = prop( .infusion_fim_cache_ev( 8 ), "designs" )[[ 1L ]]
  expect_false( identical(
    PFIM:::.pfimDesignSignature( d1 ),
    PFIM:::.pfimDesignSignature( d2 )
  ) )
} )

test_that("design signature includes arm size (population cache)", {
  ev  = .infusion_fim_cache_ev( 2 )
  d1  = prop( ev, "designs" )[[ 1L ]]
  arm = prop( d1, "arms" )[[ 1L ]]
  d2  = PFIM:::.pfimCloneS7( d1 )
  prop( prop( d2, "arms" )[[ 1L ]], "size" ) = 75
  expect_false( identical(
    PFIM:::.pfimDesignSignature( d1 ),
    PFIM:::.pfimDesignSignature( d2 )
  ) )
} )

test_that("design signature includes ODE initialConditions", {
  mk = function( ic ) {
    Arm(
      name = "a", size = 10,
      administrations = list(
        Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
      ),
      samplingTimes = list(
        SamplingTimes( outcome = "RespPK", samplings = c( 1, 2 ) )
      ),
      initialConditions = ic
    )
  }
  d1 = Design( name = "d", arms = list( mk( list( C1 = 0 ) ) ) )
  d2 = Design( name = "d", arms = list( mk( list( C1 = 1 ) ) ) )
  expect_false( identical(
    PFIM:::.pfimDesignSignature( d1 ),
    PFIM:::.pfimDesignSignature( d2 )
  ) )
} )

test_that("FIM design cache invalidates when Tinf changes", {
  local_pfim_opts( list( fim.cache = TRUE, fim.cache.hits = 0L ) )
  PFIM:::.pfimClearFimCaches()

  ev2  = .infusion_fim_cache_ev( 2 )
  ev10 = .infusion_fim_cache_ev( 10 )
  pfim_set_option( fim.cache.scope = PFIM:::.pfimProjectScopeId( projectOf( ev2 ) ) )

  r1 = PFIM:::.pfimRunEvaluationCached( ev2 )
  r2 = PFIM:::.pfimRunEvaluationCached( ev10 )

  expect_false( isTRUE( all.equal(
    getDeterminant( r1 ), getDeterminant( r2 ), tolerance = 1e-10
  ) ) )

  PFIM:::.pfimRunEvaluationCached( ev2 )
  expect_gte( pfim_get_option( "fim.cache.hits" ), 1L )
} )

test_that( "FIM cache keys include model signature", {
  ev1 = .infusion_fim_cache_ev( 2 )
  ev2 = .infusion_fim_cache_ev( 2 )
  params = prop( ev2, "modelParameters" )
  params[[ 2L ]] = PFIM::ModelParameter(
    name = "Cl",
    distribution = PFIM::LogNormal( mu = 2.3, omega = sqrt( 0.09 ) )
  )
  prop( ev2, "modelParameters" ) = params
  expect_false( identical(
    PFIM:::.pfimFimCacheKey( ev1 ),
    PFIM:::.pfimFimCacheKey( ev2 )
  ) )
} )

test_that( "FIM design cache misses when model parameters change", {
  local_pfim_opts( list( fim.cache = TRUE, fim.cache.hits = 0L ) )
  PFIM:::.pfimClearFimCaches()

  ev1 = .infusion_fim_cache_ev( 2 )
  ev2 = .infusion_fim_cache_ev( 2 )
  params = prop( ev2, "modelParameters" )
  params[[ 2L ]] = PFIM::ModelParameter(
    name = "Cl",
    distribution = PFIM::LogNormal( mu = 2.3, omega = sqrt( 0.09 ) )
  )
  prop( ev2, "modelParameters" ) = params

  pfim_set_option( fim.cache.scope = PFIM:::.pfimProjectScopeId( projectOf( ev1 ) ) )
  d1 = getDeterminant( PFIM:::.pfimRunEvaluationCached( ev1 ) )
  d2 = getDeterminant( PFIM:::.pfimRunEvaluationCached( ev2 ) )
  expect_false( isTRUE( all.equal( d1, d2, tolerance = 1e-8 ) ) )
} )

test_that( ".pfimFimCacheClear(NULL) is safe", {
  local_pfim_opts( list( fim.cache.scope = NULL ) )
  expect_silent( PFIM:::.pfimFimCacheClear( NULL ) )
} )

test_that( "run() clears cache evaluation context on exit", {
  ev = .infusion_fim_cache_ev( 2 )
  invisible( run( ev ) )
  expect_null( pfim_get_option( "fim.cache.evaluation", NULL ) )
} )

test_that( ".pfimProjectScopeId is stable for one Evaluation", {
  PFIM:::.pfimClearProjectScopeRegistry()
  ev = .infusion_fim_cache_ev( 2 )
  s1 = PFIM:::.pfimProjectScopeId( ev )
  s2 = PFIM:::.pfimProjectScopeId( ev )
  expect_equal( s1, s2 )
} )

test_that( "FIM cache keys include covariate effects", {
  ev1 = cas10_evaluation( "cov_sig1", "population" )
  ev2 = cas10_evaluation( "cov_sig2", "population" )
  covs = prop( ev2, "modelCovariates" )
  covs[[ 1L ]] = CategoricalCovariate(
    name = "Sex",
    categories = c( "M", "F" ),
    categoriesProportions = c( 0.5, 0.5 ),
    effects = list( "F" = c( "V" = log( 1.5 ) ) )
  )
  prop( ev2, "modelCovariates" ) = covs
  expect_false( identical(
    PFIM:::.pfimFimCacheKey( ev1 ),
    PFIM:::.pfimFimCacheKey( ev2 )
  ) )
} )

test_that( "FIM cache keys include all designs", {
  ev = .infusion_fim_cache_ev( 2 )
  d1 = prop( ev, "designs" )[[ 1L ]]
  d2 = PFIM:::.pfimCloneS7( d1 )
  prop( d2, "name" ) = "design2"
  prop( ev, "designs" ) = list( d1, d2 )
  ev1 = PFIM:::.pfimCloneS7( ev )
  prop( ev1, "designs" ) = list( d1 )
  expect_false( identical(
    PFIM:::.pfimFimCacheKey( ev ),
    PFIM:::.pfimFimCacheKey( ev1 )
  ) )
} )

test_that( "constraint grid registers distinct design-cache entries per cell", {
  local_pfim_opts( list( fim.cache = TRUE, constraints.maxTasks = NULL ) )
  PFIM:::.pfimEnvClear( PFIM:::.pfimFimDesignCache )
  opt = .minimal_mult_opt()
  PFIM:::.pfimFimCacheBegin( opt )
  invisible( PFIM:::generateFimsFromConstraints( opt ) )
  keys = PFIM:::.pfimCacheKeys( PFIM:::.pfimFimDesignCache )
  n_comb = length( generateSamplingTimesCombination(
    projectProp( opt, "designs" )[[ 1L ]]
  )$opt_arm )
  expect_equal( length( keys ), n_comb )
  samplings = lapply( keys, function( key ) {
    ev  = get( key, envir = PFIM:::.pfimFimDesignCache )
    arm = prop( prop( ev, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
    sort( prop( prop( arm, "samplingTimes" )[[ 1L ]], "samplings" ) )
  } )
  expect_equal( length( unique( samplings ) ), length( keys ) )
} )

test_that( ".pfimEvaluationFromDesign clones evaluation template", {
  local_pfim_opts( list( fim.cache = FALSE ) )
  template = .infusion_fim_cache_ev( 2 )
  ev_done  = run( template )
  d1 = prop( template, "designs" )[[ 1L ]]
  d2 = PFIM:::.pfimCloneS7( d1 )
  prop( d2, "name" ) = "design_variant"
  r1 = PFIM:::.pfimEvaluationFromDesign( template, d1, ev_done )
  template2 = .infusion_fim_cache_ev( 2 )
  r2 = PFIM:::.pfimEvaluationFromDesign( template2, d2, ev_done )
  expect_equal( length( prop( template, "designs" ) ), 1L )
  expect_equal( prop( prop( r1, "designs" )[[ 1L ]], "name" ), "d1" )
  expect_equal( prop( prop( r2, "designs" )[[ 1L ]], "name" ), "design_variant" )
} )
