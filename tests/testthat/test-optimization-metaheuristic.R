# PSO/PGBO batch fitness, convergence API, Fedorov singular mixtures.

test_that(".pfimMetaheuristicFitnessBatch returns one value per row", {
  old = pfim_get_option( "fim.cache" )
  pfim_set_option( fim.cache = TRUE )
  on.exit( pfim_set_option( fim.cache = old ), add = TRUE )

  opt = .minimal_continuous_opt( "PSOAlgorithm", name = "pso_batch_test" )
  .pfimFimCacheBegin( opt )
  design = pluck( projectProp( opt, "designs" ), 1L )
  arms   = prop( design, "arms" )
  layout = .buildFlatSamplingLayout( design )
  ev     = .evaluationFromOptimization( opt, design, name = "" )

  flat = layout$initial_flat
  mat  = rbind( flat, flat )
  out  = PFIM:::.pfimMetaheuristicFitnessBatch( ev, design, arms, mat, layout = layout )

  expect_length( out, 2L )
  expect_true( all( is.finite( out ) ) )
  expect_equal( out[ 1L ], out[ 2L ] )
})

test_that("PSO batch callback integrates with pso_optimize_Rcpp", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_continuous_opt( "PSOAlgorithm", name = "pso_batch_cpp" )
  opt = run( opt )
  expect_gt( getDcriterion( opt ), 0 )
})

test_that( "PSOAlgorithm exposes converged flag in optimizer outputs", {
  .skip_optimizer_run_on_cran()
  opt = suppressWarnings( run( .minimal_continuous_opt( "PSOAlgorithm", name = "pso_conv_api" ) ) )
  out = .opt_algo_status( opt )
  expect_type( out, "list" )
  expect_true( all( c( "converged", "iterations", "improved" ) %in% names( out ) ) )
  # Default tolerance=0 → stall stop disabled → converged is NA (not FALSE).
  expect_true( is.na( out$converged ) )
} )

test_that( "PGBOAlgorithm exposes converged flag in optimizer outputs", {
  .skip_optimizer_run_on_cran()
  opt = suppressWarnings( run( .minimal_continuous_opt( "PGBOAlgorithm", name = "pgbo_conv_api" ) ) )
  out = .opt_algo_status( opt )
  expect_type( out, "list" )
  expect_true( all( c( "converged", "iterations", "improved", "bestD", "initialD" ) %in% names( out ) ) )
  expect_true( is.na( out$converged ) )
  expect_gte( out$bestD, out$initialD )
} )

test_that( "PSO can converge early when tolerance is set", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_continuous_opt(
    "PSOAlgorithm",
    name = "pso_early_conv",
    optimizerParameters = list(
      maxIteration = 50L,
      populationSize = 5L,
      seed = 42L,
      personalLearningCoefficient = 2.05,
      globalLearningCoefficient = 2.05,
      tolerance = 1e-2,
      stallIterations = 2L,
      showProcess = FALSE
    )
  )
  opt = suppressWarnings( run( opt ) )
  out = .opt_algo_status( opt )
  expect_true( isTRUE( out$converged ) || isFALSE( out$converged ) )
  expect_true( out$iterations <= 50L )
} )

test_that( "SimplexAlgorithm exposes algorithmOutput; NA when tol disabled", {
  .skip_optimizer_run_on_cran()
  opt_na = suppressWarnings( run( .minimal_continuous_opt(
    "SimplexAlgorithm",
    name = "simplex_api_na",
    optimizerParameters = list(
      pctInitialSimplexBuilding = 20,
      maxIteration = 15L,
      tolerance = 0,
      showProcess = FALSE
    )
  ) ) )
  out_na = .opt_algo_status( opt_na )
  expect_true( all( c( "converged", "iterations", "bestCost" ) %in% names( out_na ) ) )
  expect_true( is.na( out_na$converged ) )

  opt = suppressWarnings( run( .minimal_continuous_opt( "SimplexAlgorithm", name = "simplex_api" ) ) )
  out = .opt_algo_status( opt )
  expect_true( isTRUE( out$converged ) || isFALSE( out$converged ) )
} )

test_that( "Fedorov-Wynn handles singular-mixture cold start without C++ errors", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_singular_smoke" )
  expect_no_warning( suppressMessages( PFIM:::generateFimsFromConstraints( opt ) ) )
  expect_no_warning( suppressWarnings( run( opt ) ) )
} )
