# Batch fitness for PSO metaheuristic.

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
} )

test_that("PSO batch callback integrates with pso_optimize_Rcpp", {
  opt = .minimal_continuous_opt( "PSOAlgorithm", name = "pso_batch_cpp" )
  opt = run( opt )
  expect_gt( getDcriterion( opt ), 0 )
} )
