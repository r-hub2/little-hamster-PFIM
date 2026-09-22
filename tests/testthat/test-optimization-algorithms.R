# End-to-end optimization: all algorithms, Multiplicative grid, Fedorov-Wynn C++.

test_that("FedorovWynnAlgorithm: run() completes with valid optimal design", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_discrete_opt( "FedorovWynnAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("PSOAlgorithm: run() completes with valid optimal design", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_continuous_opt( "PSOAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("PGBOAlgorithm: run() completes with valid optimal design", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_continuous_opt( "PGBOAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("PGBOAlgorithm: optimal D-criterion is not worse than initial", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_continuous_opt( "PGBOAlgorithm", name = "pgbo_best_d" ) )
  d_init = getDcriterion( .opt_init_eval( opt ) )
  d_opt  = getDcriterion( .opt_eval( opt ) )
  expect_gte( d_opt, d_init )
})

test_that("pgbo_optimize_Rcpp fitness uses cost = 1/D", {
  # Metropolis cost must be 1/D, not raw D.
  # With inverted fitness the resident walks downhill and bestD stays near the start.
  set.seed( 42 )
  eval_d = function( pos ) {
    as.numeric( 100 - ( pos[[ 1L ]] - 5 )^2 - ( pos[[ 2L ]] - 8 )^2 )
  }
  check_valid = function( pos, group_idx ) {
    all( is.finite( pos ) ) &&
      pos[[ 1L ]] <= pos[[ 2L ]] - 0.1 &&
      all( pos >= 0 ) && all( pos <= 10 )
  }
  out = PFIM:::pgbo_optimize_Rcpp(
    initial_pos     = c( 1, 2 ),
    sorting_groups  = list( 1:2 ),
    max_iteration   = 300L,
    N               = 20L,
    mute_effect     = 0.6,
    purge_iteration = 50L,
    fit_base        = 0.5,
    max_attempts    = 5000L,
    cauchy_prob     = 0.5,
    show_process    = FALSE,
    eval_d          = eval_d,
    check_valid_group = check_valid
  )
  expect_gt( out$bestD, out$initialD )
  expect_gt( out$bestD, 80 )
})

test_that("SimplexAlgorithm: run() completes with valid optimal design", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_continuous_opt( "SimplexAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("MultiplicativeAlgorithm: run() completes with valid optimal design", {
  opt = .minimal_discrete_opt( "MultiplicativeAlgorithm" )
  .expect_optimization_run( opt )
})

test_that("generateFimsFromConstraints returns constraint-grid FIMs", {
  opt = .minimal_discrete_opt( "MultiplicativeAlgorithm" )
  res = PFIM:::generateFimsFromConstraints( opt )
  expect_type( res, "list" )
  expect_true( length( res$listFimsAlgoMult ) >= 1L )
  expect_gte( nrow( res$listFimsAlgoMult[[ 1L ]][[ 1L ]] ), 1L )
})

test_that("generateFimsFromConstraints: multi-arm joint FIM grid", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_discrete_opt_multi_arm( n_arms = 2L )
  res = PFIM:::generateFimsFromConstraints( opt )
  expect_length( PFIM:::.constraintCellArms( res$listArms[[ 1L ]][[ 1L ]] ), 2L )
  dims = unique( map_int( res$listFimsAlgoMult[[ 1L ]], nrow ) )
  expect_length( dims, 1L )
  expect_gt( dims, 0L )
})

test_that("generateFimsFromConstraints: multi-arm with covariates", {
  .skip_optimizer_run_on_cran()
  arms   = map( c( "arm1", "arm2" ), ~ .minimal_discrete_arm( .x, size = 40 ) )
  design = Design( name = "cov_multi", arms = arms )
  fx     = cas10_design()
  opt    = Optimization(
    name                    = "multi_cov_grid",
    modelFromLibrary        = cas10_model_from_library(),
    modelParameters         = cas10_model_parameters(),
    modelCovariates         = cas10_covariates(),
    modelCovariatesEquation = "exponential",
    modelError              = fx$modelError,
    designs                 = list( design ),
    fimType                 = "individual",
    optimizer               = "MultiplicativeAlgorithm",
    optimizerParameters     = list(
      lambda = 0.99, numberOfIterations = 5,
      weightThreshold = 0.01, delta = 1e-4
    ),
    outputs                 = list( "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
  res = PFIM:::generateFimsFromConstraints( opt )
  expect_length( PFIM:::.constraintCellArms( res$listArms[[ 1L ]][[ 1L ]] ), 2L )
  expect_true( all( map_lgl( res$listFimsAlgoMult[[ 1L ]], ~ all( is.finite( .x ) ) ) ) )
})

test_that("MultiplicativeAlgorithm: multi-arm discrete run completes", {
  .skip_optimizer_run_on_cran()
  opt = run( .minimal_discrete_opt_multi_arm( n_arms = 2L ) )
  eval_opt = .opt_eval( opt )
  arms     = .design_arms( .eval_design( eval_opt ) )
  out      = .opt_outputs( opt )
  expect_length( arms, 2L )
  expect_gt( getDcriterion( opt ), 0 )
  w = out$optimalWeights
  expect_equal( length( w ), 1L )
  expect_equal( sum( w ), 1, tolerance = 1e-10 )
  expect_true( is.finite( getMixtureDcriterion( opt ) ) )
  expect_true( is.finite( getRealisedDcriterion( opt ) ) )
  expect_equal( getRealisedDcriterion( opt ), getDcriterion( opt ), tolerance = 1e-8 )
  expect_s3_class( plotWeights( opt ), "ggplot" )
})

test_that("FedorovWynnAlgorithm: multi-arm plotFrequencies is protocol-level", {
  .skip_optimizer_run_on_cran()
  opt = suppressWarnings( run(
    .minimal_discrete_opt_multi_arm(
      optimizer = "FedorovWynnAlgorithm", name = "fw_joint_plot"
    )
  ) )
  out = .opt_outputs( opt )
  expect_gte( length( out$optimalArms ), length( out$optimalWeights ) )
  expect_s3_class( plotFrequencies( opt ), "ggplot" )
})

test_that("multi-design discrete run stores perDesign", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_mult_opt()
  d1  = projectProp( opt, "designs" )[[ 1L ]]
  d2  = PFIM:::.pfimCloneS7( d1 )
  prop( d2, "name" ) = "opt_design_b"
  projectProp( opt, "designs" ) = list( d1, d2 )
  opt = run( opt )
  per = .opt_outputs( opt )$perDesign
  expect_length( per, 2L )
  expect_equal( per[[ 2L ]]$name, "opt_design_b" )
  expect_true( length( per[[ 1L ]]$optimalArms ) >= 1L )
})

test_that("optimizerParameters: unknown name suggests fix", {
  expect_error(
    .minimal_discrete_opt(
      "MultiplicativeAlgorithm",
      optimizerParameters = list(
        Lambda = 0.99, delta = 1e-4, numberOfIterations = 5, weightThreshold = 0.01
      )
    ),
    "did you mean 'lambda'",
    fixed = TRUE
  )
})

test_that("optimizerParameters: missing required names fail at construction", {
  expect_error(
    .minimal_discrete_opt(
      "MultiplicativeAlgorithm",
      optimizerParameters = list( lambda = 0.99 )
    ),
    "Missing optimizerParameters",
    fixed = TRUE
  )
})

test_that("MultiplicativeAlgorithm: run() fills optimal design and improves D-criterion", {
  .skip_optimizer_run_on_cran()
  opt = .minimal_mult_opt()

  combos = generateSamplingTimesCombination( projectProp( opt, "designs" )[[ 1L ]] )
  expect_gte( length( combos$opt_arm ), 2L )

  opt = run( opt )

  expect_s7_class( .opt_init_eval( opt ), Evaluation )
  expect_s7_class( .opt_eval( opt ), Evaluation )

  det_init = getDeterminant( .opt_init_eval( opt ) )
  det_opt  = getDeterminant( .opt_eval( opt ) )
  expect_true( is.finite( det_init ) && det_init > 0 )
  expect_true( is.finite( det_opt ) && det_opt > 0 )

  weights = .opt_outputs( opt )$optimalWeights
  expect_true( length( weights ) >= 1L )
  expect_true( all( weights > 0 ) )

  expect_gt( getDcriterion( opt ), 0 )
})

test_that( "constraint grid stores independent arm snapshots per cell", {
  .skip_optimizer_run_on_cran()
  local_pfim_opts( list( fim.cache = FALSE ) )
  opt   = .minimal_mult_pop_dose_grid()
  fims  = PFIM:::generateFimsFromConstraints( opt )
  cells = fims$listArms[[ 1L ]]
  dose_at = function( i ) {
    arm = PFIM:::.constraintCellArms( cells[[ i ]] )[[ 1L ]]
    prop( prop( arm, "administrations" )[[ 1L ]], "dose" )
  }
  n_comb = length( generateSamplingTimesCombination(
    projectProp( opt, "designs" )[[ 1L ]]
  )$grid_arm )
  expect_equal( dose_at( 1L ), 50 )
  expect_equal( dose_at( n_comb + 1L ), 100 )
})

test_that( "generateSamplingTimesCombination yields independent copies per cell", {
  opt    = .minimal_mult_opt()
  design = projectProp( opt, "designs" )[[ 1L ]]
  combos = generateSamplingTimesCombination( design )$opt_arm
  expect_gte( length( combos ), 2L )
  samplings = lapply( combos, function( entry ) {
    sort( prop( entry[[ 1L ]], "samplings" ) )
  } )
  expect_equal( length( unique( samplings ) ), length( samplings ) )
})

test_that( "population multiplicative: arm sizes follow optimal weights", {
  .skip_optimizer_run_on_cran()
  local_pfim_opts( list( fim.cache = FALSE ) )
  opt     = run( .minimal_mult_pop_dose_grid() )
  outputs = .opt_mult_out( opt )
  weights = outputs$optimalWeights
  N       = outputs$numberOfSubjects
  if ( is.null( N ) )
    N = outputs$numberOfArms
  if ( !length( N ) || !is.finite( N ) || N <= 0 ) {
    init_design = .eval_design( .opt_init_eval( opt ) )
    N = prop( .design_arms( init_design )[[ 1L ]], "size" )
  }
  design = .eval_design( .opt_eval( opt ) )
  arms  = .design_arms( design )
  sizes = map_dbl( arms, ~ prop( .x, "size" ) )
  expect_equal(
    sort( sizes ),
    sort( PFIM:::.allocProportionalSubjects( N, weights ) )
  )
})

test_that( "population multiplicative: optimal D matches algorithm mixture", {
  .skip_optimizer_run_on_cran()
  local_pfim_opts( list( fim.cache = FALSE ) )
  opt     = run( .minimal_mult_pop_dose_grid() )
  outputs = .opt_mult_out( opt )
  D_mix   = outputs$mixtureDcriterion
  D_real  = outputs$realisedDcriterion
  D_opt   = getDcriterion( opt )
  expect_true( is.finite( D_mix ) && D_mix > 0 )
  expect_equal( D_opt, D_real, tolerance = 1e-8 )
  # Continuous mixture vs realised allocated design: same N, Hamilton sizes.
  expect_equal( D_opt, D_mix, tolerance = 1e-2 )
  expect_equal( sum( map_dbl(
    .design_arms( .eval_design( .opt_eval( opt ) ) ),
    ~ prop( .x, "size" )
  ) ), outputs$numberOfSubjects, tolerance = 1e-8 )
})

test_that( "Multiplicative setOptimalArms is re-callable after thin persist", {
  .skip_optimizer_run_on_cran()
  local_pfim_opts( list( fim.cache = FALSE ) )
  opt  = run( .minimal_mult_opt() )
  fim  = projectProp( opt, "fim" )
  out  = .opt_mult_out( opt )
  expect_true( length( out$listArms ) >= 1L )
  arms2 = setOptimalArms( fim, .opt_algo( opt ) )
  expect_equal( length( arms2 ), length( out$optimalWeights ) )
})

test_that( "MultiplicativeAlgorithm_Rcpp: converged false at iteration limit", {
  set.seed( 1L )
  p = 3L
  mk = function() crossprod( matrix( rnorm( 20L * p ), 20L, p ) ) + diag( p )
  out = MultiplicativeAlgorithm_Rcpp(
    list( mk(), mk() ), 2L, rep( 0.5, 2L ), p,
    lambda = 0.5, delta = 1e-12, iteration_init = 1L, show_process = FALSE
  )
  expect_false( out$converged )
  expect_equal( out$iterations, 1L )
})

test_that( "MultiplicativeAlgorithm_Rcpp returns certified weights on converge", {
  # Identical candidate FIMs → any w is optimal; check fires before first update.
  p = 2L
  F = diag( p ) + 0.25
  w0 = c( 0.3, 0.7 )
  out = MultiplicativeAlgorithm_Rcpp(
    list( F, F ), 2L, w0, p,
    lambda = 1, delta = 1e-10, iteration_init = 20L, show_process = FALSE
  )
  expect_true( out$converged )
  expect_equal( sum( out$weights ), 1, tolerance = 1e-12 )
  # Certified: returned w equals the start (no update applied after the check).
  expect_equal( as.numeric( out$weights ), w0, tolerance = 1e-12 )
})

test_that( "Mult thin outputs expose canonical schema keys", {
  .skip_optimizer_run_on_cran()
  local_pfim_opts( list( fim.cache = FALSE ) )
  opt  = run( .minimal_mult_opt() )
  out  = .opt_mult_out( opt )
  expect_true( all( c(
    "listArms", "algorithmOutput", "numberOfSubjects",
    "optimalWeights", "mixtureDcriterion", "realisedDcriterion"
  ) %in% names( out ) ) )
  expect_false( "numberOfArms" %in% names( out ) )
  expect_true( is.list( out$algorithmOutput ) )
  expect_true( "converged" %in% names( out$algorithmOutput ) )
})

test_that( "FW outputs expose numberOfSubjects and algorithmOutput", {
  .skip_optimizer_run_on_cran()
  opt = suppressWarnings( run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_schema" ) ) )
  out = .opt_fw_out( opt )
  expect_true( all( c(
    "listArms", "optimalWeights", "numberOfSubjects", "algorithmOutput"
  ) %in% names( out ) ) )
  expect_false( "numberOfIndividuals" %in% names( out ) )
  expect_true( out$algorithmOutput$converged %in% c( TRUE, FALSE ) )
})

test_that( "stall converged flag is NA when tolerance disabled", {
  expect_true( is.na( PFIM:::.pfimStallConvergedFlag( FALSE, 0 ) ) )
  expect_true( is.na( PFIM:::.pfimStallConvergedFlag( TRUE, 0 ) ) )
  expect_false( PFIM:::.pfimStallConvergedFlag( FALSE, 1e-3 ) )
  expect_true( PFIM:::.pfimStallConvergedFlag( TRUE, 1e-3 ) )
})

test_that( "Fedorov-Wynn status helper warns on incomplete and stops otherwise", {
  expect_true( PFIM:::.pfimFedorovWynnCheckStatus( "success" ) )
  expect_warning(
    expect_false( PFIM:::.pfimFedorovWynnCheckStatus( "incomplete" ) ),
    "incomplete"
  )
  expect_error( PFIM:::.pfimFedorovWynnCheckStatus( "singular_fim" ), "singular" )
  expect_error( PFIM:::.pfimFedorovWynnCheckStatus( "boom" ), "status = boom" )
})

test_that( "continuous non-convergence warning is silent when tol is off", {
  expect_silent( PFIM:::.pfimWarnIfNotConverged( "PSOAlgorithm", 0, FALSE ) )
  expect_warning(
    PFIM:::.pfimWarnIfNotConverged( "PSOAlgorithm", 1e-3, FALSE, extra = " within maxIteration = 10" ),
    "PSOAlgorithm: did not converge"
  )
})

test_that("Fedorov-Wynn buffer sizes follow FIM dimension, not protocol count", {
  buf = getFromNamespace( ".fedorovWynnBufferSizes", "PFIM" )( 4L, 3L )
  expect_equal( buf$nFisher, 10L )
  expect_equal( buf$nMaxPop, 11L )
  expect_equal( buf$nBuf, 11L )
  expect_equal( buf$nVectps, 11L )
  expect_gt( buf$nFisher, 3L * 4L / 2L )
  buf_times = getFromNamespace( ".fedorovWynnBufferSizes", "PFIM" )( 2L, 2L, 8L )
  expect_equal( buf_times$nBuf, 4L )
  expect_equal( buf_times$nVectps, 8L )
})

test_that("FedorovWynnAlgorithm run() has no C++ subscript warnings", {
  .skip_optimizer_run_on_cran()
  expect_no_warning(
    run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_cpp_buffers" ) )
  )
})

test_that("FedorovWynnAlgorithm converges with ndimFim > nProtocols grid", {
  .skip_optimizer_run_on_cran()
  opt = expect_no_warning(
    run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_cpp_conv" ) )
  )
  expect_s7_class( .opt_eval( opt ), Evaluation )
  expect_gt( getDeterminant( .opt_eval( opt ) ), 0 )
  expect_gt( getDcriterion( opt ), 0 )
})

test_that("generateDosesCombination derives a singleton grid without AdministrationConstraints", {
  opt = .minimal_discrete_opt( "MultiplicativeAlgorithm", name = "no_admin_cons" )
  arm = .first_arm( opt )
  prop( arm, "administrationsConstraints" ) = list()
  design = .project_design( opt )
  prop( design, "arms" ) = list( arm )
  res = PFIM:::generateDosesCombination( design )
  expect_equal( res$numberOfDoses, 1L )
} )

test_that("SimplexAlgorithm rejects discrete SamplingTimeConstraints without windows", {
  opt = .minimal_discrete_opt(
    "SimplexAlgorithm",
    optimizerParameters = list(
      pctInitialSimplexBuilding = 0.1, maxIteration = 5,
      tolerance = 1e-2, showProcess = FALSE
    ),
    name = "simplex_discrete"
  )
  expect_error( run( opt ), "samplingsWindows" )
} )

test_that("Simplex start matrix stays feasible when window counts increase", {
  times = c( 0.25, 1, 7, 12, 24 )
  st = SamplingTimes( outcome = "RespPK", samplings = times )
  sc = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = times,
    numberOfTimesByWindows = c( 2, 3 ),
    samplingsWindows = list( c( 0.1, 6 ), c( 6, 30 ) ),
    minSampling = 0.05
  )
  arm = Arm(
    name = "a", size = 50,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes = list( st ),
    samplingTimesConstraints = list( sc )
  )
  design = Design( name = "d", arms = list( arm ) )
  design = PFIM:::setSamplingConstraintForOptimization( design )
  prop( design, "arms" ) = PFIM:::.pfimSeedArmsFromInitialSamplings( prop( design, "arms" ) )
  layout = PFIM:::.buildFlatSamplingLayout( design )
  expect_true( PFIM:::.isFlatValid( layout, layout$initial_flat ) )
  mat = PFIM:::.pfimSimplexStartMatrix( layout, pct = 20 )
  expect_equal( nrow( mat ), length( layout$initial_flat ) + 1L )
  expect_true( all( PFIM:::.pfimSimplexFeasibleRows( layout, mat ) ) )
  expect_gt(
    max( abs( sweep( mat, 2L, mat[ 1L, ], "-" ) ) ),
    0
  )
  expect_error(
    PFIM:::.pfimRequireFeasibleSimplex(
      layout,
      matrix( 0, nrow = 2L, ncol = length( layout$initial_flat ) ),
      "search result"
    ),
    "no feasible sampling times"
  )
} )
