# ODE eval/opt test builders.
#' Skip full optimizer `run()` on CRAN (keep fast unit tests).
#' @keywords internal
.skip_optimizer_run_on_cran = function() skip_on_cran()

#' Skip Report() tests when rmarkdown or pandoc is unavailable.
#' @keywords internal
.skip_if_no_pandoc = function() {
  skip_if_not_installed( "rmarkdown" )
  skip_if_not( rmarkdown::pandoc_available(), "pandoc not available" )
}

#' @keywords internal
local_pfim_opts = function( opts, env = parent.frame() ) {
  old = setNames( lapply( names( opts ), pfim_get_option ), names( opts ) )
  do.call( pfim_set_option, opts )
  withr::defer( do.call( pfim_set_option, old ), envir = env )
}

#' @keywords internal
.ode_solver_params = function() {
  list( atol = 1e-8, rtol = 1e-8 )
}

#' @keywords internal
.expect_valid_fim = function( evaluation, min_dim = 2L ) {
  fim = getFisherMatrix( evaluation )
  M   = fim$fisherMatrix
  expect_true( is.matrix( M ) )
  expect_equal( nrow( M ), ncol( M ) )
  expect_gte( nrow( M ), min_dim )
  expect_true( all( is.finite( M ) ) )
  expect_gt( getDeterminant( evaluation ), 0 )
}

#' @keywords internal
.ode_evaluation = function( pk_model,
                              model_parameters,
                              administration,
                              sampling_times,
                              outputs,
                              fim_type = "population",
                              initial_condition = NULL,
                              name = pk_model ) {
  arm_args = list(
    name              = "ode_arm",
    size              = 200,
    administrations   = list( administration ),
    samplingTimes     = list( sampling_times )
  )
  if ( !is.null( initial_condition ) )
    arm_args$initialConditions = initial_condition

  arm    = do.call( Arm, arm_args )
  design = Design( name = "ode_design", arms = list( arm ) )

  Evaluation(
    name                    = name,
    modelFromLibrary        = list( PKModel = pk_model ),
    modelParameters         = model_parameters,
    modelError              = list(
      Combined1( output = names( outputs )[[ 1L ]], sigmaInter = 0.5, sigmaSlope = 0.15 )
    ),
    outputs                 = outputs,
    designs                 = list( design ),
    fimType                 = fim_type,
    odeSolverParameters     = .ode_solver_params()
  )
}

#' @keywords internal
.expect_ode_model_class = function( evaluation, pk_model, expected_class ) {
  raw_eq = prop( PFIM:::pkModelLibrary, "models" )[[ pk_model ]]
  prop( evaluation, "modelEquations" ) = raw_eq
  model = defineModelType( evaluation )
  expect_s7_class( model, expected_class )
  invisible( model )
}

#' Default multiplicative settings for fast discrete tests.
#' @keywords internal
.mult_test_algo_params = function() {
  list(
    lambda             = 0.99,
    numberOfIterations = 200L,
    weightThreshold    = 0.01,
    delta              = 1e-3
  )
}

#' @keywords internal
.minimal_mult_opt = function() {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  ac    = AdministrationConstraints( outcome = "RespPK", doses = as.list( 100 ) )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 3, 4 ) )
  sc    = SamplingTimeConstraints(
    outcome                      = "RespPK",
    initialSamplings             = c( 1, 2, 3, 4 ),
    fixedTimes                   = c( 1 ),
    numberOfsamplingsOptimisable = 2
  )
  arm = Arm(
    name                       = "opt_arm",
    size                       = 30,
    administrations            = list( admin ),
    administrationsConstraints = list( ac ),
    samplingTimes              = list( st ),
    samplingTimesConstraints   = list( sc )
  )
  design = Design( name = "opt_design", arms = list( arm ) )
  params = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
  )
  err = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )

  Optimization(
    name                = "mult_opt_test",
    modelFromLibrary    = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters     = params,
    modelError          = list( err ),
    optimizer           = "MultiplicativeAlgorithm",
    optimizerParameters = .mult_test_algo_params(),
    designs             = list( design ),
    fimType             = "individual",
    fim                 = IndividualFim(),
    outputs             = list( "RespPK" )
  )
}

#' @keywords internal
.expect_optimization_run = function( opt ) {
  opt = suppressWarnings( run( opt ) )
  expect_s7_class( .opt_init_eval( opt ), Evaluation )
  expect_s7_class( .opt_eval( opt ), Evaluation )
  det_init = getDeterminant( .opt_init_eval( opt ) )
  det_opt  = getDeterminant( .opt_eval( opt ) )
  expect_true( is.finite( det_init ) && det_init > 0 )
  expect_true( is.finite( det_opt ) && det_opt > 0 )
  expect_gt( getDcriterion( opt ), 0 )
  invisible( opt )
}

#' @keywords internal
.minimal_discrete_arm = function( name, size = 30 ) {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  ac    = AdministrationConstraints( outcome = "RespPK", doses = as.list( 100 ) )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 3, 4 ) )
  sc    = SamplingTimeConstraints(
    outcome                      = "RespPK",
    initialSamplings             = c( 1, 2, 3, 4 ),
    fixedTimes                   = c( 1 ),
    numberOfsamplingsOptimisable = 2
  )
  Arm(
    name                       = name,
    size                       = size,
    administrations            = list( admin ),
    administrationsConstraints = list( ac ),
    samplingTimes              = list( st ),
    samplingTimesConstraints   = list( sc )
  )
}

#' @keywords internal
.minimal_discrete_opt_multi_arm = function( optimizer = "MultiplicativeAlgorithm",
                                            n_arms    = 2L,
                                            fimType   = "individual",
                                            name      = "multi_arm_discrete" ) {
  arms   = map( paste0( "arm", seq_len( n_arms ) ), ~ .minimal_discrete_arm( .x ) )
  design = Design( name = "multi_arm_design", arms = arms )
  params = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
  )
  err = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
  if ( identical( optimizer, "FedorovWynnAlgorithm" ) ) {
    combos = generateSamplingTimesCombination( design )
    proto_row = function( idx ) {
      unlist( lapply( combos, function( arm_combos ) {
        i = min( idx, length( arm_combos ) )
        unlist( lapply( arm_combos[[ i ]], function( st ) prop( st, "samplings" ) ),
                use.names = FALSE )
      } ), use.names = FALSE )
    }
    optimizerParameters = list(
      elementaryProtocols = list( proto_row( 1L ), proto_row( 2L ) ),
      numberOfSubjects    = 30,
      proportionsOfSubjects = c( 0.5, 0.5 ),
      showProcess         = FALSE
    )
  } else {
    optimizerParameters = .mult_test_algo_params()
  }
  Optimization(
    name                = name,
    modelFromLibrary    = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters     = params,
    modelError          = list( err ),
    optimizer           = optimizer,
    optimizerParameters = optimizerParameters,
    designs             = list( design ),
    fimType             = fimType,
    fim                 = if ( identical( fimType, "population" ) ) PopulationFim() else IndividualFim(),
    outputs             = list( "RespPK" )
  )
}

#' Single-arm population multiplicative opt with two dose levels (constraint grid).
#' @keywords internal
.minimal_mult_pop_dose_grid = function( name = "mult_pop_dose_grid" ) {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  ac    = AdministrationConstraints( outcome = "RespPK", doses = list( c( 50 ), c( 100 ) ) )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 3, 4 ) )
  sc    = SamplingTimeConstraints(
    outcome                      = "RespPK",
    initialSamplings             = c( 1, 2, 3, 4 ),
    fixedTimes                   = c( 1 ),
    numberOfsamplingsOptimisable = 2
  )
  arm = Arm(
    name                       = "grid_arm",
    size                       = 100,
    administrations            = list( admin ),
    administrationsConstraints = list( ac ),
    samplingTimes              = list( st ),
    samplingTimesConstraints   = list( sc )
  )
  design = Design( name = "dose_grid", arms = list( arm ) )
  params = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
  )
  err = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
  Optimization(
    name                = name,
    modelFromLibrary    = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters     = params,
    modelError          = list( err ),
    optimizer           = "MultiplicativeAlgorithm",
    optimizerParameters = .mult_test_algo_params(),
    designs             = list( design ),
    fimType             = "population",
    fim                 = PopulationFim(),
    outputs             = list( "RespPK" )
  )
}

#' @keywords internal
.minimal_discrete_opt = function( optimizer,
                                   optimizerParameters = NULL,
                                   name = "discrete_opt_test" ) {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  ac    = AdministrationConstraints( outcome = "RespPK", doses = as.list( 100 ) )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 3, 4 ) )
  sc    = SamplingTimeConstraints(
    outcome                      = "RespPK",
    initialSamplings             = c( 1, 2, 3, 4 ),
    fixedTimes                   = c( 1 ),
    numberOfsamplingsOptimisable = 2
  )
  arm = Arm(
    name                       = "opt_arm",
    size                       = 30,
    administrations            = list( admin ),
    administrationsConstraints = list( ac ),
    samplingTimes              = list( st ),
    samplingTimesConstraints   = list( sc )
  )
  design = Design( name = "opt_design", arms = list( arm ) )
  params = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
  )
  err = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )

  if ( is.null( optimizerParameters ) ) {
    optimizerParameters = switch( optimizer,
      MultiplicativeAlgorithm = .mult_test_algo_params(),
      FedorovWynnAlgorithm = {
        combos = generateSamplingTimesCombination( design )[[ 1L ]]
        protos = map( combos[ seq_len( min( 2L, length( combos ) ) ) ], function( st_list ) {
          sort( unlist( map( st_list, ~ prop( .x, "samplings" ) ) ) )
        } )
        list(
          elementaryProtocols     = protos,
          numberOfSubjects        = 30,
          proportionsOfSubjects   = rep( 1 / length( protos ), length( protos ) ),
          showProcess             = FALSE
        )
      },
      stop( "Unknown discrete optimizer: ", optimizer, call. = FALSE )
    )
  }

  Optimization(
    name                = name,
    modelFromLibrary    = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters     = params,
    modelError          = list( err ),
    optimizer           = optimizer,
    optimizerParameters = optimizerParameters,
    designs             = list( design ),
    fimType             = "individual",
    fim                 = IndividualFim(),
    outputs             = list( "RespPK" )
  )
}

#' @keywords internal
.minimal_continuous_opt = function( optimizer, optimizerParameters = NULL, name = "continuous_opt_test" ) {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 3, 4 ) )
  sc    = SamplingTimeConstraints(
    outcome                = "RespPK",
    initialSamplings       = c( 1, 2, 3, 4 ),
    numberOfTimesByWindows = c( 2, 2 ),
    samplingsWindows       = list( c( 1, 2 ), c( 3, 4 ) ),
    minSampling            = 0.5
  )
  arm = Arm(
    name                     = "opt_arm",
    size                     = 30,
    administrations          = list( admin ),
    samplingTimes            = list( st ),
    samplingTimesConstraints = list( sc )
  )
  design = Design( name = "opt_design", arms = list( arm ) )
  params = list(
    ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
    ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
  )
  err = Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )

  if ( is.null( optimizerParameters ) ) {
    optimizerParameters = switch( optimizer,
      PSOAlgorithm = list(
        maxIteration = 5, populationSize = 5, seed = 42,
        personalLearningCoefficient = 2.05,
        globalLearningCoefficient = 2.05,
        showProcess = FALSE
      ),
      PGBOAlgorithm = list(
        N = 10, muteEffect = 0.1, maxIteration = 5,
        purgeIteration = 2, seed = 42, showProcess = FALSE
      ),
      SimplexAlgorithm = list(
        pctInitialSimplexBuilding = 0.1, maxIteration = 15,
        tolerance = 1e-2, showProcess = FALSE
      ),
      stop( "Unknown continuous optimizer: ", optimizer, call. = FALSE )
    )
  }

  Optimization(
    name                = name,
    modelFromLibrary    = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters     = params,
    modelError          = list( err ),
    optimizer           = optimizer,
    optimizerParameters = optimizerParameters,
    designs             = list( design ),
    fimType             = "individual",
    fim                 = IndividualFim(),
    outputs             = list( "RespPK" )
  )
}
