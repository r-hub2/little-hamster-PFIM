# Shared builders for ODE evaluation and optimization tests.

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
    arm_args$initialCondition = initial_condition

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
    optimizerParameters = list(
      lambda             = 0.99,
      numberOfIterations = 30,
      weightThreshold    = 0.01,
      delta              = 1e-4
    ),
    designs             = list( design ),
    fimType             = "individual",
    fim                 = IndividualFim(),
    outputs             = list( "RespPK" )
  )
}

#' @keywords internal
.expect_optimization_run = function( opt ) {
  opt = run( opt )
  od  = prop( opt, "optimisationDesign" )
  expect_type( od, "list" )
  expect_s7_class( od$evaluationInitialDesign, Evaluation )
  expect_s7_class( od$evaluationOptimalDesign, Evaluation )
  det_init = getDeterminant( od$evaluationInitialDesign )
  det_opt  = getDeterminant( od$evaluationOptimalDesign )
  expect_true( is.finite( det_init ) && det_init > 0 )
  expect_true( is.finite( det_opt ) && det_opt > 0 )
  expect_gt( getDcriterion( opt ), 0 )
  invisible( opt )
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
      MultiplicativeAlgorithm = list(
        lambda = 0.99, numberOfIterations = 30,
        weightThreshold = 0.01, delta = 1e-4
      ),
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
        tolerance = 1e-2, seed = 42, showProcess = FALSE
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
