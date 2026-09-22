test_that( "FD gradient preserves nominal model parameters", {
  ev = Evaluation(
    name = "fd_mu_preserve",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.5, omega = 0.2 ) ),
      ModelParameter( name = "V", distribution = LogNormal( mu = 3.5, omega = 0.3 ) )
    ),
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs = list( RespPK = "C1" ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list(
          Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 6, 12 ) )
        )
      ) )
    ) ),
    fimType = "population"
  )
  model = rebuildEvalModel( ev, finiteDifference = TRUE )
  arm   = prop( prop( ev, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  mu_before = map_dbl(
    prop( model, "modelParameters" ),
    ~ prop( prop( .x, "distribution" ), "mu" )
  )

  evaluateModelGradient( model, arm )

  mu_after = map_dbl(
    prop( model, "modelParameters" ),
    ~ prop( prop( .x, "distribution" ), "mu" )
  )
  expect_equal( mu_after, mu_before )
} )

test_that( "evaluateDesign matches independent per-arm runs", {
  fx     = cas10_design()
  design = cas10_design_two_arms()
  shared_ev = Evaluation(
    name                    = "two_arm_shared",
    modelFromLibrary        = cas10_model_from_library(),
    modelParameters         = cas10_model_parameters(),
    modelCovariates         = cas10_covariates(),
    modelCovariatesEquation = "exponential",
    modelError              = fx$modelError,
    designs                 = list( design ),
    fimType                 = "individual",
    outputs                 = list( "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
  model    = rebuildEvalModel( shared_ev, finiteDifference = TRUE )
  fimProto = IndividualFim()
  shared   = PFIM:::evaluateDesign( design, model, fimProto )

  run_one_arm = function( arm ) {
    solo = Design( name = prop( design, "name" ), arms = list( arm ) )
    ev = Evaluation(
      name                    = "solo",
      modelFromLibrary        = cas10_model_from_library(),
      modelParameters         = cas10_model_parameters(),
      modelCovariates         = cas10_covariates(),
      modelCovariatesEquation = "exponential",
      modelError              = fx$modelError,
      designs                 = list( solo ),
      fimType                 = "individual",
      outputs                 = list( "RespPK" ),
      odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
    )
    m = rebuildEvalModel( ev, finiteDifference = TRUE )
    PFIM:::evaluateDesign( solo, m, IndividualFim() )
  }

  arms = prop( design, "arms" )
  purrr::walk( seq_along( arms ), function( i ) {
    M_shared = prop(
      prop( prop( shared, "evaluationArms" )[[ i ]], "evaluationFim" ),
      "fisherMatrix"
    )
    M_solo = prop(
      prop( prop( run_one_arm( arms[[ i ]] ), "evaluationArms" )[[ 1L ]], "evaluationFim" ),
      "fisherMatrix"
    )
    expect_equal( M_shared, M_solo, tolerance = 1e-8 )
  } )
} )

test_that( "checkValiditySamplingConstraint rejects infeasible windows", {
  arm = Arm(
    name = "a", size = 40,
    administrations = list(
      Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
    ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 3 ) )
    ),
    samplingTimesConstraints = list(
      SamplingTimeConstraints(
        outcome                = "RespPK",
        samplingsWindows       = list( c( 0, 0.5 ) ),
        numberOfTimesByWindows = c( 3L ),
        minSampling            = 1
      )
    )
  )
  design = Design( name = "bad", arms = list( arm ) )
  expect_error(
    checkValiditySamplingConstraint( design ),
    "sampling times constraint is not possible"
  )
} )

test_that( "generateSamplingTimesCombination validates combn inputs", {
  arm = Arm(
    name = "a", size = 40,
    administrations = list(
      Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
    ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 1, 2 ) )
    ),
    samplingTimesConstraints = list(
      SamplingTimeConstraints(
        outcome                      = "RespPK",
        initialSamplings             = c( 1, 2, 3, 4 ),
        fixedTimes                   = c( 1, 2 ),
        numberOfsamplingsOptimisable = 1
      )
    )
  )
  design = Design( name = "d", arms = list( arm ) )
  expect_error(
    generateSamplingTimesCombination( design ),
    "cannot be less than the number of fixed times"
  )
} )

test_that( "MultiplicativeAlgorithm_Rcpp reports singular FIM", {
  p   = 2L
  bad = matrix( c( 1, 1, 1, 1 ), 2, 2 )
  res = MultiplicativeAlgorithm_Rcpp(
    fisherMatrices = list( bad, bad ),
    n_fim = 2L,
    weights = c( 0.5, 0.5 ),
    p = p,
    lambda = 0.99,
    delta = 1e-3,
    iteration_init = 5L,
    show_process = FALSE
  )
  expect_true( isTRUE( res$singularFim ) )
  expect_false( isTRUE( res$converged ) )
} )

test_that( "continuous optimizers preserve distinct initial and optimal evaluations", {
  algos = c( "PSOAlgorithm", "PGBOAlgorithm", "SimplexAlgorithm" )
  purrr::walk( algos, function( algo ) {
    opt = suppressWarnings( run(
      .minimal_continuous_opt( algo, name = paste0( "init_opt_", algo ) )
    ) )
    od = prop( opt, "optimisationDesign" )
    expect_false( identical( od$evaluationInitialDesign, od$evaluationOptimalDesign ) )
    init_ev = od$evaluationInitialDesign
    opt_ev  = od$evaluationOptimalDesign
    init_arm = prop( prop( init_ev, "evaluationDesign" )[[ 1L ]], "arms" )[[ 1L ]]
    opt_arm  = prop( prop( opt_ev,  "evaluationDesign" )[[ 1L ]], "arms" )[[ 1L ]]
    expect_false( identical( init_arm, opt_arm ) )
  } )
} )

test_that( "pso_optimize_Rcpp rejects non-positive population size", {
  expect_error(
    PFIM:::pso_optimize_Rcpp(
      n_pop_in = 0L, max_iter = 1L, initial_pos = c( 1, 2 ),
      windows_list = list( matrix( c( 0, 10 ), 1 ), matrix( c( 0, 10 ), 1 ) ),
      sorting_groups = list( c( 1L, 2L ) ),
      phi1 = 2.05, phi2 = 2.05, constriction = 0.7, show_process = FALSE,
      eval_fitness = function( x ) sum( x^2 )
    ),
    "n_pop_in must be positive"
  )
} )

test_that( "prop<- writes Optimization project fields", {
  opt = .minimal_mult_opt()
  expect_equal( prop( opt, "name" ), "mult_opt_test" )
  prop( opt, "name" ) = "via_prop"
  expect_equal( prop( opt, "name" ), "via_prop" )
  prop( opt, "fimType" ) = "individual"
  expect_equal( prop( opt, "fimType" ), "individual" )
} )

test_that( "Administration outcome may be a state or a mapped outputs alias", {
  mk = function( admin_outcome, samp_outcome ) {
    Evaluation(
      name = "user_ode_admin",
      modelEquations = list( Deriv_Cc = "-k*Cc", Deriv_E = "Rin - kout*E" ),
      modelParameters = list(
        ModelParameter( name = "k", distribution = LogNormal( mu = 0.1, omega = 0.2 ) ),
        ModelParameter( name = "Rin", distribution = LogNormal( mu = 1, omega = 0.2 ) ),
        ModelParameter( name = "kout", distribution = LogNormal( mu = 0.2, omega = 0.2 ) )
      ),
      modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
      outputs = list( RespPK = "Cc", RespPD = "E" ),
      designs = list( Design(
        name = "d",
        arms = list( Arm(
          name = "a", size = 10,
          administrations = list(
            Administration( outcome = admin_outcome, timeDose = 0, dose = 1 )
          ),
          samplingTimes = list(
            SamplingTimes( outcome = samp_outcome, samplings = c( 1, 2 ) )
          ),
          initialConditions = list( Cc = 0, E = 0 )
        ) )
      ) ),
      fimType = "population"
    )
  }
  expect_silent( PFIM:::.pfimValidateProjectOutcomes( mk( "RespPK", "Cc" ) ) )
  expect_silent( PFIM:::.pfimValidateProjectOutcomes( mk( "Cc", "RespPK" ) ) )
  expect_silent( PFIM:::.pfimValidateProjectOutcomes( mk( "Cc", "Cc" ) ) )
  err = tryCatch(
    PFIM:::.pfimValidateProjectOutcomes( mk( "NotAState", "Cc" ) ),
    error = function( e ) conditionMessage( e )
  )
  expect_match( err, "^PFIM:.*Administration" )
  expect_false( grepl( "In index", err, fixed = TRUE ) )
} )

test_that( "ODE IC expressions reject session symbols", {
  env = PFIM:::.odeIcEvalEnv( c( V = 10 ), list( dose_Cc = 100 ) )
  expect_equal( PFIM:::.pfimEvalIcExpr( "dose_Cc/V", env, "Cc" ), 10 )
  expect_error(
    PFIM:::.pfimEvalIcExpr( "missing_symbol", env, "Cc" ),
    "unknown name"
  )
  expect_error(
    PFIM:::.pfimEvalIcExpr( "Sys.setenv(PFIM_IC_PROBE = '1')", env, "Cc" ),
    "unknown name"
  )
} )

test_that( ".pfimCloseOwnedGraphicsDevices leaves pre-existing devices", {
  f = tempfile( fileext = ".pdf" )
  grDevices::pdf( f )
  id = grDevices::dev.cur()
  on.exit( {
    if ( grDevices::dev.cur() > 1L )
      try( grDevices::dev.off(), silent = TRUE )
    unlink( f )
  } )
  open = grDevices::dev.list()
  PFIM:::.pfimCloseOwnedGraphicsDevices( open )
  expect_equal( grDevices::dev.cur(), id )
} )

test_that( "library 2-cpt ODE names keep catalogue states when only RespPK is mapped", {
  ev = Evaluation(
    name = "mm2_names",
    modelFromLibrary = list( PKModel = "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2" ),
    modelParameters = list(
      ModelParameter( name = "Vm",  distribution = LogNormal( mu = 0.08, omega = 0.1 ) ),
      ModelParameter( name = "Km",  distribution = LogNormal( mu = 0.40, omega = 0.1 ) ),
      ModelParameter( name = "V1",  distribution = LogNormal( mu = 10,   omega = 0.1 ) ),
      ModelParameter( name = "V2",  distribution = LogNormal( mu = 20,   omega = 0.1 ) ),
      ModelParameter( name = "k12", distribution = LogNormal( mu = 0.5,  omega = 0.1 ) ),
      ModelParameter( name = "k21", distribution = LogNormal( mu = 0.3,  omega = 0.1 ) )
    ),
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs = list( RespPK = "C1" ),
    designs = list( Design( name = "d", arms = list( Arm(
      name = "a", size = 10,
      administrations = list( Administration( outcome = "RespPK", timeDose = 0, Tinf = 1, dose = 100 ) ),
      samplingTimes = list( SamplingTimes( outcome = "RespPK", samplings = c( 1, 2 ) ) ),
      initialConditions = list( C1 = 0, C2 = 0 )
    ) ) ) ),
    fimType = "population"
  )
  nms = PFIM:::.derivativeNamesFromCompartments( ev, 2L, c( "Deriv_C1", "Deriv_C2" ) )
  expect_equal( nms, c( "Deriv_C1", "Deriv_C2" ) )
  expect_false( anyNA( nms ) )
} )

test_that( "categorical covariates reject unknown effects, negative proportions, unknown IOV labels", {
  expect_error(
    CategoricalCovariate(
      name = "Sex", categories = c( "M", "F" ),
      categoriesProportions = c( 1.2, -0.2 ),
      effects = list( F = c( V = log( 1.2 ) ) )
    ),
    "non-negative"
  )
  expect_error(
    CategoricalCovariate(
      name = "Sex", categories = c( "M", "F" ),
      categoriesProportions = c( 0.5, 0.5 ),
      effects = list( Femme = c( V = log( 1.2 ) ) )
    ),
    "not in categories"
  )
  expect_error(
    CategoricalCovariateWithIOV(
      name = "Trt", categories = c( "A", "B" ),
      sequences = list( c( "A", "C" ) ),
      sequencesProportions = 1,
      effects = list( B = c( V = log( 1.1 ) ) )
    ),
    "not in categories"
  )
} )

test_that( "passive time-dependent PD uses t_<own outcome>", {
  ev = Evaluation(
    name = "passive_t",
    modelEquations = list(
      RespPK = "dose_RespPK/V*exp(-k*t)",
      RespPD = "S0 + kprog*t"
    ),
    modelParameters = list(
      ModelParameter( name = "V",     distribution = LogNormal( mu = 10, omega = 0.2 ) ),
      ModelParameter( name = "k",     distribution = LogNormal( mu = 0.2, omega = 0.2 ) ),
      ModelParameter( name = "S0",    distribution = LogNormal( mu = 1, omega = 0.2 ) ),
      ModelParameter( name = "kprog", distribution = LogNormal( mu = 0.1, omega = 0.2 ) )
    ),
    modelError = list(
      Constant( output = "RespPK", sigmaInter = 0.1 ),
      Constant( output = "RespPD", sigmaInter = 0.1 )
    ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 20,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 4 ) ),
          SamplingTimes( outcome = "RespPD", samplings = c( 1, 2, 4 ) )
        )
      ) )
    ) ),
    fimType = "population"
  )
  ev = run( ev )
  expect_gt( getDeterminant( ev ), 0 )
} )

test_that( "ImmediateBaselineLinear catalogue evaluates with PK", {
  ev = Evaluation(
    name = "baseline_lin",
    modelFromLibrary = list(
      PKModel = "Linear1BolusSingleDose_kV",
      PDModel = "ImmediateBaselineLinear_S0kprog"
    ),
    modelParameters = list(
      ModelParameter( name = "k",     distribution = LogNormal( mu = 0.2, omega = 0.2 ) ),
      ModelParameter( name = "V",     distribution = LogNormal( mu = 10, omega = 0.2 ) ),
      ModelParameter( name = "S0",    distribution = LogNormal( mu = 1, omega = 0.2 ) ),
      ModelParameter( name = "kprog", distribution = LogNormal( mu = 0.1, omega = 0.2 ) )
    ),
    modelError = list(
      Constant( output = "RespPK", sigmaInter = 0.1 ),
      Constant( output = "RespPD", sigmaInter = 0.1 )
    ),
    outputs = list( RespPK = "RespPK", RespPD = "RespPD" ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 20,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 4 ) ),
          SamplingTimes( outcome = "RespPD", samplings = c( 1, 2, 4 ) )
        )
      ) )
    ) ),
    fimType = "population"
  )
  ev = run( ev )
  expect_gt( getDeterminant( ev ), 0 )
} )

test_that( "dose_RespPK without a matching administration is a PFIM error", {
  ev = Evaluation(
    name = "dose_mismatch",
    modelEquations = list( Deriv_C1 = "-Vm*C1/(Km+C1) + dose_RespPK/V*ka*exp(-ka*t)" ),
    modelParameters = list(
      ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.2 ) ),
      ModelParameter( name = "V",  distribution = LogNormal( mu = 10, omega = 0.2 ) ),
      ModelParameter( name = "Vm", distribution = LogNormal( mu = 1, omega = 0.2 ) ),
      ModelParameter( name = "Km", distribution = LogNormal( mu = 1, omega = 0.2 ) )
    ),
    modelError = list( Constant( output = "C1", sigmaInter = 0.1 ) ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 10,
        administrations = list( Administration( outcome = "C1", timeDose = 0, dose = 100 ) ),
        samplingTimes = list( SamplingTimes( outcome = "C1", samplings = c( 1, 2 ) ) ),
        initialConditions = list( C1 = 0 )
      ) )
    ) ),
    fimType = "population"
  )
  expect_error( run( ev ), "^PFIM:.*dose_" )
} )

test_that( "declared zero covariate effect stays in the population FIM", {
  expect_warning(
    Covariate(
      name = "Sex", categories = c( "M", "F" ),
      categoriesProportions = c( 0.5, 0.5 ),
      effects = list( F = c( V = 0 ) )
    ),
    "stay in the FIM"
  )
  sex = suppressWarnings( Covariate(
    name = "Sex", categories = c( "M", "F" ),
    categoriesProportions = c( 0.5, 0.5 ),
    effects = list( F = c( V = 0 ) )
  ) )
  ev = Evaluation(
    name = "zero_beta",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
      ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
    ),
    modelCovariates = list( sex ),
    modelCovariatesEquation = "exponential",
    modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs = list( "RespPK" ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list(
          Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 6, 12 ) )
        )
      ) )
    ) ),
    fimType = "population"
  )
  out = run( ev )
  rn = rownames( getFisherMatrix( out )$fisherMatrix )
  expect_true( any( grepl( "beta_V_Sex_F|\u03b2_V_Sex_F", rn ) ) )
} )

test_that( "continuous optimizer warns when initialSamplings differ from Arm samplingTimes", {
  st = SamplingTimes( outcome = "RespPK", samplings = c( 0.3, 1.5, 8, 14, 24 ) )
  sc = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = c( 0.5, 2, 8, 12, 24 ),
    numberOfTimesByWindows = c( 2, 3 ),
    samplingsWindows = list( c( 0.1, 6 ), c( 6, 30 ) ),
    minSampling = 0.05
  )
  arm = Arm(
    name = "scB", size = 50,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes = list( st ),
    samplingTimesConstraints = list( sc )
  )
  expect_warning(
    PFIM:::.pfimSeedArmsFromInitialSamplings( list( arm ) ),
    "initialSamplings"
  )
} )

test_that( "checkValiditySamplingConstraint rejects arm times outside windows", {
  sc = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = c( 0.5, 2, 8, 12, 24 ),
    numberOfTimesByWindows = c( 2, 3 ),
    samplingsWindows = list( c( 0.1, 6 ), c( 6, 30 ) ),
    minSampling = 0.05
  )
  arm = Arm(
    name = "a", size = 50,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.25, 1, 3, 12, 24 ) )
    ),
    samplingTimesConstraints = list( sc )
  )
  expect_error(
    checkValiditySamplingConstraint( Design( name = "d", arms = list( arm ) ) ),
    "do not place"
  )
} )
