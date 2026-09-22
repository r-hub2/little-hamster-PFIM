# Integration reference cases (Evaluation_tests / Optimisation_tests).

eval_opt_reference_root = function() {
  env = Sys.getenv( "PFIM_TESTS_REF_ROOT", unset = "" )
  if ( nzchar( env ) && dir.exists( env ) )
    return( normalizePath( env, winslash = "/", mustWork = TRUE ) )

  pkg = Sys.getenv( "PFIM_PACKAGE_ROOT", unset = system.file( package = "PFIM" ) )
  candidates = c(
    normalizePath( file.path( pkg, "..", "..", "..", "tests_PFIM", "resultats_de_references" ),
                   winslash = "/", mustWork = FALSE ),
    normalizePath( file.path( pkg, "..", "..", "tests_PFIM", "resultats_de_references" ),
                   winslash = "/", mustWork = FALSE )
  )
  hits = candidates[ dir.exists( candidates ) ]
  if ( length( hits ) ) hits[[ 1L ]] else NA_character_
}

fim_summary = function( fim, d ) {
  M = prop( fim, "fisherMatrix" )
  list(
    d            = as.numeric( d ),
    det          = det( M ),
    trace        = sum( diag( M ) ),
    ncol         = ncol( M ),
    m11          = M[ 1L, 1L ],
    m_nn         = M[ ncol( M ), ncol( M ) ],
    fisherMatrix = M
  )
}

fim_summary_from_evaluation = function( evaluation ) {
  fim = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  fim_summary( fim, PFIM:::Dcriterion( fim ) )
}

fim_summary_from_optimization = function( optimization ) {
  fim = PFIM:::.getOptimalFim( optimization )
  fim_summary( fim, getDcriterion( optimization ) )
}

parse_report_fim_summaries = function( html_path ) {
  if ( !file.exists( html_path ) )
    return( list() )

  txt   = paste( readLines( html_path, warn = FALSE ), collapse = "\n" )
  parts = strsplit( txt, "Determinant, D-criterion", fixed = TRUE )[[ 1L ]][ -1L ]
  lapply( parts, function( part ) {
    seg = sub( "(?s)^.*?<tbody>(.*?)</tbody>.*", "\\1", part, perl = TRUE )
    if ( identical( seg, part ) )
      return( NULL )

    nums = suppressWarnings( as.numeric( unlist( regmatches(
      seg,
      gregexpr( "[0-9]+\\.[0-9]+(?:e[+-]?[0-9]+)?", seg, ignore.case = TRUE, perl = TRUE )
    ) ) ) )
    nums = nums[ is.finite( nums ) ]
    if ( length( nums ) < 2L )
      return( NULL )

    list( det = nums[ 1L ], d = nums[ 2L ] )
  } )
}

build_pk_analytic_admin1_evaluation = function( fim_type = "population" ) {
  modelFromLibrary = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" )
  modelParameters = list(
    ModelParameter( name = "V",  distribution = LogNormal( mu = 8,    omega = sqrt( 0.020 ) ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 0.13,  omega = sqrt( 0.06 ) ) ),
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1.6,   omega = sqrt( 0.7 ) ) )
  )
  modelError = list( Combined1( output = "RespPK", sigmaInter = 0.6, sigmaSlope = 0.07 ) )
  administration = Administration( outcome = "RespPK", timeDose = c( 0, 20 ), dose = c( 100, 100 ) )
  samplingTimes = SamplingTimes(
    outcome = "RespPK",
    samplings = c( 0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120 )
  )
  arm = Arm(
    name = "BrasTest", size = 32,
    administrations = list( administration ),
    samplingTimes   = list( samplingTimes )
  )
  Evaluation(
    name = "pk_analytic", modelFromLibrary = modelFromLibrary,
    modelParameters = modelParameters, modelError = modelError,
    outputs = list( "RespPK" ), designs = list( Design( name = "design1", arms = list( arm ) ) ),
    fimType = fim_type, odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )
}

build_pk_infusion_ss_pop_evaluation = function() {
  modelFromLibrary = list( PKModel = "Linear1InfusionSteadyState_ClVtau" )
  modelParameters = list(
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.09 ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = 0.09 ) )
  )
  administration = Administration( outcome = "RespPK", Tinf = 5, tau = 5, dose = 20 )
  samplingTimes = SamplingTimes(
    outcome = "RespPK",
    samplings = c( 0, 1, 2, 5, 7, 8, 10, 12, 14, 15, 16, 20, 21, 30, 40, 50, 60, 70, 80, 100 )
  )
  arm = Arm(
    name = "BrasTest", size = 40,
    administrations = list( administration ),
    samplingTimes   = list( samplingTimes )
  )
  Evaluation(
    name = "", modelFromLibrary = modelFromLibrary, modelParameters = modelParameters,
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.1, sigmaSlope = 0.1 ) ),
    outputs = list( "RespPK" ), designs = list( Design( name = "design1", arms = list( arm ) ) ),
    fimType = "population", odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )
}

build_mult_pk_ode_pop_optimization = function() {
  administration = Administration( outcome = "C1", timeDose = 0, dose = 50 )
  samplingTimes  = SamplingTimes( outcome = "C1", samplings = c( 0.5, 1, 24, 36, 120 ) )
  arm = Arm(
    name = "BrasTest1", size = 100,
    administrations = list( administration ),
    samplingTimes   = list( samplingTimes ),
    administrationsConstraints = list(
      AdministrationConstraints( outcome = "C1", doses = list( 50, 100, 150 ) )
    ),
    samplingTimesConstraints = list(
      SamplingTimeConstraints(
        outcome = "C1",
        initialSamplings = c( 0.5, 1, 2, 6, 96, 9, 12, 24, 36, 48, 72, 120 ),
        numberOfsamplingsOptimisable = 5,
        fixedTimes = c( 0.5, 96, 72 )
      )
    ),
    initialCondition = list( C1 = "dose_C1" )
  )
  Optimization(
    name = "PKPD_ODE_multi_doses_populationFIM",
    modelEquations = list( Deriv_C1 = "-k*C1" ),
    modelParameters = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.082, omega = sqrt( 0.25 ) ) )
    ),
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.6, sigmaSlope = 0.07 ) ),
    optimizer = "MultiplicativeAlgorithm",
    optimizerParameters = list(
      lambda = 0.99, numberOfIterations = 1000,
      weightThreshold = 0.001, delta = 1e-4, showProcess = FALSE
    ),
    designs = list( Design( name = "design1", arms = list( arm ), numberOfArms = 100 ) ),
    fimType = "population",
    outputs = list( RespPK = "C1" ),
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )
}

eval_opt_reference_cases = function() {
  ref_root = eval_opt_reference_root()
  mult_html = if ( !is.na( ref_root ) ) {
    file.path(
      ref_root,
      "Optimisation_tests/discrete/MultiplicativeAlgorithm",
      "multiplicative_Algorithm_PK_ode_dose_not_in_eqs_results/popFIM.html"
    )
  } else {
    NA_character_
  }

  list(
    list(
      id = "eval_pk_analytic_population_admin1",
      run = function() run( build_pk_analytic_admin1_evaluation( "population" ) ),
      summarize = fim_summary_from_evaluation,
      reference = list(
        d = 688.3841322, det = 5.042510457e+22, trace = 100772.7156,
        ncol = 8L, m11 = 21.67653169, m_nn = 36976.73506
      ),
      tol = list( d = 1e-6, det = 1e-6, trace = 0.05, m11 = 1e-4, m_nn = 0.05 )
    ),
    list(
      id = "eval_pk_analytic_individual_admin1",
      run = function() run( build_pk_analytic_admin1_evaluation( "individual" ) ),
      summarize = fim_summary_from_evaluation,
      reference = list(
        d = 75.21689308, det = 2407559074, trace = 15544.06304,
        ncol = 5L, m11 = 6.832029551, m_nn = 1483.888392
      ),
      tol = list( d = 1e-6, det = 1e-4, trace = 0.05, m11 = 1e-4, m_nn = 0.05 )
    ),
    list(
      id = "eval_pk_analytic_Bayesian_admin1",
      run = function() run( build_pk_analytic_admin1_evaluation( "Bayesian" ) ),
      summarize = fim_summary_from_evaluation,
      reference = list(
        # FO Bayesian prior Omega^{-1}.
        d = 140.7949664, det = 2791009.956, trace = 770.2487651,
        ncol = 3L, m11 = 487.2498913, m_nn = 29.26213695
      ),
      tol = list( d = 1e-6, det = 1e-4, trace = 0.05, m11 = 1e-4, m_nn = 1e-4 )
    ),
    list(
      id = "eval_pk_infusion_ss_pop_admin1",
      run = function() run( build_pk_infusion_ss_pop_evaluation() ),
      summarize = fim_summary_from_evaluation,
      reference = list(
        d = 3121.345726, det = 9.248072958e+20, trace = 353651.4695,
        ncol = 6L, m11 = 32.08687132, m_nn = 87799.86237
      ),
      tol = list( d = 1e-6, det = 1e-6, trace = 0.1, m11 = 1e-4, m_nn = 0.1 )
    ),
    list(
      id = "opt_mult_pk_ode_pop_initial",
      slow = TRUE,
      expect_converge_warning = TRUE,
      run = function() {
        opt = run( build_mult_pk_ode_pop_optimization() )
        prop( opt, "optimisationDesign" )$evaluationInitialDesign
      },
      summarize = fim_summary_from_evaluation,
      reference = list(
        d = 6810.422, det = 2.151276e+15, trace = NA_real_,
        ncol = 4L, m11 = NA_real_, m_nn = NA_real_
      ),
      tol = list( d = 1e-3, det = 1e-6 ),
      ref_html = mult_html, html_index = 1L
    ),
    list(
      id = "opt_mult_pk_ode_pop_optimal",
      slow = TRUE,
      expect_converge_warning = TRUE,
      run = function() run( build_mult_pk_ode_pop_optimization() ),
      summarize = fim_summary_from_optimization,
      reference = list(
        d = 7671.727, det = 3.463957e+15, trace = 130446.2,
        ncol = 4L, m11 = 58929.43, m_nn = 69647.11
      ),
      # Integer Hamilton allocation shifts D by ~1e-4 relative; keep loose tol.
      tol = list( d = 5e-3, det = 1e-3, trace = 1.0, m11 = 1.0, m_nn = 1.0 )
    )
  )
}
