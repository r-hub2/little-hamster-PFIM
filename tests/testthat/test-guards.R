# Edge-case and constructor guards for PFIM 7.1.

.bayesian_fixed_omega_eval = function() {
  arm = Arm(
    name            = "bay_fix_plot",
    size            = 40,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 8 ) ) )
  )
  Evaluation(
    name             = "bay_fixed_omega_plot",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters  = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
      ModelParameter(
        name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ),
        fixedOmega = TRUE
      )
    ),
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 ) ),
    designs    = list( Design( name = "d", arms = list( arm ) ) ),
    fimType    = "Bayesian",
    outputs    = list( "RespPK" )
  )
}

test_that( "Bayesian SE/RSE plot aligns with fixed-omega FIM columns", {
  ev = run( .bayesian_fixed_omega_eval() )
  # fixedOmega keeps η_V in the Bayesian FIM (MAP still has that random effect).
  expect_equal( ncol( getFisherMatrix( ev )$fisherMatrix ), 2L )

  px = PFIM:::.bayesianSeRsePlotData( ev, "SE" )
  expect_equal( nrow( px$data ), 2L )
  expect_equal( px$data$Parameter, c( "k", "V" ) )
  expect_true( all( is.finite( px$data$SE ) ) )

  expect_length( tablesForReport( prop( ev, "fim" ), ev )$SEAndRSETable, 1L )
  shrink = prop( PFIM:::setEvaluationFim( prop( ev, "fim" ), ev ), "shrinkage" )
  expect_equal( ncol( shrink ), 2L )
} )

test_that( "unknown fimType errors in Evaluation constructor", {
  expect_error(
    Evaluation(
      name             = "bad_fim",
      modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
      modelParameters  = list(
        ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
        ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
      ),
      modelError = list( Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 ) ),
      designs = list( Design(
        name = "d",
        arms = list( Arm(
          name = "a", size = 10,
          administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
          samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 1, 2 ) ) )
        ) )
      ) ),
      fimType = "not_a_fim",
      outputs = list( "RespPK" )
    ),
    "Invalid fimType"
  )
} )

test_that( "kblAlign matches column counts", {
  expect_equal( PFIM:::.kblAlign( 8L, 3L ), c( "l", "l", "l", "c", "c", "c", "c", "c" ) )
  expect_equal( PFIM:::.kblAlign( 5L, 1L ), c( "l", "c", "c", "c", "c" ) )
  expect_equal( PFIM:::.kblAlign( 1L, 1L ), "l" )
} )

test_that( "C++ kernels reject empty / mismatched dimensions", {
  expect_error(
    computePopFimCombo_Rcpp(
      matrix( 0, 0, 2 ), numeric( 0 ), numeric( 0 ), numeric( 0 ),
      diag( 2 ), integer( 0 ), list(), FALSE
    ),
    "at least one parameter"
  )

  expect_error(
    computePopFimCombo_Rcpp(
      matrix( 1, 1, 2 ), 1, 0.25, 0,
      matrix( 1, 2, 3 ), 2L, list(), FALSE
    ),
    "square"
  )

  expect_error(
    computePopFimCombo_Rcpp(
      matrix( 1, 1, 3 ), 1, 0.25, 0.1,
      diag( 3 ), c( 2L, 0L, 1L ), list(), TRUE
    ),
    "positive"
  )

  p = matrix( c( 1, 2, 3, 4 ), nrow = 2 )
  expect_error(
    fun_amoeba_Rcpp( p, 1, 1e-4, 2L, function( x, ... ) 1, NULL, NULL, FALSE ),
    "length\\(y\\) must equal nrow\\(p\\)"
  )

  expect_error(
    fun_amoeba_Rcpp( p, c( 1, 2 ), 1e-4, -1L, function( x, ... ) 1, NULL, NULL, FALSE ),
    "itmax"
  )

  expect_error(
    pgbo_optimize_Rcpp(
      c( 1, 2 ), list( integer( 0 ) ), 1L, 1L, 0.1, 1L, 0.5, 5L, 0.5, FALSE,
      function( x ) 1, function( x, g ) TRUE
    ),
    "empty groups"
  )

  expect_error(
    pgbo_optimize_Rcpp(
      c( 1, 2 ), list( c( 1L, 99L ) ), 1L, 1L, 0.1, 1L, 0.5, 5L, 0.5, FALSE,
      function( x ) 1, function( x, g ) TRUE
    ),
    "1\\.\\.length\\(pos\\)"
  )

  expect_error(
    pso_optimize_Rcpp(
      2L, 1L, c( 1, 2 ),
      list( matrix( c( 0, 10 ), 1 ), matrix( c( 0, 10 ), 1 ) ),
      list( integer( 0 ) ),
      0.5, 0.5, 0.5, FALSE,
      function( x ) 1
    ),
    "empty groups"
  )

  expect_error(
    pso_optimize_Rcpp(
      2L, 1L, c( 1, 2 ),
      list( matrix( c( 0, 10 ), 1 ), matrix( c( 0, 10 ), 1 ) ),
      list( c( 1L, 99L ) ),
      0.5, 0.5, 0.5, FALSE,
      function( x ) 1
    ),
    "1\\.\\.n_dims"
  )

  expect_error(
    pso_optimize_Rcpp(
      2L, 1L, c( 1, 2 ),
      list( matrix( c( 0, 10 ), 1 ) ),
      list( c( 1L, 2L ) ),
      0.5, 0.5, 0.5, FALSE,
      function( x ) 1
    ),
    "windows_list length"
  )

  expect_error(
    pso_optimize_Rcpp(
      2L, 1L, c( 1, 2 ),
      list( matrix( 0, 1, 1 ), matrix( c( 0, 10 ), 1 ) ),
      list( c( 1L, 2L ) ),
      0.5, 0.5, 0.5, FALSE,
      function( x ) 1
    ),
    ">= 2 columns"
  )
} )

test_that( "fedorovWynnBufferSizes covers n_times for vectps", {
  buf = PFIM:::.fedorovWynnBufferSizes( 2L, 2L, 8L )
  expect_equal( buf$nBuf, 4L )
  expect_equal( buf$nVectps, 8L )
} )

.pfim_fw_minimal_call = function( protdep, freqdep, delta = 1e-4 ) {
  n_cand  = 2L
  n_times = 2L
  p       = 1L
  n_buf   = 4L
  sample_times = matrix( c( 1, 2, 3, 4 ), nrow = 2L, byrow = TRUE )
  # p=1: scalar FIMs — protocol 2 (5) strictly dominates protocol 1 (2).
  info = matrix( c( 2.0, 5.0 ), nrow = 2L, ncol = 1L )
  protocols = list( n_cand, n_times, p, 0L, sample_times, info )
  FedorovWynnAlgorithm_Rcpp(
    protocols    = protocols,
    ndimen       = as.integer( c( n_cand, p, 20L ) ),
    nbprot       = as.integer( n_cand ),
    numprot      = integer( n_buf ),
    freq         = numeric( n_buf ),
    nbdata       = integer( n_buf ),
    vectps       = numeric( n_times ),
    fisher       = numeric( 1L ),
    error        = integer( 1L ),
    protdep      = as.integer( protdep ),
    freqdep      = as.numeric( freqdep ),
    show_process = FALSE,
    delta        = delta
  )
}

test_that( "Fedorov rejects duplicate initial candidates and sizes optimal times by support", {
  expect_error(
    .pfim_fw_minimal_call( c( 2L, 1L, 1L ), c( 0.5, 0.5 ) ),
    "duplicate initial candidate"
  )

  out = .pfim_fw_minimal_call( c( 1L, 1L ), 1.0 )
  n_active = sum( out$numprot > 0L )
  expect_gt( n_active, 0L )
  expect_equal( nrow( out$optimal_sampling_times ), n_active )
  expect_equal( ncol( out$optimal_sampling_times ), 2L )
} )

test_that( "Fedorov p=1 replaces support when alpha-star equals 1", {
  # Seed on weak protocol 1; optimal is protocol 2 with weight 1.
  out = .pfim_fw_minimal_call( c( 1L, 1L ), 1.0 )
  expect_true( isTRUE( out$converged ) )
  expect_identical( out$status, "success" )
  active = which( out$numprot > 0L & out$freq > 0 )
  expect_equal( out$numprot[ active ], 2L )
  expect_equal( sum( out$freq[ active ] ), 1, tolerance = 1e-10 )
} )

test_that( "Fedorov runs Wynn before optimality check with poor initial weights", {
  # Both candidates start in support with bad weights (0.9 on weak, 0.1 on strong).
  # No outsider to add; Wynn must rebalance then meet the optimality test.
  out = .pfim_fw_minimal_call( c( 2L, 1L, 2L ), c( 0.9, 0.1 ) )
  expect_true( isTRUE( out$converged ) )
  expect_identical( out$status, "success" )
  active = which( out$numprot > 0L & out$freq > 1e-8 )
  expect_true( 2L %in% out$numprot[ active ] )
  w2 = out$freq[ out$numprot == 2L ]
  expect_gt( sum( w2 ), 0.99 )
} )

test_that( "Fedorov delta rejects negative values", {
  expect_error(
    .pfim_fw_minimal_call( c( 1L, 1L ), 1.0, delta = -1 ),
    "delta must be a finite non-negative"
  )
} )

test_that( "sensitivity x-lab helper uses plotmath caption", {
  lab = PFIM:::.pfimSensitivityXLab( "h", "design_1", "armA", "RespPK", "mu_ka" )
  expect_true( grepl( "atop\\(", lab, perl = TRUE ) )
  expect_true( grepl( "Design:", lab, fixed = TRUE ) )
  expect_type( PFIM:::.pfimParsePlotmath( lab ), "language" )
} )

test_that( "constructors reject NA tau, outcome, gamma, and sizes", {
  expect_error(
    Administration( outcome = "RespPK", timeDose = 0, dose = 1, tau = NA_real_ ),
    "tau must be a non-negative finite"
  )
  expect_error(
    Administration( outcome = NA_character_, timeDose = 0, dose = 1 ),
    "outcome must be a non-empty string"
  )
  expect_error(
    SamplingTimes( outcome = NA_character_, samplings = 1 ),
    "outcome must be a non-empty string"
  )
  expect_error(
    ModelParameter( name = "k", gamma = NA_real_, distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    "gamma must be non-negative"
  )
  expect_error(
    LogNormal( mu = 1, omega = NA_real_ ),
    "omega must be non-negative"
  )
  expect_error(
    Covariate( name = NA_character_ ),
    "name must be a non-empty string"
  )
  expect_error(
    Arm( name = "a", size = NA_real_ ),
    "size must be non-negative"
  )
  expect_error(
    Arm( name = NA_character_, size = 10 ),
    "name must be a non-empty string"
  )
  expect_error(
    ModelParameter( name = NA_character_, distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    "name must be a non-empty string"
  )
  expect_error(
    LogNormal( mu = NA_real_, omega = 0.3 ),
    "mu must be finite"
  )
  expect_error(
    Combined1( output = NA_character_, sigmaInter = 0.1, sigmaSlope = 0.1 ),
    "output must be a non-empty string"
  )
  expect_error(
    Combined1( output = "y", sigmaInter = NA_real_, sigmaSlope = 0.1 ),
    "sigmaInter and sigmaSlope must be non-negative"
  )
  expect_error(
    Combined1( output = "y", sigmaInter = 0.1, sigmaSlope = 0.1, cError = NA_real_ ),
    "cError must be positive"
  )
  expect_error(
    Constant( output = "y", sigmaInter = 0.1, sigmaSlope = NA_real_ ),
    "sigmaSlope must be 0"
  )
} )

test_that( "Arm initialCondition is an alias of initialConditions", {
  ic = list( C1 = 0, C2 = 1 )
  arm = Arm( name = "a", size = 10, initialCondition = ic )
  expect_equal( prop( arm, "initialConditions" ), ic )
  expect_error(
    Arm( name = "a", size = 10, initialConditions = ic, initialCondition = ic ),
    "not both"
  )
} )
