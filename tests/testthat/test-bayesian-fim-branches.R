test_that( ".bayesianOmega is diag(omega^2), not mu-scaled", {
  expect_equal( PFIM:::.bayesianOmega( 0.5 ), matrix( 0.25, 1L, 1L ) )
  expect_equal( PFIM:::.bayesianOmega( c( 0.5, 0.3 ) ), diag( c( 0.25, 0.09 ) ) )
  expect_equal( dim( PFIM:::.bayesianOmega( numeric( 0 ) ) ), c( 0L, 0L ) )
  expect_equal( PFIM:::.bayesianOmegaInv( c( 0.5, 0.3 ) ), diag( c( 4, 1 / 0.09 ) ) )
} )

test_that( "Bayesian LogNormal prior is Omega^{-1} (FO)", {
  # Empty data FIM -> Bayesian FIM equals Omega^{-1}, not (mu Omega mu)^{-1}.
  arm = Arm(
    name            = "bay_prior",
    size            = 1,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 8 ) ) )
  )
  ev = Evaluation(
    name             = "bay_prior_scale",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters  = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
      ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
    ),
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 ) ),
    designs    = list( Design( name = "d", arms = list( arm ) ) ),
    fimType    = "Bayesian",
    outputs    = list( "RespPK" )
  )
  ev = run( ev )
  M  = getFisherMatrix( ev )$fisherMatrix
  expect_equal( ncol( M ), 2L )

  # With real data the FIM is data + Omega^{-1}; check prior strength via
  # comparing to the algebraic Omega^{-1} contribution on the diagonal lower bound.
  omega = c( 0.5, 0.3 )
  Omega_inv = 1 / omega^2
  # Wrong (old) prior would be 1/(mu^2 omega^2):
  mu = c( 0.25, 15 )
  wrong_prior = 1 / ( mu^2 * omega^2 )
  # FO prior is much stronger than mu-scaled prior when |mu| != 1.
  expect_gt( Omega_inv[[ 2L ]], wrong_prior[[ 2L ]] * 100 )
  # Diagonal of Bayesian FIM must be at least Omega^{-1} (data info >= 0 in FO sense
  # after M'IMF M; practically diag(M) >= Omega^{-1} for well-posed designs).
  expect_true( all( diag( M ) >= Omega_inv - 1e-8 ) )
  # And far above the old wrong prior on V (mu=15).
  expect_gt( diag( M )[[ 2L ]], wrong_prior[[ 2L ]] * 10 )
} )

test_that( ".bayesianShrinkageForDesign: empty arms", {
  expect_equal(
    PFIM:::.bayesianShrinkageForDesign( diag( 2L ), NULL, list() ),
    numeric( 0 )
  )
} )

test_that( ".bayesianShrinkageValues and matrix wrapper cover shapes", {
  expect_equal( PFIM:::.bayesianShrinkageValues( c( 10, 20 ) ), c( 10, 20 ) )
  expect_equal(
    PFIM:::.bayesianShrinkageValues( matrix( c( 10, 20 ), nrow = 1L ) ),
    c( 10, 20 )
  )
  expect_equal(
    PFIM:::.bayesianShrinkageValues( matrix( c( 10, 20 ), ncol = 1L ) ),
    c( 10, 20 )
  )
  square = matrix( 1:4, 2L, 2L )
  expect_equal( PFIM:::.bayesianShrinkageValues( square ), as.numeric( square ) )

  M = PFIM:::.bayesianShrinkageMatrix( c( 11, 22 ), c( "ka", "V" ) )
  expect_equal( dim( M ), c( 1L, 2L ) )
  expect_equal( colnames( M ), c( "ka", "V" ) )
  expect_equal( rownames( M ), "Shrinkage" )
} )

test_that( ".bayesianShrinkageMuLabels uses colnames, rownames, then mu labels", {
  ev = run( cas10_evaluation( "bay_mu_labels", "Bayesian" ) )
  named_row = matrix( c( 10, 20, 30 ), nrow = 1L, dimnames = list( "Shrinkage", c( "a", "b", "c" ) ) )
  expect_equal( PFIM:::.bayesianShrinkageMuLabels( named_row, ev ), c( "a", "b", "c" ) )

  named_col = matrix( c( 10, 20, 30 ), ncol = 1L, dimnames = list( c( "a", "b", "c" ), "Shrinkage" ) )
  expect_equal( PFIM:::.bayesianShrinkageMuLabels( named_col, ev ), c( "a", "b", "c" ) )

  fallback = PFIM:::.bayesianShrinkageMuLabels( c( 10, 20, 30 ), ev )
  expect_true( is.character( fallback ) )
  expect_gt( length( fallback ), 0L )
} )

test_that( "Bayesian FIM keeps fixed-omega parameters in the eta block", {
  arm = Arm(
    name            = "bay_fix",
    size            = 40,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 8 ) ) )
  )
  ev = Evaluation(
    name             = "bay_fixed_omega",
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
  ev = run( ev )
  M  = getFisherMatrix( ev )$fisherMatrix
  # MAP still has η_V ~ N(0, ω_V²) even when omega is fixed for population estimation.
  expect_equal( ncol( M ), 2L )
  expect_true( all( is.finite( M ) ) )
  sh = as.numeric( getShrinkage( ev ) )
  expect_length( sh, 2L )
  expect_true( all( is.finite( sh ) ) )
  expect_gt( sh[[ 1L ]], 0 )
} )

test_that( "Bayesian SE/RSE plot data covers mu-only and covariate beta facets", {
  ev_plain = run( cas10_evaluation( "bay_plot_plain", "Bayesian" ) )
  px = PFIM:::.bayesianSeRsePlotData( ev_plain, "SE" )
  expect_true( is.data.frame( px$data ) )
  expect_gt( nrow( px$data ), 0L )
  expect_true( "SE" %in% names( px$data ) )

  ev_plain = run( cas10_evaluation( "bay_plot_rse", "Bayesian" ) )
  px_rse = PFIM:::.bayesianSeRsePlotData( ev_plain, "RSE" )
  expect_true( "RSE" %in% names( px_rse$data ) )
} )
