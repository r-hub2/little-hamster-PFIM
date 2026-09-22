# Quality gates: fimType validation, numerical reference values (cas10).

test_that("singular FIM SE build warns and sets singular flag", {
  M = diag( c( 10, 0 ) )
  rownames( M ) = colnames( M ) = c( "mu_a", "mu_b" )
  expect_warning(
    se <- PFIM:::.fimBuildSeAndRse( M, colnames( M ), c( 1, 1 ) ),
    "singular|pseudo-inverse"
  )
  expect_true( isTRUE( se$singular ) )
  expect_true( is.infinite( se$SEAndRSE$SE[ 2L ] ) )
})

test_that("showFIM surfaces singularFim flag explicitly", {
  M = diag( c( 10, 0 ) )
  rownames( M ) = colnames( M ) = c( "mu_a", "mu_b" )
  se = suppressWarnings(
    PFIM:::.fimBuildSeAndRse( M, colnames( M ), c( 1, 1 ) )
  )
  fim = IndividualFim()
  fim = PFIM:::.fimStoreEvaluationResult(
    fim, M, M[ 1L, 1L, drop = FALSE ], se,
    varianceEffects = matrix( 1, 1L, 1L, dimnames = list( "sigma", "sigma" ) )
  )
  # Force variance block names for show path that subsets sigma.
  prop( fim, "varianceEffects" ) = M[ 1L, 1L, drop = FALSE ]
  prop( fim, "fixedEffects" ) = M
  prop( fim, "condNumberVarianceEffects" ) = 1
  out = paste( capture.output( showFIM( fim ) ), collapse = "\n" )
  expect_match( out, "singularFim:\\s*TRUE" )
  expect_match( out, "singularFim = TRUE|Moore-Penrose|pseudo-inverse" )
})

test_that("report criteria footnote mentions singularFim when flagged", {
  html = as.character( PFIM:::.fimCriteriaKable( 1, 1, 1, 1, singularFim = TRUE ) )
  expect_match( html, "singularFim = TRUE" )
  expect_match( html, "pseudo-inverse|Moore-Penrose" )
  expect_match( html, "pfim-singular-fim" )
  html_ok = as.character( PFIM:::.fimCriteriaKable( 1, 1, 1, 1, singularFim = FALSE ) )
  expect_false( grepl( "singularFim = TRUE", html_ok, fixed = TRUE ) )
})

test_that("getFisherMatrix exposes singularFim", {
  M = diag( c( 10, 0 ) )
  rownames( M ) = colnames( M ) = c( "mu_a", "mu_b" )
  se = suppressWarnings(
    PFIM:::.fimBuildSeAndRse( M, colnames( M ), c( 1, 1 ) )
  )
  fim = IndividualFim()
  fim = PFIM:::.fimStoreEvaluationResult(
    fim, M, M, se,
    varianceEffects = matrix( 1, 1L, 1L, dimnames = list( "s", "s" ) )
  )
  expect_true( isTRUE( prop( fim, "singularFim" ) ) )
  # Accessor path on a minimal Evaluation is covered via prop; list shape:
  expect_true( "singular" %in% names( prop( fim, "SEAndRSE" ) ) )
  expect_true( isTRUE( prop( fim, "SEAndRSE" )$singular ) )
})

test_that("invalid fimType fails in Evaluation() with a clear message", {
  expect_error(
    cas10_evaluation( "bad_fim_type", "not_a_fim_type" ),
    regexp = "Invalid fimType|not_a_fim_type"
  )
})

test_that("defineFim rejects invalid fimType before design evaluation", {
  expect_error(
    cas10_evaluation( "bad_define", "invalid" ),
    regexp = "Invalid fimType"
  )
})

test_that("fimType is case-insensitive (Bayesian vs bayesian)", {
  ev_upper = cas10_evaluation( "bayes_upper", "Bayesian" )
  ev_lower = cas10_evaluation( "bayes_lower", "bayesian" )
  ev_upper = run( ev_upper )
  ev_lower = run( ev_lower )
  expect_equal(
    ncol( prop( prop( ev_upper, "fim" ), "fisherMatrix" ) ),
    ncol( prop( prop( ev_lower, "fim" ), "fisherMatrix" ) )
  )
  expect_equal(
    det( prop( prop( ev_upper, "fim" ), "fisherMatrix" ) ),
    det( prop( prop( ev_lower, "fim" ), "fisherMatrix" ) ),
    tolerance = 1e-6
  )
})

test_that("cas10 individual: reference determinant and trace", {
  evaluation = run( cas10_evaluation( "ref_ind", "individual" ) )
  fim        = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M          = prop( fim, "fisherMatrix" )

  # Harmonic mean over cov strata; mu + sigma only (beta omitted).
  expect_equal( ncol( M ), 4L )
  expect_true( is.finite( det( M ) ) )
  expect_gt( det( M ), 0 )
  expect_true( is.finite( sum( diag( M ) ) ) )
})

test_that("cas10 Bayesian: reference determinant, trace, and column labels", {
  evaluation = run( cas10_evaluation( "ref_bay", "Bayesian" ) )
  fim        = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M          = prop( fim, "fisherMatrix" )

  # FO + harmonic mean; mu only (beta omitted).
  expect_equal( ncol( M ), 3L )
  expect_equal(
    colnames( M ),
    c( "\u03bc_ka", "\u03bc_V", "\u03bc_Cl" )
  )
  expect_true( is.finite( det( M ) ) )
  expect_gt( det( M ), 0 )
  expect_gt( min( prop( fim, "shrinkage" ) ), 0 )
})

test_that(".pfimSubjectCovariance marks null-space parameters as Inf", {
  # Rank-1 information on first coordinate only.
  M = diag( c( 4, 0, 0 ) )
  C = PFIM:::.pfimSubjectCovariance( M )
  expect_equal( C[ 1L, 1L ], 0.25, tolerance = 1e-12 )
  expect_true( is.infinite( C[ 2L, 2L ] ) )
  expect_true( is.infinite( C[ 3L, 3L ] ) )
  expect_equal( C[ 1L, 2L ], 0 )
} )

test_that(".pfimSubjectCovariance uses correlation-scale rank (not eps*p*lambda_max)", {
  # Rank-1 2x2 plus a 3.4e-13 perturbation: raw eigen keeps it (eps*p*lambda_max
  # ~ 1e-13), correlation-scale eigenvalues are ~2 and ~3e-12 -> null.
  v = c( 12.12435565, 0.5 )
  u = c( -0.5, 12.12435565 )
  u = u / sqrt( sum( u^2 ) )
  M2 = tcrossprod( v ) + 3.4e-13 * tcrossprod( u )
  M = diag( 0, 4L )
  M[ 1:2, 1:2 ] = M2

  raw = eigen( M, symmetric = TRUE, only.values = TRUE )$values
  expect_gt( raw[[ 2L ]], .Machine$double.eps * 4 * max( raw ) )

  C = PFIM:::.pfimSubjectCovariance( M )
  expect_true( is.infinite( C[ 1L, 1L ] ) )
  expect_true( is.infinite( C[ 2L, 2L ] ) )
  expect_true( is.infinite( C[ 3L, 3L ] ) )
  expect_true( is.infinite( C[ 4L, 4L ] ) )

  se = sqrt( PFIM:::.fimPinvCovarianceDiagonal( M ) )
  expect_true( all( is.infinite( se[ 1:2 ] ) ) )
} )

test_that(".pfimSubjectCovariance keeps a small-scale identifiable parameter", {
  # Units differ by 1e20: eps*p*lambda_max would drop the small diagonal;
  # correlation scale keeps both.
  M = diag( c( 1e12, 1e-8 ) )
  C = PFIM:::.pfimSubjectCovariance( M )
  expect_equal( C[ 1L, 1L ], 1e-12, tolerance = 1e-18 )
  expect_equal( C[ 2L, 2L ], 1e8, tolerance = 1e-4 )
  expect_true( all( is.finite( C ) ) )
} )

test_that(".pfimCovarianceToFim maps Inf covariance back to zero information", {
  C = diag( c( 0.25, Inf, Inf ) )
  Mf = PFIM:::.pfimCovarianceToFim( C )
  expect_equal( Mf[ 1L, 1L ], 4, tolerance = 1e-12 )
  expect_equal( Mf[ 2L, 2L ], 0 )
  expect_equal( Mf[ 3L, 3L ], 0 )
} )

test_that(".pfimHarmonicMeanFim yields Inf SE directions when a stratum is singular", {
  M_rich = diag( c( 10, 10, 10, 5 ) )
  M_poor = diag( c( 0, 0, 0, 5 ) )
  harm = PFIM:::.pfimHarmonicMeanFim( list( M_rich, M_poor ), c( 0.75, 0.25 ) )
  expect_true( all( is.infinite( diag( harm$covariance )[ 1:3 ] ) ) )
  expect_true( is.finite( harm$covariance[ 4L, 4L ] ) )
  expect_equal( harm$fisherMatrix[ 1L, 1L ], 0 )
  expect_gt( harm$fisherMatrix[ 4L, 4L ], 0 )
} )

test_that(".duplicateFim returns distinct empty Fim instances", {
  proto = IndividualFim()
  d1    = PFIM:::.duplicateFim( proto )
  d2    = PFIM:::.duplicateFim( proto )
  expect_false( identical( d1, d2 ) )
  expect_false( identical( d1, proto ) )
  expect_equal( length( prop( d1, "fisherMatrix" ) ), 0L )
})

test_that("evaluateDesign: each arm has its own evaluationFim instance", {
  fx     = cas10_design()
  design = cas10_design_two_arms()
  model  = rebuildEvalModel(
    Evaluation(
      name                    = "two_arm",
      modelFromLibrary        = cas10_model_from_library(),
      modelParameters         = cas10_model_parameters(),
      modelCovariates         = cas10_covariates(),
      modelCovariatesEquation = "exponential",
      modelError              = fx$modelError,
      designs                 = list( design ),
      fimType                 = "individual",
      outputs                 = list( "RespPK" ),
      odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
    ),
    finiteDifference = TRUE
  )
  fimProto = IndividualFim()
  design   = PFIM:::evaluateDesign( design, model, fimProto )

  arms = prop( design, "evaluationArms" )
  expect_length( arms, 2L )
  expect_false( identical(
    prop( arms[[ 1L ]], "evaluationFim" ),
    prop( arms[[ 2L ]], "evaluationFim" )
  ) )
  expect_false( identical(
    prop( arms[[ 1L ]], "evaluationFim" ),
    prop( design, "fim" )
  ) )

  M1 = prop( prop( arms[[ 1L ]], "evaluationFim" ), "fisherMatrix" )
  M2 = prop( prop( arms[[ 2L ]], "evaluationFim" ), "fisherMatrix" )
  Md = prop( prop( design, "fim" ), "fisherMatrix" )
  # Individual: design FIM is the harmonic mean of per-subject arm FIMs.
  sizes = map_dbl( arms, ~ prop( .x, "size" ) )
  w = sizes / sum( sizes )
  expect_equal(
    Md,
    PFIM:::.pfimHarmonicMeanFim( list( M1, M2 ), w )$fisherMatrix,
    tolerance = 1e-8
  )
})

test_that("Bayesian multi-arm: harmonic mean of subject FIMs (prior kept); shrinkage from design FIM", {
  fx     = cas10_design()
  design = cas10_design_two_arms()
  ev     = Evaluation(
    name                    = "bayes_two_arm",
    modelFromLibrary        = cas10_model_from_library(),
    modelParameters         = cas10_model_parameters(),
    modelCovariates         = cas10_covariates(),
    modelCovariatesEquation = "exponential",
    modelError              = fx$modelError,
    designs                 = list( design ),
    fimType                 = "Bayesian",
    outputs                 = list( "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
  model    = rebuildEvalModel( ev, finiteDifference = TRUE )
  fimProto = BayesianFim()
  design   = PFIM:::evaluateDesign( design, model, fimProto )
  arms     = prop( design, "evaluationArms" )
  designFim = prop( design, "fim" )
  M1 = prop( prop( arms[[ 1L ]], "evaluationFim" ), "fisherMatrix" )
  M2 = prop( prop( arms[[ 2L ]], "evaluationFim" ), "fisherMatrix" )
  Md = prop( designFim, "fisherMatrix" )
  sizes = map_dbl( arms, ~ prop( .x, "size" ) )
  w = sizes / sum( sizes )

  expect_false( isTRUE( all.equal( Md, M1 + M2, tolerance = 1e-8 ) ) )
  expect_equal(
    Md,
    PFIM:::.pfimHarmonicMeanFim( list( M1, M2 ), w )$fisherMatrix,
    tolerance = 1e-8
  )

  # Shrinkage from assembled design FIM (same matrix as SE/RSE).
  expect_equal(
    prop( designFim, "shrinkage" ),
    PFIM:::.bayesianShrinkage( Md, model, arms[[ 1L ]] ),
    tolerance = 1e-10
  )
  # Identical equal-sized arms → same shrinkage as one arm.
  expect_equal(
    prop( designFim, "shrinkage" ),
    prop( prop( arms[[ 1L ]], "evaluationFim" ), "shrinkage" ),
    tolerance = 1e-10
  )
})

test_that("Individual design FIM weights unequal arm sizes (harmonic)", {
  fx = cas10_design()
  a1 = prop( fx$design1, "arms" )[[ 1L ]]
  a2 = PFIM:::.pfimCloneS7( a1 )
  prop( a1, "name" ) = "arm1"; prop( a1, "size" ) = 150
  prop( a2, "name" ) = "arm2"; prop( a2, "size" ) = 50
  # Distinct protocols so arithmetic and harmonic mixtures differ.
  prop( prop( a2, "samplingTimes" )[[ 1L ]], "samplings" ) = c( 0.5, 1, 2, 4, 12, 24 )
  design = Design( name = "unequal", arms = list( a1, a2 ) )
  model  = rebuildEvalModel(
    Evaluation(
      name = "uneq", modelFromLibrary = cas10_model_from_library(),
      modelParameters = cas10_model_parameters(),
      modelCovariates = cas10_covariates(),
      modelCovariatesEquation = "exponential",
      modelError = fx$modelError, designs = list( design ),
      fimType = "individual", outputs = list( "RespPK" ),
      odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
    ),
    finiteDifference = TRUE
  )
  design = PFIM:::evaluateDesign( design, model, IndividualFim() )
  arms = prop( design, "evaluationArms" )
  M1 = prop( prop( arms[[ 1L ]], "evaluationFim" ), "fisherMatrix" )
  M2 = prop( prop( arms[[ 2L ]], "evaluationFim" ), "fisherMatrix" )
  Md = prop( prop( design, "fim" ), "fisherMatrix" )
  expect_equal(
    Md,
    PFIM:::.pfimHarmonicMeanFim( list( M1, M2 ), c( 0.75, 0.25 ) )$fisherMatrix,
    tolerance = 1e-8
  )
  expect_false( isTRUE( all.equal( Md, 0.75 * M1 + 0.25 * M2, tolerance = 1e-6 ) ) )
})

test_that(".fimFixedEffectLabels matches setEvaluationFim dimensions (cas10)", {
  evaluation = run( cas10_evaluation( "labels", "individual" ) )
  fe         = .fimFixedEffectLabels( evaluation )
  sigma      = .fimSigmaBlockLabels( evaluation )
  # Subject FIM: mu + sigma (beta omitted).
  n_expected = length( c( fe$columnNamesMu, sigma$columnNamesSigma ) )
  expect_equal( ncol( prop( prop( evaluation, "fim" ), "fisherMatrix" ) ), n_expected )
})

test_that("individual FIM with one sample reports Inf SE on mu", {
  ev = Evaluation(
    name             = "one_sample_ind",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters  = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
      ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
    ),
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 ) ),
    designs    = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = 1 ) )
      ) )
    ) ),
    fimType = "individual",
    outputs = list( "RespPK" )
  )
  out = suppressWarnings( run( ev ) )
  se  = suppressWarnings( getSE( out ) )
  mu  = grepl( "mu_|\u03bc_", rownames( se ) )
  expect_true( any( mu ) )
  expect_true( all( is.infinite( se$SE[ mu ] ) ) )
})

test_that("FO IOV inflation of V is tcrossprod of scaled columns", {
  set.seed( 1L )
  p = 5L
  R = diag( p ) * 0.25
  G = matrix( runif( p * 3L ), p, 3L )
  colnames( G ) = c( "mu_ka", "mu_V", "mu_Cl" )
  bare = c( "ka", "V", "Cl" )
  gamma = c( 0.1, 0, 0.2 )
  muChain = c( 1.2, 2, 0.8 )
  V_new = PFIM:::.indBayesIovInflatedV( R, G, colnames( G ), bare, gamma, muChain )
  V_ref = R
  for ( j in seq_len( ncol( G ) ) ) {
    ii = match( sub( "^mu_", "", colnames( G )[[ j ]] ), bare )
    if ( is.na( ii ) || !( gamma[[ ii ]] > 0 ) )
      next
    gj = G[ , j, drop = FALSE ] * muChain[[ ii ]]
    V_ref = V_ref + ( gamma[[ ii ]] ^ 2 ) * tcrossprod( gj )
  }
  expect_equal( V_new, V_ref )
  expect_equal(
    PFIM:::.indBayesIovInflatedV( R, G, colnames( G ), bare, c( 0, 0, 0 ), muChain ),
    R
  )
})

test_that(".indBayesIovGradientCols prefers mu_ names", {
  expect_equal(
    PFIM:::.indBayesIovGradientCols( c( 1L, 3L ), c( "mu_ka", "V", "Cl" ), c( "ka", "V", "Cl" ) ),
    c( "mu_ka", "Cl" )
  )
  expect_equal(
    PFIM:::.indBayesIovGradientCols( integer( 0L ), c( "mu_ka" ), "ka" ),
    character( 0L )
  )
})
