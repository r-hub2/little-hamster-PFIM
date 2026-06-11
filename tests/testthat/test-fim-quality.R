# Quality gates: fimType validation, numerical reference values (cas10).

test_that("invalid fimType fails in run() with a clear message", {
  ev = cas10_evaluation( "bad_fim_type", "not_a_fim_type" )
  expect_error(
    run( ev ),
    regexp = "Invalid fimType|not_a_fim_type"
  )
})

test_that("defineFim rejects invalid fimType before design evaluation", {
  ev = cas10_evaluation( "bad_define", "invalid" )
  expect_error(
    defineFim( ev ),
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
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M          = prop( fim, "fisherMatrix" )

  expect_equal( det( M ), 5.830774, tolerance = 1e-4 )
  expect_equal( sum( diag( M ) ), 11936.63, tolerance = 0.15 )
  expect_equal( ncol( M ), 6L )
})

test_that("cas10 Bayesian: reference determinant, trace, and column labels", {
  evaluation = run( cas10_evaluation( "ref_bay", "Bayesian" ) )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  M          = prop( fim, "fisherMatrix" )

  expect_equal( det( M ), 2.197336e12, tolerance = 1e8 )
  expect_equal( sum( diag( M ) ), 3062.558, tolerance = 0.15 )
  expect_equal( ncol( M ), 5L )
  expect_equal(
    colnames( M ),
    c( "\u03bc_ka", "\u03bc_V", "\u03bc_Cl", "\u03b2_V_Sex_F", "\u03b2_Cl_Treatment_B" )
  )
  expect_gt( min( prop( fim, "shrinkage" ) ), 0 )
})

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
  design   = evaluateDesign( design, model, fimProto )

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
  expect_equal( Md, M1 + M2, tolerance = 1e-8 )
})

test_that(".fimFixedEffectLabels matches setEvaluationFim dimensions (cas10)", {
  evaluation = run( cas10_evaluation( "labels", "individual" ) )
  fe         = .fimFixedEffectLabels( evaluation )
  sigma      = .fimSigmaBlockLabels( evaluation )
  n_expected = length( c( fe$columnNamesMu, fe$columnNamesBeta, sigma$columnNamesSigma ) )
  expect_equal( ncol( prop( prop( evaluation, "fim" ), "fisherMatrix" ) ), n_expected )
})
