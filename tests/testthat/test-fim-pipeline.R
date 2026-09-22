# Test FIM pipeline

test_that("replaceVariablesLibraryOfModels preserves dose_RespPK and Tinf_RespPK", {
  text = paste(
    "dose_RespPK/V * ka * exp( -ka * t ) - Cl/V * RespPK",
    "Tinf_RespPK",
    sep = "; "
  )
  out = replaceVariablesLibraryOfModels( text, "RespPK", "C1" )
  expect_match( out, "dose_RespPK" )
  expect_match( out, "Tinf_RespPK" )
  expect_match( out, "Cl/V \\* C1" )
} )

test_that("remapPkpdLibraryText aligns dose tokens with compartment administration", {
  admin = Administration( outcome = "C1", timeDose = 0, dose = 10, Tinf = 2 )
  arm = Arm(
    name = "arm",
    size = 10,
    administrations = list( admin ),
    samplingTimes = list( SamplingTimes( outcome = "RespPK", samplings = c( 1, 2 ) ) ),
    initialConditions = list( C1 = 0, C2 = 0 )
  )
  ev = Evaluation(
    name = "dose_remap",
    modelFromLibrary = list( PKModel = "Linear1InfusionSingleDose_ClV" ),
    modelParameters = list(
      ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.1 ) ),
      ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = 0.1 ) )
    ),
    outputs = list( "RespPK" = "C1", "RespPD" = "C2" ),
    designs = list( Design( name = "d", arms = list( arm ) ) ),
    fimType = "population"
  )
  text = "dose_RespPK/Tinf_RespPK/Cl * (1 - exp(-Cl/V * t))"
  out = PFIM:::remapPkpdLibraryText( text, ev )
  expect_match( out, "dose_C1" )
  expect_match( out, "Tinf_C1" )
} )

test_that("remapOdePkLibraryText remaps dose_RespPK when C1 is administered", {
  admin = Administration( outcome = "C1", timeDose = 0, dose = 100 )
  arm = Arm(
    name = "arm",
    size = 10,
    administrations = list( admin ),
    samplingTimes = list( SamplingTimes( outcome = "C1", samplings = c( 1, 2 ) ) ),
    initialConditions = list( C1 = 0 )
  )
  ev = Evaluation(
    name = "ode_dose_remap",
    modelFromLibrary = list( PKModel = "MichaelisMenten1FirstOrderSingleDose_kaVmKmV" ),
    modelParameters = list(
      ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.1 ) ),
      ModelParameter( name = "V",  distribution = LogNormal( mu = 15, omega = 0.1 ) ),
      ModelParameter( name = "Vm", distribution = LogNormal( mu = 0.08, omega = 0.1 ) ),
      ModelParameter( name = "Km", distribution = LogNormal( mu = 0.4, omega = 0.1 ) )
    ),
    outputs = list( RespPK = "C1" ),
    designs = list( Design( name = "d", arms = list( arm ) ) ),
    fimType = "population"
  )
  text = "-Vm*C1/(Km+C1) + dose_RespPK/V*ka*exp(-ka*t)"
  out = PFIM:::remapOdePkLibraryText( text, ev )
  expect_match( out, "dose_C1" )
  expect_false( grepl( "dose_RespPK", out, fixed = TRUE ) )
} )

test_that(".evaluationFromProject preserves fim from optimization project", {
  opt = .minimal_mult_opt()
  projectProp( opt, "fim" )     = BayesianFim()
  projectProp( opt, "fimType" ) = "bayesian"
  ev = PFIM:::.evaluationFromProject( opt )
  expect_s7_class( prop( ev, "fim" ), BayesianFim )
} )

test_that("Bayesian optimization templates use evaluationForPlot for plotShrinkage", {
  templates = c(
    "OptimizationFedorovWynnAlgorithmBayesianFIM.Rmd",
    "OptimizationMultiplicativeAlgorithmBayesianFIM.Rmd",
    "OptimizationPGBOAlgorithmBayesianFIM.Rmd",
    "OptimizationPSOAlgorithmBayesianFIM.Rmd",
    "OptimizationSimplexAlgorithmBayesianFIM.Rmd"
  )
  base = system.file( "rmarkdown", "templates", "skeleton", package = "PFIM" )
  purrr::walk( templates, function( tmpl ) {
    lines = readLines( file.path( base, tmpl ), warn = FALSE )
    expect_true(
      any( grepl( "evaluationForPlot", lines, fixed = TRUE ) ),
      info = tmpl
    )
    expect_false(
      any( grepl( "plotShrinkage\\(.*pfimproject", lines ) ),
      info = tmpl
    )
  } )
} )

test_that("Bayesian optimization report renders when rmarkdown is available", {
  .skip_if_no_pandoc()
  skip_on_cran()

  opt = .minimal_mult_opt()
  projectProp( opt, "fim" )     = BayesianFim()
  projectProp( opt, "fimType" ) = "bayesian"
  opt = run( opt )

  outdir  = file.path( tempdir(), "pfim-bayes-opt-report-test" )
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  outfile = "bayes_opt_smoke.html"

  expect_error(
    Report( opt, outputPath = outdir, outputFile = outfile, plotOptions = list() ),
    NA
  )
  expect_true( file.exists( file.path( outdir, outfile ) ) )
} )
