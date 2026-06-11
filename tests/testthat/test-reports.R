# HTML report templates and Report() smoke tests.

test_that("all skeleton Rmd templates exist with exact case", {
  templates = c(
    "EvaluationPopulationFIM.Rmd",
    "EvaluationIndividualFIM.Rmd",
    "EvaluationBayesianFIM.Rmd",
    "OptimizationMultiplicativeAlgorithmPopulationFIM.Rmd",
    "OptimizationMultiplicativeAlgorithmIndividualFIM.Rmd",
    "OptimizationMultiplicativeAlgorithmBayesianFIM.Rmd",
    "OptimizationFedorovWynnAlgorithmPopulationFIM.Rmd",
    "OptimizationFedorovWynnAlgorithmIndividualFIM.Rmd",
    "OptimizationFedorovWynnAlgorithmBayesianFIM.Rmd",
    "OptimizationSimplexAlgorithmPopulationFIM.Rmd",
    "OptimizationSimplexAlgorithmIndividualFIM.Rmd",
    "OptimizationSimplexAlgorithmBayesianFIM.Rmd",
    "OptimizationPSOAlgorithmPopulationFIM.Rmd",
    "OptimizationPSOAlgorithmIndividualFIM.Rmd",
    "OptimizationPSOAlgorithmBayesianFIM.Rmd",
    "OptimizationPGBOAlgorithmPopulationFIM.Rmd",
    "OptimizationPGBOAlgorithmIndividualFIM.Rmd",
    "OptimizationPGBOAlgorithmBayesianFIM.Rmd"
  )
  base = system.file( "rmarkdown", "templates", "skeleton", package = "PFIM" )
  missing = templates[ !file.exists( file.path( base, templates ) ) ]
  expect_length( missing, 0L )
})

test_that("optimization report data includes constraints table key", {
  opt = .minimal_mult_opt()
  combos = generateSamplingTimesCombination(
    projectProp( opt, "designs" )[[ 1L ]]
  )
  expect_gte( length( combos$opt_arm ), 2L )
  opt = run( opt )
  od = prop( opt, "optimisationDesign" )
  expect_s7_class( od$evaluationOptimalDesign, Evaluation )
  expect_gt( getDcriterion( opt ), 0 )
})

test_that("Report() renders Fedorov-Wynn optimization HTML when rmarkdown is available", {
  skip_if_not_installed( "rmarkdown" )
  skip_on_cran()

  opt = run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_report_smoke" ) )
  outdir  = file.path( tempdir(), "pfim-fw-report-test" )
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  outfile = "fw_opt_smoke.html"

  expect_error(
    Report( opt, outputPath = outdir, outputFile = outfile, plotOptions = list() ),
    NA
  )
  expect_true( file.exists( file.path( outdir, outfile ) ) )
})

test_that("SamplingTimeConstraints uses numberOfsamplingsOptimisable", {
  sc = SamplingTimeConstraints(
    outcome                      = "RespPK",
    initialSamplings             = c( 1, 2, 3, 4 ),
    fixedTimes                   = 1,
    numberOfsamplingsOptimisable = 3
  )
  expect_equal( prop( sc, "numberOfsamplingsOptimisable" ), 3 )
})

test_that("Covariate is exported after library(PFIM)", {
  expect_true( "Covariate" %in% getNamespaceExports( "PFIM" ) )
  cov = Covariate(
    name = "CYP2C9",
    categories = c( "Wild", "Others" ),
    categoriesProportions = c( 0.6, 0.4 ),
    effects = list( "Others" = c( "Cl" = log( 0.5 ) ) )
  )
  expect_s7_class( cov, CategoricalCovariate )
})

test_that("Report() renders evaluation HTML when rmarkdown is available", {
  skip_if_not_installed( "rmarkdown" )
  skip_on_cran()

  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 1, 2, 4, 8 ) )
  arm   = Arm(
    name            = "rpt_arm",
    size            = 50,
    administrations = list( admin ),
    samplingTimes   = list( st )
  )
  design = Design( name = "rpt_design", arms = list( arm ) )
  ev = Evaluation(
    name             = "report_eval_smoke",
    modelFromLibrary = list( PKModel = "Linear1BolusSingleDose_kV" ),
    modelParameters  = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) ),
      ModelParameter( name = "V", distribution = LogNormal( mu = 15, omega = 0.3 ) )
    ),
    modelError       = list(
      Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 )
    ),
    outputs          = list( "RespPK" ),
    designs          = list( design ),
    fimType          = "individual"
  )
  ev = run( ev )

  outdir  = file.path( tempdir(), "pfim-report-test" )
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  outfile = "eval_smoke.html"

  expect_error(
    Report( ev, outputPath = outdir, outputFile = outfile, plotOptions = list() ),
    NA
  )
  expect_true( file.exists( file.path( outdir, outfile ) ) )
})
