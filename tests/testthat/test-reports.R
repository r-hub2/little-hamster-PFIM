# HTML report templates and Report() output.

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
  .skip_if_no_pandoc()
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
  sc2 = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = c( 1, 2, 3, 4 ),
    numberOfSamplingsOptimisable = 2
  )
  expect_equal( prop( sc2, "numberOfsamplingsOptimisable" ), 2 )
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

test_that("tablesForReport individual FIM omits beta; keeps mu+sigma (cas10)", {
  evaluation = run( cas10_evaluation( "cas10_tables", "individual" ) )
  fim        = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  tabs       = tablesForReport( fim, evaluation )
  fe         = as.matrix( prop( fim, "fixedEffects" ) )
  expect_equal( nrow( fe ), 3L )
  expect_length( tabs$fixedEffectsTable, 1L )
  expect_length( tabs$SEAndRSETable, 1L )
  html = paste( as.character( tabs$varianceEffectsTable ), collapse = "" )
  expect_match( html, "\\\\sigma_\\{inter_RespPK\\}" )
  expect_false( grepl( "\\{\\{inter", html, fixed = TRUE ) )
  expect_match( html, "table-striped" )
})

test_that("tablesForReport Bayesian FIM aligns shrinkage with mu rows (cas10)", {
  evaluation = run( cas10_evaluation( "cas10_bayes_tables", "Bayesian" ) )
  fim        = prop( evaluation, "fim" )
  tabs       = tablesForReport( fim, evaluation )
  fim        = PFIM:::setEvaluationFim( fim, evaluation )
  expect_length( tabs$fixedEffectsTable, 1L )
  expect_length( tabs$SEAndRSETable, 1L )
  expect_equal( ncol( prop( fim, "fisherMatrix" ) ), 3L )
  shrink = prop( fim, "shrinkage" )
  expect_equal( nrow( shrink ), 1L )
  expect_equal( ncol( shrink ), 3L )
  expect_equal( colnames( shrink ), .fimFixedEffectLabels( evaluation )$columnNamesMu )
  html = paste( as.character( tabs$SEAndRSETable ), collapse = "" )
  expect_true( grepl( "Shrinkage", html, fixed = TRUE ) )
  p = plotShrinkage( fim, evaluation )
  expect_s3_class( p, "ggplot" )
  expect_equal( p$data$Parameter, c( "ka", "V", "Cl" ) )
  expect_length( p$data$Shrinkage, 3L )
  expect_false( any( p$data$Parameter == "Shrinkage" ) )
})

test_that("Report() renders cas10 individual evaluation HTML when rmarkdown is available", {
  .skip_if_no_pandoc()
  skip_on_cran()

  evaluation = run( cas10_evaluation( "cas10_report", "individual" ) )
  outdir     = file.path( tempdir(), "pfim-cas10-report-test" )
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  outfile    = "cas10_eval.html"

  expect_error(
    Report( evaluation, outputPath = outdir, outputFile = outfile, plotOptions = list() ),
    NA
  )
  expect_true( file.exists( file.path( outdir, outfile ) ) )
})

test_that("Report() renders evaluation HTML when rmarkdown is available", {
  .skip_if_no_pandoc()
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

test_that( "plotSE population FIM includes beta bars when covariates are present", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim        = prop( evaluation, "fim" )
  se_n       = nrow( prop( PFIM:::setEvaluationFim( fim, evaluation ), "SEAndRSE" )$SEAndRSE )
  p          = plotSE( evaluation )
  expect_s3_class( p, "ggplot" )
  expect_equal( nrow( p$data ), se_n )
  expect_true( any( p$data$cat == "'SE'~'  '~beta" ) )
  expect_equal( unique( p$data$cat )[ 1L ], "'SE'~'  '~mu" )
})

test_that( "tablesForReport population FIM aligns SE/RSE rows with covariates", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim        = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  tabs       = tablesForReport( fim, evaluation )
  se_n       = nrow( prop( fim, "SEAndRSE" )$SEAndRSE )
  expect_length( tabs$SEAndRSETable, 1L )
  expect_error( plotRSE( evaluation ), NA )
  expect_equal( nrow( plotSE( evaluation )$data ), se_n )
})

test_that( "SE/RSE facet labels use plotmath Greek on all platforms", {
  expect_equal( PFIM:::.pfimSeRseFacetLabel( "SE",  "mu"    ), "'SE'~'  '~mu" )
  expect_equal( PFIM:::.pfimSeRseFacetLabel( "RSE", "beta"  ), "'RSE'~'  '~beta" )
  expect_equal( PFIM:::.pfimSeRseFacetLabel( "SE",  "omega" ), "'SE'~'  '~omega^2" )
  expect_equal( PFIM:::.pfimSeRseFacetLabel( "RSE", "sigma" ), "'RSE'~'  '~sigma" )

  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  p_se       = plotSE( evaluation )
  p_rse      = plotRSE( evaluation )
  expect_false( any( grepl( "^SE mu$|^SE beta$", p_se$data$cat ) ) )
  expect_true( all( grepl( "^'SE'~'  '~", p_se$data$cat ) ) )
  expect_true( all( grepl( "^'RSE'~'  '~", p_rse$data$cat ) ) )
  expect_s3_class( p_se + ggplot2::theme(), "ggplot" )
  # ggplot() construction does not evaluate layers; ggplot_build() does, and
  # catches ggplot2 3.5 unused-argument errors on facet_wrap(space=).
  expect_silent( ggplot2::ggplot_build( p_se ) )
  expect_silent( ggplot2::ggplot_build( p_rse ) )
})

test_that( "SE/RSE and base plot themes use 16 pt text", {
  se_theme = PFIM:::.pfimSeRseTheme()
  expect_equal( se_theme$axis.text.x$size, 16 )
  expect_equal( se_theme$axis.text.y$size, 16 )
  expect_equal( se_theme$axis.title.x$size, 16 )
  expect_equal( se_theme$strip.text.x$size, 16 )
  base_theme = PFIM:::.pfimBaseTheme()
  expect_equal( base_theme$axis.text.x$size, 16 )
  expect_equal( base_theme$axis.text.x$angle, 90 )
})

test_that( "report plot printer walks nested ggplot lists", {
  p = ggplot2::ggplot() + ggplot2::geom_blank()
  grDevices::pdf( nullfile() )
  on.exit( grDevices::dev.off(), add = TRUE )
  expect_silent( PFIM:::.pfimPrintReportPlots( list( a = list( b = p ) ) ) )
  expect_silent( PFIM:::.pfimPrintReportPlots( p ) )
  expect_silent( PFIM:::.pfimPrintReportPlots( list() ) )
})

test_that( "report skeleton plot chunks use 10 x 5 in figures", {
  templates = list.files(
    system.file( "rmarkdown", "templates", "skeleton", package = "PFIM" ),
    pattern = "\\.Rmd$", full.names = TRUE
  )
  expect_gt( length( templates ), 0L )
  bad_wide = templates[
    vapply( templates, function( f ) {
      any( grepl( "fig.width = 1[24]|fig.height = 8", readLines( f, warn = FALSE ) ) )
    }, logical( 1L ) )
  ]
  expect_length( bad_wide, 0L )
})

test_that( "tablesForReport population FIM aligns SE/RSE rows with covariates and IOV", {
  cas        = Filter( function( x ) x$id == "cas4_withCov_withIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim        = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  tabs       = tablesForReport( fim, evaluation )
  se_n       = nrow( prop( fim, "SEAndRSE" )$SEAndRSE )
  expect_length( tabs$SEAndRSETable, 1L )
  expect_error( plotRSE( evaluation ), NA )
  expect_equal( nrow( plotSE( evaluation )$data ), se_n )
})

test_that( "covariate report tables are built for population cov evaluation", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  modelCov   = prop( evaluation, "modelCovariates" )
  expect_gt( length( modelCov ), 0L )
  tbl = PFIM:::.buildCovariatesKable( modelCov )
  expect_false( is.null( tbl ) )
  flags = PFIM:::.covReportFlags( TRUE, NULL, tbl )
  expect_true( flags$showCovariates )
  mp = getModelParametersData( prop( evaluation, "modelParameters" )[[ 1L ]] )
  expect_equal( mp$mu_fixed, "FALSE" )
})

test_that( "covReportFlags omit empty covariate-test sections", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  ct_tables  = PFIM:::.buildCovariateTestSection( evaluation )
  flags      = PFIM:::.covReportFlags( TRUE, ct_tables, PFIM:::.buildCovariatesKable( prop( evaluation, "modelCovariates" ) ) )

  expect_true( flags$showCovariates )
  expect_true( flags$showCovTestSignificance )
  expect_true( flags$showCovTestNonRelevance )
  expect_false( flags$showCovTestRelevance )

  cas1 = Filter( function( x ) x$id == "cas1_noCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  ev1  = run( build_cov_iov_evaluation( cas1 ) )
  f1   = PFIM:::.covReportFlags( FALSE, NULL )
  expect_false( f1$showCovariates )
  expect_false( f1$showCovTestSignificance )
  expect_false( f1$showCovTestNonRelevance )
  expect_false( f1$showCovTestRelevance )
})

test_that( "population FIM report criteria table uses D-criterion label", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  fim        = PFIM:::setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  tabs       = tablesForReport( fim, evaluation )
  html       = paste( as.character( tabs$FIMCriteriaTable ), collapse = "" )
  expect_match( html, "D-criterion" )
  expect_false( grepl( "d-criterion", html, fixed = TRUE ) )
})

test_that( "pfimKableAsis returns non-empty HTML for covariate tables", {
  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  cov_asis   = PFIM:::.pfimKableAsis(
    PFIM:::.buildCovariatesKable( prop( evaluation, "modelCovariates" ) )
  )
  ct_asis    = PFIM:::.pfimKableAsis(
    getCovariateTestTables( covariateTest( evaluation ) )$nonRelevance
  )
  expect_s3_class( cov_asis, "knit_asis" )
  expect_s3_class( ct_asis, "knit_asis" )
  expect_true( grepl( "<table", as.character( cov_asis ), fixed = TRUE ) )
  expect_true( grepl( "<table", as.character( ct_asis ), fixed = TRUE ) )
})

test_that( "model error report table omits cError column", {
  err  = Constant( output = "RespPK", sigmaInter = 0.1 )
  html = as.character( PFIM:::.buildModelErrorKable( list( err ) ) )
  expect_false( grepl( "> NA <", html, fixed = TRUE ) )
  expect_false( grepl( "c_\\{error\\}", html ) )
  expect_match( html, "sigma_\\{inter\\}" )
  expect_match( html, "> 0.1 <" )
})

test_that( "sensitivity plot labels use plotmath Greek for mu and beta", {
  expect_equal( PFIM:::.pfimParamPlotmathLabel( "mu_ka" ), "mu[ka]" )
  expect_equal( PFIM:::.pfimParamPlotmathLabel( "\u03bc_ka" ), "mu[ka]" )
  expect_equal(
    PFIM:::.pfimParamPlotmathLabel( "beta_Cl_Treatment_B" ),
    "beta['Cl_Treatment_B']"
  )
  expect_equal( PFIM:::.pfimSensitivityYLab( "mu_ka" ), "frac(df, d*mu[ka])" )

  cas        = Filter( function( x ) x$id == "cas3_withCov_noIOV", cov_iov_evaluation_cases() )[[ 1L ]]
  evaluation = run( build_cov_iov_evaluation( cas ) )
  si         = plotSensitivityIndices( evaluation, list() )
  arm_plots  = si[[ 1L ]][[ 1L ]][[ 1L ]]
  beta_keys = names( arm_plots )[
    startsWith( names( arm_plots ), "beta_" ) |
      startsWith( names( arm_plots ), "\u03b2_" )
  ]
  expect_length( beta_keys, 0L )
  p = arm_plots[[ "mu_ka" ]]
  expect_s3_class( p, "ggplot" )
  y_lab = paste( as.character( p$labels$y ), collapse = " " )
  expect_false( grepl( "df/dmu", y_lab, fixed = TRUE ) )
  expect_match( y_lab, "mu\\[ka\\]" )
})

test_that( "show(Evaluation) dispatches via showFIM for individual FIM", {
  ev  = run( cas10_evaluation( "show_ind", "individual" ) )
  out = capture.output( show( ev ) )
  expect_true( any( grepl( "Individual Fisher Matrix", out ) ) )
  expect_false( any( grepl( "Population Fisher Matrix", out ) ) )
} )

test_that( "show(Evaluation) dispatches via showFIM for Bayesian FIM", {
  ev  = run( cas10_evaluation( "show_bay", "Bayesian" ) )
  out = capture.output( show( ev ) )
  expect_true( any( grepl( "Bayesian Fisher Matrix", out ) ) )
  expect_false( any( grepl( "Population Fisher Matrix", out ) ) )
} )
