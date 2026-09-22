# Unit tests for numerical helpers and CRAN packaging invariants.

test_that("Dcriterion on empty or non-matrix FIM is 0", {
  expect_equal( PFIM:::Dcriterion( IndividualFim() ), 0 )
  expect_equal( PFIM:::.fimDcriterionFromMatrix( numeric( 0 ) ), 0 )
  expect_equal( PFIM:::.fimDcriterionFromMatrix( 1:6 ), 0 )
})

test_that("fisherMatrix must be a square matrix when non-empty", {
  fim = IndividualFim()
  expect_error(
    { prop( fim, "fisherMatrix" ) = as.numeric( 1:4 ) },
    "must be a matrix"
  )
  expect_error(
    { prop( fim, "fisherMatrix" ) = matrix( as.numeric( 1:6 ), 2L, 3L ) },
    "must be square"
  )
  prop( fim, "fisherMatrix" ) = diag( 2 )
  expect_equal( prop( fim, "fisherMatrix" ), diag( 2 ) )
})

test_that("residual error labels are sigma (SD), not sigma-squared", {
  expect_equal( PFIM:::.greekConsole[[ "sigma" ]], "\u03c3_" )
  expect_equal( PFIM:::.greekPlotmath[[ "sigma" ]], "sigma" )
  expect_equal(
    PFIM:::.pfimParamNameAscii( "\u03c3_inter_RespPK" ),
    "sigma_inter_RespPK"
  )
  expect_equal(
    PFIM:::.pfimSeRseFacetLabel( "RSE", "sigma" ),
    "'RSE'~'  '~sigma"
  )
})

test_that("getDcriterion on identity FIM equals 1", {
  fim = IndividualFim()
  prop( fim, "fisherMatrix" ) = diag( 3L )
  expect_equal( PFIM:::Dcriterion( fim ), 1 )
})

test_that("getDeterminant matches Dcriterion^(p) on the FIM", {
  fim = IndividualFim()
  M   = diag( c( 4, 9, 16 ) )
  prop( fim, "fisherMatrix" ) = M
  expect_equal( PFIM:::.fimDeterminant( M ), PFIM:::Dcriterion( fim )^3 )
})

test_that("condition number is Inf for indefinite FIM", {
  M = matrix( c( 1, 2, 2, -1 ), 2, 2 )
  expect_equal( PFIM:::.conditionNumber( M ), Inf )
})

test_that("condition number is Inf for rank-deficient PSD FIM", {
  M = matrix( c( 1, 0, 0, 0 ), 2, 2 )
  expect_equal( PFIM:::.conditionNumber( M ), Inf )
})

test_that("LibraryOfPKModels class is not shadowed by default instance", {
  expect_true( is.function( LibraryOfPKModels ) )
  expect_s7_class( PFIM:::pkModelLibrary, LibraryOfPKModels )
  expect_s7_class( PFIM:::pdModelLibrary, LibraryOfPDModels )
  expect_gt( length( prop( PFIM:::pkModelLibrary, "models" ) ), 0L )
})

test_that("Administration aligns single timeDose for multiple doses", {
  admin = Administration( outcome = "RespPK", timeDose = c( 0 ), dose = c( 1, 2 ) )
  aligned = getFromNamespace( ".alignAdministrationDosing", "PFIM" )( admin )
  expect_equal( aligned$timeDose, c( 0, 0 ) )
  expect_equal( aligned$dose, c( 1, 2 ) )
  expect_s7_class(
    Administration( outcome = "RespPK", tau = 24, dose = 5500 ),
    Administration
  )
})

test_that("Administration validator rejects incompatible timeDose and dose lengths", {
  expect_error(
    Administration( outcome = "RespPK", timeDose = c( 0, 1 ), dose = c( 1, 2, 3 ) ),
    "length\\(timeDose\\)"
  )
})

test_that("Arm validator rejects negative size", {
  expect_error( Arm( name = "a1", size = -1 ), "size must be non-negative" )
})

test_that("Arm validator rejects non-Administration administrations", {
  expect_error(
    Arm( name = "a1", administrations = list( list( dose = 1 ) ) ),
    "Administration"
  )
})

test_that("Design validator rejects non-Arm arms", {
  expect_error(
    Design( name = "d1", arms = list( list( name = "a1" ) ) ),
    "Arm"
  )
})

test_that("SamplingTimeConstraints accepts scalar minSampling with multiple windows", {
  expect_s7_class(
    SamplingTimeConstraints(
      outcome                      = "RespPK",
      samplingsWindows             = list( c( 1, 2 ), c( 3, 4 ) ),
      numberOfTimesByWindows       = c( 2, 2 ),
      minSampling                  = 0.5
    ),
    SamplingTimeConstraints
  )
})

test_that("SamplingTimeConstraints validator rejects inconsistent window vectors", {
  expect_error(
    SamplingTimeConstraints(
      outcome                      = "RespPK",
      samplingsWindows             = list( c( 0, 1 ), c( 2, 3 ) ),
      numberOfTimesByWindows       = c( 1, 2, 3 ),
      minSampling                  = 0.5
    ),
    "numberOfTimesByWindows length"
  )
})

test_that("SamplingTimeConstraints requires windows in increasing time order", {
  expect_error(
    SamplingTimeConstraints(
      outcome = "RespPK",
      initialSamplings = c( 0.5, 2, 3, 6, 8 ),
      numberOfTimesByWindows = c( 2, 3 ),
      samplingsWindows = list( c( 4, 12 ), c( 0.1, 4 ) ),
      minSampling = 0.05
    ),
    "increasing order"
  )
})

test_that("SamplingTimeConstraints validator rejects too many optimisable samplings", {
  expect_error(
    SamplingTimeConstraints(
      outcome                      = "RespPK",
      initialSamplings             = c( 1, 2 ),
      numberOfsamplingsOptimisable = 3
    ),
    "numberOfsamplingsOptimisable"
  )
})

test_that("SamplingTimeConstraints rejects initialSamplings that miss the windows", {
  expect_error(
    SamplingTimeConstraints(
      outcome = "RespPK",
      initialSamplings = c( 0.25, 1, 3, 12, 24 ),
      numberOfTimesByWindows = c( 2, 3 ),
      samplingsWindows = list( c( 0.1, 6 ), c( 6, 30 ) ),
      minSampling = 0.05
    ),
    "do not place"
  )
})

test_that("contiguous windows assign a shared endpoint to the earlier window", {
  expect_s7_class(
    SamplingTimeConstraints(
      outcome = "RespPK",
      initialSamplings = c( 0.5, 2, 4, 6, 8 ),
      numberOfTimesByWindows = c( 3, 2 ),
      samplingsWindows = list( c( 0.1, 4 ), c( 4, 12 ) ),
      minSampling = 0.05
    ),
    SamplingTimeConstraints
  )
  expect_s7_class(
    SamplingTimeConstraints(
      outcome = "RespPK",
      initialSamplings = c( 1, 2, 4, 6, 12, 20 ),
      numberOfTimesByWindows = c( 2, 2, 2 ),
      samplingsWindows = list( c( 0, 4 ), c( 4, 12 ), c( 12, 24 ) ),
      minSampling = 0.05
    ),
    SamplingTimeConstraints
  )
  expect_error(
    SamplingTimeConstraints(
      outcome = "RespPK",
      initialSamplings = c( 0.5, 2, 13 ),
      numberOfTimesByWindows = c( 2, 1 ),
      samplingsWindows = list( c( 0.1, 4 ), c( 4, 12 ) ),
      minSampling = 0.05
    ),
    "outside every sampling window"
  )
  sc = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = c( 0.5, 2, 4, 6, 8 ),
    numberOfTimesByWindows = c( 3, 2 ),
    samplingsWindows = list( c( 0.1, 4 ), c( 4, 12 ) ),
    minSampling = 0.05
  )
  arm = Arm(
    name = "a", size = 10,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6, 8 ) )
    ),
    samplingTimesConstraints = list( sc )
  )
  chk = checkSamplingTimeConstraintsForMetaheuristic(
    sc, arm, c( 0.5, 2, 4, 6, 8 ), "RespPK"
  )
  expect_true( chk$constraintWindowsLength )
  expect_true( chk$constraintMinimalSampling )
})

test_that("SamplingTimeConstraints rejects initialSamplings that violate minSampling", {
  expect_error(
    SamplingTimeConstraints(
      outcome = "RespPK",
      initialSamplings = c( 1.0, 1.1, 3.5, 3.6 ),
      numberOfTimesByWindows = c( 2, 2 ),
      samplingsWindows = list( c( 1, 2 ), c( 3, 4 ) ),
      minSampling = 0.5
    ),
    "minSampling"
  )
})

test_that("SamplingTimes validator rejects non-finite times", {
  expect_s7_class( SamplingTimes(), SamplingTimes )
  expect_error(
    SamplingTimes( outcome = "RespPK", samplings = c( 1, Inf ) ),
    "samplings must be finite"
  )
})

test_that("AdministrationConstraints validator rejects negative doses", {
  expect_s7_class( AdministrationConstraints(), AdministrationConstraints )
  expect_error(
    AdministrationConstraints( outcome = "RespPK", doses = list( -1 ) ),
    "non-negative numeric"
  )
})

test_that("Distribution validator rejects negative omega", {
  expect_error( LogNormal( mu = 1, omega = -0.1 ), "omega must be non-negative" )
})

.pfimCtorBody = function( x ) paste( deparse( body( x ) ), collapse = "\n" )

test_that("new_object uses positional parent (CRAN S7 .parent and devel _parent)", {
  # Installed constructors (always available under R CMD check; no R/ sources).
  expect_match( .pfimCtorBody( Combined1 ), "new_object\\(\\s*ModelError\\(" )
  expect_match( .pfimCtorBody( Combined2 ), "new_object\\(\\s*ModelError\\(" )
  expect_match( .pfimCtorBody( Constant ), "new_object\\(\\s*ModelError\\(" )
  expect_match( .pfimCtorBody( Proportional ), "new_object\\(\\s*ModelError\\(" )
  expect_match( .pfimCtorBody( Additive ), "new_object\\(\\s*CovariateModelEquation\\(" )
  expect_match( .pfimCtorBody( Exponential ), "new_object\\(\\s*CovariateModelEquation\\(" )
  expect_match( .pfimCtorBody( CategoricalCovariate ), "new_object\\(\\s*Covariate\\(" )
  expect_match(
    .pfimCtorBody( CategoricalCovariateWithIOV ),
    "new_object\\(\\s*Covariate\\("
  )
  expect_match( .pfimCtorBody( ModelError ), "new_object\\(\\s*S7_object\\(" )
  expect_match( .pfimCtorBody( Covariate ), "new_object\\(\\s*S7_object\\(" )
  expect_match( .pfimCtorBody( CovariateModelEquation ), "new_object\\(\\s*S7_object\\(" )
  expect_match( .pfimCtorBody( SamplingTimeConstraints ), "new_object\\(\\s*S7_object\\(" )
  expect_match( .pfimCtorBody( Evaluation ), "new_object\\(\\s*PFIMProject\\(" )
  expect_match( .pfimCtorBody( Optimization ), "new_object\\(\\s*S7_object\\(" )
  ctor_src = paste(
    purrr::map_chr(
      list(
        Combined1, Combined2, Constant, Proportional, ModelError,
        Additive, Exponential, CovariateModelEquation, Covariate,
        CategoricalCovariate, CategoricalCovariateWithIOV,
        SamplingTimeConstraints, Evaluation, Optimization
      ),
      .pfimCtorBody
    ),
    collapse = "\n"
  )
  expect_false( grepl( "new_object\\(\\s*\\.(parent|_parent)\\s*=", ctor_src ) )
  expect_false( grepl( "new_object\\(\\s*ModelError\\s*,", ctor_src ) )
  expect_false( grepl( "new_object\\(\\s*Covariate\\s*,", ctor_src ) )
  expect_false( grepl( "new_object\\(\\s*CovariateModelEquation\\s*,", ctor_src ) )
  expect_false( grepl( "new_object\\(\\s*SamplingTimeConstraints\\s*,", ctor_src ) )
  expect_false( grepl( "new_object\\(\\s*CategoricalCovariate\\s*,", ctor_src ) )
} )

test_that("R/ sources do not pass .parent / _parent to new_object", {
  r_dir = testthat::test_path( "../../R" )
  skip_if_not( dir.exists( r_dir ), "package R/ sources not available" )
  hits = list.files( r_dir, pattern = "\\.[Rr]$", full.names = TRUE ) |>
    purrr::keep( ~ grepl(
      "new_object\\s*\\(\\s*\\.(parent|_parent)\\s*=",
      paste( readLines( .x, warn = FALSE ), collapse = "\n" )
    ) ) |>
    purrr::map_chr( basename )
  expect_equal( hits, character() )
} )

test_that("S7 constructors build instances under the #409 parent rule", {
  expect_s7_class( ModelError( output = "y" ), ModelError )
  expect_s7_class(
    Combined1( output = "y", sigmaInter = 0.1, sigmaSlope = 0.1 ),
    Combined1
  )
  expect_s7_class(
    Combined2( output = "y", sigmaInter = 0.1, sigmaSlope = 0.1 ),
    Combined2
  )
  expect_s7_class( Constant( output = "y", sigmaInter = 0.1 ), Constant )
  expect_s7_class( Proportional( output = "y", sigmaSlope = 0.1 ), Proportional )
  expect_s7_class( CovariateModelEquation(), CovariateModelEquation )
  expect_s7_class( Additive(), Additive )
  expect_s7_class( Exponential(), Exponential )
  expect_s7_class( Covariate( name = "wt" ), Covariate )
  expect_s7_class(
    CategoricalCovariate(
      name = "sex",
      categories = c( "F", "M" ),
      categoriesProportions = c( 0.5, 0.5 )
    ),
    CategoricalCovariate
  )
  expect_s7_class(
    CategoricalCovariateWithIOV(
      name = "trt",
      categories = c( "A", "B" ),
      sequences = list( c( "A", "B" ), c( "B", "A" ) ),
      sequencesProportions = c( 0.5, 0.5 )
    ),
    CategoricalCovariateWithIOV
  )
  expect_s7_class( SamplingTimeConstraints(), SamplingTimeConstraints )
} )

test_that("ggplot2 facet_wrap space support is cached in .pfimSession", {
  pfim_reset_session()
  expect_false(
    exists( "ggplot2.supportsFacetSpace", envir = PFIM:::.pfimSession, inherits = FALSE )
  )
  first = PFIM:::.pfimGgplot2SupportsFacetSpace()
  expect_true(
    exists( "ggplot2.supportsFacetSpace", envir = PFIM:::.pfimSession, inherits = FALSE )
  )
  expect_identical( PFIM:::.pfimGgplot2SupportsFacetSpace(), first )
  expect_type( first, "logical" )
} )

test_that("Optimization uses composition (project slot)", {
  opt = Optimization(
    name = "opt_test",
    fimType = "individual",
    fim = IndividualFim(),
    designs = list()
  )
  expect_equal( projectProp( opt, "fimType" ), "individual" )
  expect_s7_class( prop( opt, "project" ), PFIMProject )
  expect_false( S7::S7_inherits( opt, PFIMProject ) )
})

test_that("Covariate is exported from NAMESPACE", {
  expect_true( "Covariate" %in% getNamespaceExports( "PFIM" ) )
})

test_that("pipeline generics are not on the public API", {
  ns = getNamespaceExports( "PFIM" )
  expect_true( "run" %in% ns )
  expect_true( "Report" %in% ns )
  expect_true( "getDcriterion" %in% ns )
  expect_true( "plotSensitivityIndices" %in% ns )
  expect_true( "FedorovWynnAlgorithm" %in% ns )
  expect_true( "MultiplicativeAlgorithm" %in% ns )
  expect_false( "Dcriterion" %in% ns )
  expect_false( "generateFimsFromConstraints" %in% ns )
  expect_false( "optimizeDesign" %in% ns )
  expect_false( "setEvaluationFim" %in% ns )
  expect_false( "evaluateDesign" %in% ns )
  expect_false( "pfim_cache_stats" %in% ns )
})

test_that("pfim_set_option rejects unknown names and unnamed args", {
  expect_error( pfim_set_option( not.a.real.option = TRUE ), "unknown option" )
  expect_error( pfim_set_option( 1 ), "named arguments" )
  empty_name = list( FALSE )
  names( empty_name ) = ""
  expect_error( do.call( pfim_set_option, empty_name ), "named arguments" )
  pfim_set_option( verbose = FALSE )
  expect_false( isTRUE( pfim_get_option( "verbose" ) ) )
})

test_that( ".pfimIsNonEmptyScalar rejects NA and empty", {
  expect_false( PFIM:::.pfimIsNonEmptyScalar( NA_character_ ) )
  expect_false( PFIM:::.pfimIsNonEmptyScalar( "" ) )
  expect_false( PFIM:::.pfimIsNonEmptyScalar( character( 0 ) ) )
  expect_true( PFIM:::.pfimIsNonEmptyScalar( "arm1" ) )
  expect_true( PFIM:::.pfimIsBlankScalar( NA_character_ ) )
  expect_true( PFIM:::.pfimIsBlankScalar( "" ) )
  expect_false( PFIM:::.pfimIsBlankScalar( character( 0 ) ) )
  expect_false( PFIM:::.pfimIsBlankScalar( "arm1" ) )
} )

test_that( "NA FIM cache scope is ignored", {
  old = pfim_get_option( "fim.cache.scope", NULL )
  on.exit( pfim_set_option( fim.cache.scope = old ), add = TRUE )
  pfim_set_option( fim.cache.scope = NA_character_ )
  expect_identical( PFIM:::.pfimActiveCacheScope(), "global" )
  expect_null( PFIM:::.pfimFimCacheClear( NA_character_ ) )
  expect_null( PFIM:::.pfimKableAsis( NULL ) )
} )

test_that("show is not exported from PFIM (use methods::show)", {
  expect_false( "show" %in% getNamespaceExports( "PFIM" ) )
  expect_true( "show" %in% getNamespaceExports( "methods" ) )
  expect_true( "showFIM" %in% getNamespaceExports( "PFIM" ) )
  # Internal S7 external-generic binding may exist; must not be exported.
  expect_error( PFIM::show, "not an exported object|could not find" )
  expect_true( is.function( methods::show ) )
})

test_that( "testthat condition helpers do not take assignment via named `=`", {
  files = list.files( testthat::test_path(), pattern = "[.]R$", full.names = TRUE )
  hits = .pfimTestFilesWithNamedAssign( files )
  expect_equal( length( hits ), 0L, info = paste( basename( hits ), collapse = ", " ) )
})

test_that("tost is exported and documented as covariateTest alias", {
  expect_true( "tost" %in% getNamespaceExports( "PFIM" ) )
  expect_identical( getFromNamespace( "tost", "PFIM" ), covariateTest )
  rd = testthat::test_path( "../../man/covariateTest.Rd" )
  skip_if_not( file.exists( rd ), "source Rd unavailable in installed-package check" )
  lines = readLines( rd, warn = FALSE )
  expect_true( any( grepl( "^\\\\alias\\{tost\\}", lines ) ) )
})

test_that("generateDosesCombination enumerates dose grids", {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 10 )
  ac = AdministrationConstraints( outcome = "RespPK", doses = as.list( c( 10, 20 ) ) )
  arm = Arm(
    name = "a1",
    size = 10,
    administrations = list( admin ),
    administrationsConstraints = list( ac )
  )
  design = Design( name = "d1", arms = list( arm ) )
  out = PFIM:::generateDosesCombination( design )
  expect_equal( out$numberOfDoses, 2L )
  expect_equal( out$a1$RespPK, c( 10, 20 ) )
})

test_that( "proportional subject allocation sums exactly to total", {
  a = PFIM:::.allocProportionalSubjects( 100, c( 1, 1, 1 ) )
  expect_equal( sum( a ), 100 )
  expect_equal( a, c( 34, 33, 33 ) )
  expect_equal(
    PFIM:::.allocProportionalSubjects( 300, c( 1, 2, 1 ) ),
    c( 75, 150, 75 )
  )
  b = PFIM:::.allocProportionalSubjects( 100, c( 0.3003, 0.6997 ) )
  expect_equal( sum( b ), 100 )
  expect_equal( b, c( 30, 70 ) )
})

test_that( "population multiplicative N uses numberOfSubjects", {
  algo = MultiplicativeAlgorithm()
  prop( algo, "multiplicativeAlgorithmOutputs" ) = list(
    numberOfSubjects = 100,
    weightsIndex = c( 25L, 49L ),
    optimalWeights = c( 0.3003, 0.6997 )
  )
  cell = list( entries = list( list( arm = Arm( name = "x", size = 50 ) ) ) )
  cells = setNames( rep( list( cell ), 49 ), as.character( seq_len( 49 ) ) )
  N = PFIM:::.multiplicativePopulationNTotal( algo, cells, c( 25L, 49L ) )
  expect_equal( N, 100 )
  expect_equal(
    PFIM:::.allocProportionalSubjects( N, c( 0.3003, 0.6997 ) ),
    c( 30, 70 )
  )
})
