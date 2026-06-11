# Unit tests for numerical helpers and CRAN packaging invariants.

test_that("getDcriterion on identity FIM equals 1", {
  fim = IndividualFim()
  prop( fim, "fisherMatrix" ) = diag( 3L )
  expect_equal( PFIM:::Dcriterion( fim ), 1 )
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
