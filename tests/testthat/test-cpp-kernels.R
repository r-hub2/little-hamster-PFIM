# Lightweight checks on compiled optimizer kernels (no separate C++ test harness).

.pfim_init_c_path = function() {
  candidates = c(
    testthat::test_path( "../../src/init.c" ),
    file.path( system.file( package = "PFIM" ), "..", "src", "init.c" )
  )
  hits = candidates[ file.exists( candidates ) ]
  if ( length( hits ) ) hits[[ 1L ]] else NA_character_
}

.pfim_registered_call_symbols = function() {
  init_path = .pfim_init_c_path()
  init = readLines( init_path, warn = FALSE )
  hits = grep( '^\\s*\\{"_PFIM_', init, value = TRUE )
  sub( '^\\s*\\{"([^"]+)".*', "\\1", hits )
}

.pfim_rcpp_call_symbols = function() {
  rcpp_path = testthat::test_path( "../../R/RcppExports.R" )
  lines = grep( " <- function", readLines( rcpp_path, warn = FALSE ), value = TRUE )
  paste0( "_PFIM_", sub( " <- function.*", "", lines ) )
}

test_that( "init.c registers every Rcpp .Call export", {
  skip_if_not( !is.na( .pfim_init_c_path() ), "src/init.c absent (installed tarball)" )
  expect_setequal( .pfim_registered_call_symbols(), .pfim_rcpp_call_symbols() )
} )

test_that("pso_optimize_Rcpp batch callback returns one cost per particle", {
  n_dims = 2L
  initial = c( 1, 4 )
  windows = list(
    matrix( c( 0, 10 ), nrow = 1 ),
    matrix( c( 0, 10 ), nrow = 1 )
  )
  groups = list( c( 1L, 2L ) )

  batch_state = new.env( parent = emptyenv() )
  batch_state$hits = 0L
  eval_batch = function( pos ) {
    batch_state$hits = batch_state$hits + 1L
    expect_equal( ncol( pos ), n_dims )
    expect_true( nrow( pos ) >= 1L )
    rep( 0.5, nrow( pos ) )
  }

  res = PFIM:::pso_optimize_Rcpp(
    n_pop_in           = 3L,
    max_iter           = 1L,
    initial_pos        = initial,
    windows_list       = windows,
    sorting_groups     = groups,
    phi1               = 2.05,
    phi2               = 2.05,
    constriction       = 0.7298,
    show_process       = FALSE,
    eval_fitness       = function( x ) 1,
    eval_fitness_batch = eval_batch
  )

  expect_true( batch_state$hits >= 2L )
  expect_length( res$globalBestDesign, n_dims )
  expect_true( is.finite( res$globalBestCost ) )
} )

test_that( "pso_optimize_Rcpp keeps fixed dimensions when windows are empty", {
  n_dims = 2L
  initial = c( 3.5, 7.25 )
  windows = list(
    matrix( numeric( 0 ), nrow = 0, ncol = 2 ),
    matrix( c( 0, 10 ), nrow = 1 )
  )
  groups = list( c( 1L, 2L ) )

  res = PFIM:::pso_optimize_Rcpp(
    n_pop_in           = 2L,
    max_iter           = 1L,
    initial_pos        = initial,
    windows_list       = windows,
    sorting_groups     = groups,
    phi1               = 2.05,
    phi2               = 2.05,
    constriction       = 0.7298,
    show_process       = FALSE,
    eval_fitness       = function( x ) sum( x^2 ),
    eval_fitness_batch = NULL
  )

  expect_equal( res$globalBestDesign[ 1L ], initial[ 1L ] )
  expect_true( res$globalBestDesign[ 2L ] >= 0 && res$globalBestDesign[ 2L ] <= 10 )
} )
