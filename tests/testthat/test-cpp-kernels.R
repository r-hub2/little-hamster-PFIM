# Lightweight checks on compiled optimizer kernels (no separate C++ test harness).

test_that("pso_optimize_Rcpp batch callback returns one cost per particle", {
  n_dims = 2L
  initial = c( 1, 4 )
  windows = list(
    matrix( c( 0, 10 ), nrow = 1 ),
    matrix( c( 0, 10 ), nrow = 1 )
  )
  groups = list( c( 1L, 2L ) )

  batch_hits = 0L
  eval_batch = function( pos ) {
    batch_hits <<- batch_hits + 1L
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

  expect_true( batch_hits >= 2L )
  expect_length( res$globalBestDesign, n_dims )
  expect_true( is.finite( res$globalBestCost ) )
} )
