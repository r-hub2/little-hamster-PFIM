# ODE sim-time cache (C++).

test_that( "pfimCacheHash_Rcpp is stable for the same parts", {
  parts = c( "bolusIc", "0.5,2,4,6", "RespPK:0:100" )
  expect_equal( pfimCacheHash_Rcpp( parts ), pfimCacheHash_Rcpp( parts ) )
  expect_type( pfimCacheHash_Rcpp( parts ), "character" )
} )

test_that( "pfimEnvCacheKey_Rcpp shortens long keys", {
  short = paste( rep( "a", 100L ), collapse = "" )
  long  = paste( rep( "b", 600L ), collapse = "" )
  expect_equal( pfimEnvCacheKey_Rcpp( short ), short )
  expect_true( nchar( pfimEnvCacheKey_Rcpp( long ) ) < nchar( long ) )
} )

test_that( "pfimOdeSimTimesCached_Rcpp matches sorted unique grid", {
  pfimOdeSimTimesCacheClear_Rcpp()
  raw = c( 0, 0.5, 2, 4, 6 )
  ev  = data.frame( time = c( 0, 1, 3 ), value = c( 100, 50, 50 ) )
  key = pfimCacheHash_Rcpp( c(
    "doseEvent", paste( raw, collapse = "," ),
    paste( ev$time, ev$value, sep = ":", collapse = "," )
  ) )

  ref = sort( unique( c( raw, ev$time ) ) )
  ref = ref[ c( TRUE, diff( ref ) > 1e-8 ) ]

  out = pfimOdeSimTimesCached_Rcpp( key, raw, ev$time, enabled = FALSE )
  expect_equal( as.numeric( out ), ref, tolerance = 1e-12 )

  pfimOdeSimTimesCacheClear_Rcpp()
  out1 = pfimOdeSimTimesCached_Rcpp( key, raw, ev$time, enabled = TRUE )
  out2 = pfimOdeSimTimesCached_Rcpp( key, raw, ev$time, enabled = TRUE )
  expect_equal( as.numeric( out1 ), ref, tolerance = 1e-12 )
  expect_equal( as.numeric( out2 ), ref, tolerance = 1e-12 )
  expect_equal( pfimOdeSimTimesCacheSize_Rcpp(), 1L )
} )

test_that( "pfimOdeSimTimesCached_Rcpp respects maxEntries LRU cap", {
  pfimOdeSimTimesCacheClear_Rcpp()
  raw = c( 0, 0.5, 2, 4, 6 )
  ev  = data.frame( time = c( 0, 1, 3 ), value = c( 100, 50, 50 ) )
  k1  = pfimCacheHash_Rcpp( "ode_lru_1" )
  k2  = pfimCacheHash_Rcpp( "ode_lru_2" )
  k3  = pfimCacheHash_Rcpp( "ode_lru_3" )

  pfimOdeSimTimesCached_Rcpp( k1, raw, ev$time, enabled = TRUE, max_entries = 2L )
  pfimOdeSimTimesCached_Rcpp( k2, raw, ev$time, enabled = TRUE, max_entries = 2L )
  pfimOdeSimTimesCached_Rcpp( k3, raw, ev$time, enabled = TRUE, max_entries = 2L )
  expect_equal( pfimOdeSimTimesCacheSize_Rcpp(), 2L )

  pfimOdeSimTimesCached_Rcpp( k1, raw, ev$time, enabled = TRUE, max_entries = 2L )
  expect_equal( pfimOdeSimTimesCacheSize_Rcpp(), 2L )
} )
