path = system.file( "fixtures/eval-opt-references.R", package = "PFIM" )
if ( !nzchar( path ) )
  stop( "inst/fixtures/eval-opt-references.R not found." )
source( path, local = FALSE )

`%||%` = function( x, y ) if ( is.null( x ) || ( length( x ) == 0L && !is.list( x ) ) ) y else x

expect_fim_summary_matches_reference = function( stats, ref, tol, label = "" ) {
  prefix = if ( nzchar( label ) ) paste0( label, ": " ) else ""

  expect_equal( stats$d, ref$d, tolerance = tol$d %||% 1e-6,
                info = paste0( prefix, "D-criterion" ) )
  if ( !is.na( ref$det ) )
    expect_equal( stats$det, ref$det, tolerance = tol$det %||% 1e-6,
                  info = paste0( prefix, "det(FIM)" ) )
  if ( !is.na( ref$trace ) )
    expect_equal( stats$trace, ref$trace, tolerance = tol$trace %||% 0.1,
                  info = paste0( prefix, "trace(FIM)" ) )
  expect_equal( stats$ncol, ref$ncol, info = paste0( prefix, "ncol(FIM)" ) )
  if ( !is.na( ref$m11 ) )
    expect_equal( stats$m11, ref$m11, tolerance = tol$m11 %||% 1e-4,
                  info = paste0( prefix, "FIM[1,1]" ) )
  if ( !is.na( ref$m_nn ) )
    expect_equal( stats$m_nn, ref$m_nn, tolerance = tol$m_nn %||% 0.1,
                  info = paste0( prefix, "FIM[n,n]" ) )
  invisible( stats )
}
