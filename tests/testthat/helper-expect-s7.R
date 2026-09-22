# Compatible with testthat >= 3.0 (expect_s7_class needs >= 3.2.3).
expect_s7_class = function( object, class ) {
  fn = get0( "expect_s7_class", envir = asNamespace( "testthat" ), inherits = FALSE )
  if ( !is.null( fn ) ) {
    fn( object, class )
  } else {
    expect_true( S7::S7_inherits( object, class ) )
  }
}
