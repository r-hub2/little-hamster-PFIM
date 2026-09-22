path = system.file( "fixtures/cov-iov-poster.R", package = "PFIM" )
if ( !nzchar( path ) )
  stop( "inst/fixtures/cov-iov-poster.R not found." )
source( path, local = FALSE )
