# Frozen PFIM gold numbers (inst/fixtures/gold-fim.R). Do not edit to silence tests.

.pfimLoadGold = function() {
  path = system.file( "fixtures", "gold-fim.R", package = "PFIM" )
  if ( !nzchar( path ) || !file.exists( path ) )
    path = file.path( "inst", "fixtures", "gold-fim.R" )
  e = new.env( parent = baseenv() )
  sys.source( path, envir = e )
  e$.pfimGold
}

.pfimGold = .pfimLoadGold()
