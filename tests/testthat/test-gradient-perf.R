test_that("gradient performance options preserve FIM results (ODE 1cpt)", {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) )
  arm   = Arm(
    name = "arm", size = 40,
    administrations = list( admin ),
    samplingTimes   = list( st )
  )
  design = Design( name = "d", arms = list( arm ) )
  mp = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  err = Constant( output = "RespPK", sigmaInter = 0.1 )
  ev = Evaluation(
    name = "perf_test",
    modelFromLibrary = list( PKModel = "Linear1FirstOrderSingleDose_kaClV" ),
    modelParameters  = mp,
    modelError       = list( err ),
    outputs          = list( "RespPK" ),
    designs          = list( design ),
    fimType          = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )

  opts_off = list(
    perf.adminCache    = FALSE,
    perf.fdCache       = FALSE,
    perf.odeTimesCache = FALSE
  )
  opts_on = list(
    perf.adminCache    = TRUE,
    perf.fdCache       = TRUE,
    perf.odeTimesCache = TRUE
  )

  run_baseline = function( opts ) {
    old = setNames( lapply( names( opts ), pfim_get_option ), names( opts ) )
    do.call( pfim_set_option, opts )
    on.exit( do.call( pfim_set_option, old ), add = TRUE )
    out = run( ev )
    list(
      det = getDeterminant( out ),
      rse = getRSE( out )$RSE
    )
  }

  base = run_baseline( opts_off )
  fast = run_baseline( opts_on )

  expect_equal( fast$det, base$det, tolerance = 1e-6 )
  expect_equal( fast$rse, base$rse, tolerance = 1e-6 )
})
