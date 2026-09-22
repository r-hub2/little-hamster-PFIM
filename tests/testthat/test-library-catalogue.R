# Catalogue equations: free symbols match the name suffix; closed forms parse.

.catalogue_symbols = function( expr ) {
  txt = unlist( expr, use.names = FALSE )
  unique( unlist( lapply( parse( text = txt ), all.names ), use.names = FALSE ) )
}

.catalogue_eval = function( expr, values ) {
  env = list2env( values, parent = baseenv() )
  txt = unlist( expr, use.names = FALSE )
  unname( vapply(
    txt,
    function( s ) as.numeric( eval( parse( text = s )[[ 1L ]], envir = env ) ),
    numeric( 1L )
  ) )
}

test_that( "Linear2InfusionSingleDose_ClQV1V2 inlines alpha/beta and matches micro-constants", {
  clq = Linear2InfusionSingleDose_ClQV1V2()
  kk  = Linear2InfusionSingleDose_kk12k21V()
  syms = .catalogue_symbols( clq )
  expect_false( any( c( "alpha", "beta" ) %in% syms ) )
  expect_true( all( c( "Cl", "Q", "V1", "V2" ) %in% syms ) )

  clq_vals = list(
    Cl = 2, Q = 1, V1 = 10, V2 = 20,
    dose_RespPK = 100, Tinf_RespPK = 1, t = 0.5
  )
  kk_vals = list(
    k = 2 / 10, k12 = 1 / 10, k21 = 1 / 20, V = 10,
    dose_RespPK = 100, Tinf_RespPK = 1, t = 0.5
  )
  expect_equal(
    .catalogue_eval( clq$duringInfusion, clq_vals ),
    .catalogue_eval( kk$duringInfusion, kk_vals ),
    tolerance = 1e-10
  )
  clq_vals$t = 2
  kk_vals$t = 2
  expect_equal(
    .catalogue_eval( clq$afterInfusion, clq_vals ),
    .catalogue_eval( kk$afterInfusion, kk_vals ),
    tolerance = 1e-10
  )
} )

test_that( "Linear2FirstOrderSteadyState_kakk12k21Vtau uses micro-constants only", {
  micro = Linear2FirstOrderSteadyState_kakk12k21Vtau()
  clq   = Linear2FirstOrderSteadyState_kaClQV1V2tau()
  syms = .catalogue_symbols( micro )
  expect_false( any( c( "alpha", "beta", "Cl", "Q", "V1", "V2" ) %in% syms ) )
  expect_true( all( c( "ka", "k", "k12", "k21", "V", "tau" ) %in% unique( syms ) ) )

  micro_vals = list(
    ka = 1, k = 0.2, k12 = 0.1, k21 = 0.05, V = 10,
    tau = 12, dose_RespPK = 100, t = 2
  )
  clq_vals = list(
    ka = 1, Cl = 2, Q = 1, V1 = 10, V2 = 20,
    tau = 12, dose_RespPK = 100, t = 2
  )
  expect_equal(
    .catalogue_eval( micro, micro_vals ),
    .catalogue_eval( clq, clq_vals ),
    tolerance = 1e-10
  )
} )

test_that( "Linear2InfusionSteadyState closed forms inline alpha/beta and parse", {
  kk  = Linear2InfusionSteadyState_kk12k21Vtau()
  clq = Linear2InfusionSteadyState_ClQV1V2tau()
  expect_false( any( c( "alpha", "beta" ) %in% .catalogue_symbols( kk ) ) )
  expect_false( any( c( "alpha", "beta" ) %in% .catalogue_symbols( clq ) ) )

  kk_vals = list(
    k = 0.2, k12 = 0.1, k21 = 0.05, V = 10,
    tau = 12, dose_RespPK = 100, Tinf_RespPK = 1, t = 0.5
  )
  clq_vals = list(
    Cl = 2, Q = 1, V1 = 10, V2 = 20,
    tau = 12, dose_RespPK = 100, Tinf_RespPK = 1, t = 0.5
  )
  expect_equal(
    .catalogue_eval( kk$duringInfusion, kk_vals ),
    .catalogue_eval( clq$duringInfusion, clq_vals ),
    tolerance = 1e-10
  )
  kk_vals$t = 2
  clq_vals$t = 2
  expect_equal(
    .catalogue_eval( kk$afterInfusion, kk_vals ),
    .catalogue_eval( clq$afterInfusion, clq_vals ),
    tolerance = 1e-10
  )
} )

test_that( "MichaelisMenten 1-cpt bolus vignette name is alias of VmKm", {
  pk = prop( pkModelLibrary, "models" )
  expect_identical(
    pk[[ "MichaelisMenten1BolusSingleDose_VmKmV" ]],
    pk[[ "MichaelisMenten1BolusSingleDose_VmKm" ]]
  )
} )

test_that( "unknown library model name errors", {
  arm = Arm(
    name = "a", size = 10,
    administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
    samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 1, 2 ) ) )
  )
  ev = Evaluation(
    name = "bad_lib",
    modelFromLibrary = list( PKModel = "NotALibraryModel" ),
    modelParameters = list(
      ModelParameter( name = "k", distribution = LogNormal( mu = 0.25, omega = 0.5 ) )
    ),
    modelError = list( Combined1( output = "RespPK", sigmaInter = 0.5, sigmaSlope = 0.15 ) ),
    designs = list( Design( name = "d", arms = list( arm ) ) ),
    fimType = "population",
    outputs = list( "RespPK" )
  )
  expect_error(
    PFIM:::defineModelEquationsFromLibraryOfModel( ev ),
    "was not found"
  )
} )

test_that( "MichaelisMenten 1-cpt infusion is the vignette ODE, old name is alias", {
  pk = prop( pkModelLibrary, "models" )
  canonical = pk[[ "MichaelisMenten1InfusionSingleDose_VmKmV" ]]
  alias     = pk[[ "MichaelisMenten1InfusionSingleDose_VmKmk12k21V1V2" ]]
  expect_identical( canonical, alias )
  expect_equal( names( canonical$duringInfusion ), "Deriv_C1" )
  expect_equal( names( canonical$afterInfusion ), "Deriv_C1" )
  expect_equal(
    canonical$duringInfusion$Deriv_C1,
    "-Vm*C1/(Km+C1) + dose_RespPK/V/Tinf_RespPK"
  )
  expect_equal( canonical$afterInfusion$Deriv_C1, "-Vm*C1/(Km+C1)" )
  syms = .catalogue_symbols( canonical )
  expect_false( any( c( "Cl", "k12", "k21", "V1", "V2", "C2" ) %in% syms ) )
  expect_true( all( c( "Vm", "Km", "V", "C1" ) %in% unique( syms ) ) )
} )

test_that( "MichaelisMenten 2-cpt infusion C1 keeps peripheral transfer", {
  pk = prop( pkModelLibrary, "models" )
  eq = pk[[ "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2" ]]
  expect_match( eq$duringInfusion$Deriv_C1, "k12", perl = TRUE )
  expect_match( eq$afterInfusion$Deriv_C1, "k12", perl = TRUE )
  expect_match( eq$duringInfusion$Deriv_C1, "k21", perl = TRUE )
} )

test_that( "MichaelisMenten 2-cpt ODEs use V1 (catalogue suffix), not V", {
  pk = prop( pkModelLibrary, "models" )
  keys = c(
    "MichaelisMenten2BolusSingleDose_VmKmk12k21V1V2",
    "MichaelisMenten2FirstOrderSingleDose_kaVmKmk12k21V1V2",
    "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2"
  )
  for ( key in keys ) {
    eq = pk[[ key ]]
    txt = paste( unlist( eq, use.names = FALSE ), collapse = " " )
    expect_false( grepl( "(^|[^[:alnum:]_])V([^[:alnum:]_]|$)", txt ) )
    expect_true( grepl( "(^|[^[:alnum:]_])V1([^[:alnum:]_]|$)", txt ) )
    expect_true( grepl( "(^|[^[:alnum:]_])V2([^[:alnum:]_]|$)", txt ) )
    syms = .catalogue_symbols( eq )
    expect_false( "V" %in% syms )
    expect_true( all( c( "V1", "V2" ) %in% unique( syms ) ) )
  }
} )

test_that( "quadratic PD uses Aquad and immediate full Imax keys exist", {
  pd = prop( pdModelLibrary, "models" )
  q = pd[[ "ImmediateDrugImaxQuadratic_S0AlinAquad" ]][[ "RespPD" ]]
  expect_true( "Aquad" %in% .catalogue_symbols( q ) )
  expect_false( grepl( "Alin \\* \\(RespPK\\)\\^2", q ) )
  expect_true( "ImmediateDrugFullImax_S0C50" %in% names( pd ) )
  expect_true( "ImmediateDrugSigmoidFullImax_S0C50gamma" %in% names( pd ) )
  expect_false( "Imax" %in% .catalogue_symbols( pd[[ "ImmediateDrugFullImax_S0C50" ]] ) )
} )

test_that( "sigmoid turnover keys list gamma and use ^", {
  pd = prop( pdModelLibrary, "models" )
  canonical = c(
    "TurnoverRinSigmoidEmax_RinEmaxCC50koutEgamma",
    "TurnoverRinSigmoidImax_RinImaxCC50koutEgamma",
    "TurnoverRinSigmoidFullImax_RinCC50koutEgamma",
    "TurnoverkoutSigmoidFullImax_RinCC50koutEgamma"
  )
  legacy = c(
    "TurnoverRinSigmoidEmax_RinEmaxCC50koutE",
    "TurnoverRinSigmoidImax_RinImaxCC50koutE",
    "TurnoverRinSigmoidFullImax_RinCC50koutE",
    "TurnoverkoutSigmoidFullImax_RinCC50koutE"
  )
  purrr::walk2( canonical, legacy, function( new_name, old_name ) {
    expect_true( new_name %in% names( pd ) )
    expect_identical( pd[[ new_name ]], pd[[ old_name ]] )
    eq = pd[[ new_name ]][[ "Deriv_E" ]]
    expect_false( grepl( "**", eq, fixed = TRUE ) )
    expect_true( "gamma" %in% .catalogue_symbols( eq ) )
  } )
} )

test_that( "Linear2 bolus wraps beta so B/beta is not B/0.5 * (...)", {
  clq = Linear2BolusSingleDose_ClQV1V2()
  kk  = Linear2BolusSingleDose_kk12k21V()
  expect_false( grepl( "/0.5\\*", clq ) )
  expect_false( grepl( "/0.5\\*", kk ) )
  expect_true( grepl( "/\\(0.5\\*", clq, perl = TRUE ) )
  expect_false( any( c( "alpha", "beta" ) %in% .catalogue_symbols( clq ) ) )
  clq_vals = list( Cl = 2, Q = 1, V1 = 10, V2 = 20, dose_RespPK = 100, t = 1 )
  kk_vals  = list( k = 0.2, k12 = 0.1, k21 = 0.05, V = 10, dose_RespPK = 100, t = 1 )
  expect_equal(
    .catalogue_eval( clq, clq_vals ),
    .catalogue_eval( kk, kk_vals ),
    tolerance = 1e-10
  )
} )

test_that( "Linear1 infusion SS Cl form uses /Cl, not ((Cl/V)*V)", {
  pk = prop( pkModelLibrary, "models" )
  eq = pk[[ "Linear1InfusionSteadyState_ClVtau" ]]
  txt = unlist( eq, use.names = FALSE )
  expect_false( any( grepl( "\\(Cl/V\\)\\s*\\*\\s*V", txt ) ) )
  expect_true( grepl( "/Cl", eq$duringInfusion$RespPK, fixed = TRUE ) )
} )

test_that( "exponential baseline increase rises and decrease decays", {
  pd = prop( pdModelLibrary, "models" )
  expect_equal(
    pd[[ "ImmediateBaselineExponentialincrease_S0kprog" ]][[ "RespPD" ]],
    "S0*(1-exp(-kprog*t))"
  )
  expect_equal(
    pd[[ "ImmediateBaselineExponentialdecrease_S0kprog" ]][[ "RespPD" ]],
    "S0*exp(-kprog*t)"
  )
} )
