test_that( "model feature detection uses word boundaries for tau", {
  eqs_tau_param = list( RespPK = "ka * tau_absorption * t" )
  eqs_steady    = list( RespPK = "dose_RespPK/V * exp(-k*t)/(1-exp(-k*tau))" )
  eqs_param_tau = list( RespPK = "dose_RespPK/V * exp(-tau * t)" )
  ic = character( 0L )

  expect_false( PFIM:::.detectModelFeatures( eqs_tau_param, ic )$hasTau )
  expect_false( PFIM:::.detectModelFeatures( eqs_steady, ic )$hasTau )
  expect_true( PFIM:::.detectModelFeatures( eqs_steady, ic, hasAdminTau = TRUE )$hasTau )
  # Parameter named tau in equations without dosing-interval admin → not SS.
  expect_false( PFIM:::.detectModelFeatures( eqs_param_tau, ic, hasAdminTau = FALSE )$hasTau )
  expect_true( PFIM:::.detectModelFeatures( eqs_param_tau, ic, hasAdminTau = TRUE )$hasTau )
} )

test_that( "TurnoverkoutEmax differs from TurnoverRinEmax", {
  pd   = prop( pdModelLibrary, "models" )
  rin  = pd[[ "TurnoverRinEmax_RinEmaxCC50koutE" ]][[ "Deriv_E" ]]
  kout = pd[[ "TurnoverkoutEmax_RinEmaxCC50koutE" ]][[ "Deriv_E" ]]
  expect_false( identical( rin, kout ) )
  expect_match( kout, "kout\\*\\(", perl = TRUE )
} )

test_that( "S7 validators reject non-physical error and parameter values", {
  expect_error(
    Combined1( output = "y", sigmaInter = -0.1, sigmaSlope = 0.1 ),
    "non-negative"
  )
  expect_error(
    ModelParameter(
      name = "Cl",
      distribution = LogNormal( mu = 0, omega = -0.1 )
    ),
    "omega must be non-negative"
  )
} )
