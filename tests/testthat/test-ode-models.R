# ODE model dispatch and end-to-end FIM evaluation (library PK models with Deriv_*).

test_that("ODE dispatch: ModelODEDoseNotInEquations (1-cpt MM bolus)", {
  pk_model = "MichaelisMenten1BolusSingleDose_VmKm"
  params = list(
    ModelParameter( name = "Vm", distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km", distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) )
  )
  ev = .ode_evaluation(
    pk_model          = pk_model,
    model_parameters  = params,
    administration    = Administration( outcome = "C1", timeDose = 0, dose = 100 ),
    sampling_times    = SamplingTimes(
      outcome = "C1",
      samplings = c( 0, 0.5, 1, 2, 6, 9, 12, 24 )
    ),
    outputs           = list( RespPK = "C1" ),
    initial_condition = list( C1 = 0 ),
    name              = "ode_mm_bolus"
  )
  .expect_ode_model_class( ev, pk_model, ModelODEDoseNotInEquations )
  ev = run( ev )
  .expect_valid_fim( ev )
})

test_that("ODE dispatch: ModelODEDoseInEquations (1-cpt MM first-order)", {
  pk_model = "MichaelisMenten1FirstOrderSingleDose_kaVmKmV"
  params = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1.0, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 15.00, omega = sqrt( 0.25 ) ) ),
    ModelParameter( name = "Vm", distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km", distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) )
  )
  ev = .ode_evaluation(
    pk_model          = pk_model,
    model_parameters  = params,
    administration    = Administration( outcome = "RespPK", timeDose = 0, dose = 100 ),
    sampling_times    = SamplingTimes(
      outcome = "RespPK",
      samplings = c( 0, 0.33, 1.5, 3, 5, 8, 11, 12 )
    ),
    outputs           = list( RespPK = "C1" ),
    initial_condition = list( C1 = 0 ),
    name              = "ode_mm_fo"
  )
  .expect_ode_model_class( ev, pk_model, ModelODEDoseInEquations )
  ev = run( ev )
  .expect_valid_fim( ev )
})

test_that("ODE MM first-order accepts Administration(outcome = C1)", {
  pk_model = "MichaelisMenten1FirstOrderSingleDose_kaVmKmV"
  params = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1.0, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 15.00, omega = sqrt( 0.25 ) ) ),
    ModelParameter( name = "Vm", distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km", distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) )
  )
  ev = .ode_evaluation(
    pk_model          = pk_model,
    model_parameters  = params,
    administration    = Administration( outcome = "C1", timeDose = 0, dose = 100 ),
    sampling_times    = SamplingTimes( outcome = "C1", samplings = c( 0.33, 1.5, 3, 5, 8 ) ),
    outputs           = list( RespPK = "C1" ),
    initial_condition = list( C1 = 0 ),
    name              = "ode_mm_fo_c1"
  )
  ev = run( ev )
  .expect_valid_fim( ev )
})

test_that("ODE dispatch: ModelODEBolus (dose in initial condition)", {
  pk_model = "MichaelisMenten1BolusSingleDose_VmKm"
  params = list(
    ModelParameter( name = "Vm", distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km", distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 15, omega = sqrt( 0.25 ) ) )
  )
  ev = .ode_evaluation(
    pk_model          = pk_model,
    model_parameters  = params,
    administration    = Administration( outcome = "C1", timeDose = 0, dose = 100 ),
    sampling_times    = SamplingTimes( outcome = "C1", samplings = c( 1, 4, 12 ) ),
    outputs           = list( RespPK = "C1" ),
    initial_condition = list( C1 = "dose_C1/V" ),
    name              = "ode_mm_bolus_ic"
  )
  .expect_ode_model_class( ev, pk_model, ModelODEBolus )
})

test_that("ODE dispatch: ModelODEDoseNotInEquations (2-cpt MM bolus)", {
  pk_model = "MichaelisMenten2BolusSingleDose_VmKmk12k21V1V2"
  params = list(
    ModelParameter( name = "Vm",  distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km",  distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) ),
    ModelParameter( name = "V1",  distribution = LogNormal( mu = 10, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "V2",  distribution = LogNormal( mu = 20, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "k12", distribution = LogNormal( mu = 0.5, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "k21", distribution = LogNormal( mu = 0.3, omega = sqrt( 0.20 ) ) )
  )
  ev = .ode_evaluation(
    pk_model          = pk_model,
    model_parameters  = params,
    administration    = Administration( outcome = "C1", timeDose = 0, dose = 100 ),
    sampling_times    = SamplingTimes(
      outcome = "C1",
      samplings = c( 0, 0.5, 1, 2, 6, 9, 12, 24 )
    ),
    outputs           = list( RespPK = "C1" ),
    initial_condition = list( C1 = 0, C2 = 0 ),
    name              = "ode_mm_2cpt"
  )
  .expect_ode_model_class( ev, pk_model, ModelODEDoseNotInEquations )
  ev = run( ev )
  .expect_valid_fim( ev )
})

test_that("ODE dispatch: ModelODEInfusionDoseInEquation (MM infusion)", {
  pk_model = "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2"
  params = list(
    ModelParameter( name = "Vm",  distribution = LogNormal( mu = 0.08, omega = sqrt( 0.10 ) ) ),
    ModelParameter( name = "Km",  distribution = LogNormal( mu = 0.40, omega = sqrt( 0.30 ) ) ),
    ModelParameter( name = "V1",  distribution = LogNormal( mu = 10, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "V2",  distribution = LogNormal( mu = 20, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "k12", distribution = LogNormal( mu = 0.5, omega = sqrt( 0.20 ) ) ),
    ModelParameter( name = "k21", distribution = LogNormal( mu = 0.3, omega = sqrt( 0.20 ) ) )
  )
  ev = .ode_evaluation(
    pk_model          = pk_model,
    model_parameters  = params,
    administration    = Administration(
      outcome = "RespPK", Tinf = 2, timeDose = 0, dose = 50
    ),
    sampling_times    = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 6, 12 ) ),
    outputs           = list( RespPK = "C1" ),
    initial_condition = list( C1 = 0, C2 = 0 ),
    name              = "ode_mm_infusion"
  )
  model = .expect_ode_model_class( ev, pk_model, ModelODEInfusionDoseInEquation )
  # OOP: infusion ODE is a ModelODE (numerical parent), not ModelInfusion.
  expect_true( S7::S7_inherits( model, ModelODE ) )
  expect_true( S7::S7_inherits( model, ModelODEInfusion ) )
  expect_false( S7::S7_inherits( model, ModelInfusion ) )
  expect_true( .pfimIsOdeModel( model ) )
  expect_true( .pfimIsInfusionModel( model ) )
})

test_that( "ODE evaluation with default odeSolverParameters = list() runs", {
  pk_model = "Linear1FirstOrderSingleDose_kaClV"
  params = list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = 0.3 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = 0.3 ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  )
  ev = Evaluation(
    name = "ode_default_tol",
    modelFromLibrary = list( PKModel = pk_model ),
    modelParameters  = params,
    modelError       = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
    outputs          = list( RespPK = "RespPK" ),
    designs = list( Design(
      name = "d",
      arms = list( Arm(
        name = "a", size = 40,
        administrations = list( Administration( outcome = "RespPK", timeDose = 0, dose = 100 ) ),
        samplingTimes   = list( SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6 ) ) )
      ) )
    ) ),
    fimType = "population",
    odeSolverParameters = list()
  )
  ev = run( ev )
  .expect_valid_fim( ev )
} )

test_that("Linear2BolusSteadyState uses geometric dose accumulation", {
  fns = list(
    Linear2BolusSteadyState_ClQV1V2tau,
    Linear2BolusSteadyState_kk12k21Vtau
  )
  purrr::walk( fns, function( fn ) {
    expr = fn()
    expect_match( expr, "1-exp\\(-", perl = TRUE )
    expect_false( grepl( "exp\\(1-exp", expr, perl = TRUE ) )
  } )
})

test_that( ".pfimAlignOdeStates follows equation order, not IC list order", {
  ic = PFIM:::.pfimAlignOdeStates( c( Ad = 1, Ac = 2 ), c( "Ac", "Ad" ) )
  expect_equal( names( ic ), c( "Ac", "Ad" ) )
  expect_equal( unname( ic ), c( 2, 1 ) )
  ic_missing = PFIM:::.pfimAlignOdeStates( c( Ac = 3 ), c( "Ac", "Ad" ) )
  expect_equal( unname( ic_missing ), c( 3, 0 ) )
  expect_error(
    PFIM:::.pfimAlignOdeStates( c( A = 0, BB = 50 ), c( "A", "B" ) ),
    "not in the ODE states"
  )
} )

test_that( "ODE bolus y follows Deriv_* order when the dose is on the second state", {
  fo_eval = function( equations, initial_conditions ) {
    Evaluation(
      name = "fo_abs_order",
      modelEquations = equations,
      modelParameters = list(
        ModelParameter( name = "ka", distribution = LogNormal( mu = 1, omega = sqrt( 0.2 ) ) ),
        ModelParameter(
          name = "k1", distribution = LogNormal( mu = 0.25, omega = 0 ),
          fixedMu = TRUE, fixedOmega = TRUE
        )
      ),
      modelError = list( Constant( output = "RespPK", sigmaInter = 0.1 ) ),
      outputs = list( RespPK = "Ac" ),
      designs = list( Design( name = "d", arms = list( Arm(
        name = "a", size = 200,
        administrations = list(
          Administration( outcome = "Ad", timeDose = 0, dose = 100 )
        ),
        samplingTimes = list(
          SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 2, 4, 8, 12 ) )
        ),
        initialConditions = initial_conditions
      ) ) ) ),
      fimType = "population",
      odeSolverParameters = .ode_solver_params()
    )
  }
  eqs_central_first = list(
    Deriv_Ac = "ka*Ad - k1*Ac",
    Deriv_Ad = "-ka*Ad"
  )
  eqs_depot_first = list(
    Deriv_Ad = "-ka*Ad",
    Deriv_Ac = "ka*Ad - k1*Ac"
  )
  simulate_abs = function( equations, initial_conditions ) {
    ev = fo_eval( equations, initial_conditions )
    model = PFIM:::rebuildEvalModel( ev, finiteDifference = FALSE )
    arm = prop( prop( ev, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
    model = defineModelAdministration( model, arm )
    PFIM:::.odeSimulateBolus( model, arm, "doseEvent" )
  }
  ev_setup = fo_eval( eqs_central_first, list( Ad = 0, Ac = 0 ) )
  model = PFIM:::rebuildEvalModel( ev_setup, finiteDifference = FALSE )
  arm = prop( prop( ev_setup, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  model = defineModelAdministration( model, arm )
  expect_equal( prop( model, "variableNames" ), c( "Ac", "Ad" ) )
  expect_equal( names( prop( model, "initialConditions" ) ), c( "Ac", "Ad" ) )

  sim_central = simulate_abs( eqs_central_first, list( Ac = 0, Ad = 0 ) )
  sim_depot = simulate_abs( eqs_depot_first, list( Ac = 0, Ad = 0 ) )
  sim_ic_perm = simulate_abs( eqs_central_first, list( Ad = 0, Ac = 0 ) )
  expect_equal( sim_central$Ac, sim_depot$Ac, tolerance = 1e-6 )
  expect_equal( sim_central$Ad, sim_depot$Ad, tolerance = 1e-6 )
  expect_equal( sim_central$Ac, sim_ic_perm$Ac, tolerance = 1e-10 )

  ev_central = run( fo_eval( eqs_central_first, list( Ac = 0, Ad = 0 ) ) )
  ev_depot = run( fo_eval( eqs_depot_first, list( Ac = 0, Ad = 0 ) ) )
  ev_ic_perm = run( fo_eval( eqs_central_first, list( Ad = 0, Ac = 0 ) ) )
  .expect_valid_fim( ev_central )
  .expect_valid_fim( ev_depot )
  expect_equal( getDcriterion( ev_central ), getDcriterion( ev_depot ), tolerance = 1e-4 )
  expect_equal( getDcriterion( ev_central ), getDcriterion( ev_ic_perm ), tolerance = 1e-8 )
} )
