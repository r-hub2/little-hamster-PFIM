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
    ModelParameter( name = "V",   distribution = LogNormal( mu = 10, omega = sqrt( 0.20 ) ) ),
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
    ModelParameter( name = "V",   distribution = LogNormal( mu = 10, omega = sqrt( 0.20 ) ) ),
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
  .expect_ode_model_class( ev, pk_model, ModelODEInfusionDoseInEquation )
})
