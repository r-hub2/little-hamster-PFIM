# ODE bolus PKPD: run() + Report() (needs pandoc locally / CI).

test_that( "ODE bolus PKPD: run() then Report()", {
  .skip_if_no_pandoc()
  skip_on_cran()

  modelEquations = list(
    Deriv_Cc = "-Cl/V * Cc",
    Deriv_E  = "Rin * ( 1-(Imax*(Cc/V) )/( (Cc/V)+C50 ) )-kout*E"
  )
  modelParameters = list(
    ModelParameter( name = "V",    distribution = LogNormal( mu = 8,    omega = sqrt( 0.02 ) ) ),
    ModelParameter( name = "Cl",   distribution = LogNormal( mu = 0.13, omega = sqrt( 0.06 ) ) ),
    ModelParameter( name = "Rin",  distribution = LogNormal( mu = 5.4,  omega = sqrt( 0.2 ) ) ),
    ModelParameter( name = "kout", distribution = LogNormal( mu = 0.06, omega = sqrt( 0.02 ) ) ),
    ModelParameter( name = "Imax", distribution = LogNormal( mu = 1.0,  omega = sqrt( 0.1 ) ) ),
    ModelParameter( name = "C50",  distribution = LogNormal( mu = 1.2,  omega = sqrt( 0.01 ) ) )
  )
  modelError = list(
    Combined1( output = "RespPK", sigmaInter = 0.6, sigmaSlope = 0.07 ),
    Constant( output = "RespPD", sigmaInter = 4 )
  )
  administrationRespPK = Administration(
    outcome = "Cc", timeDose = c( 0, 20 ), dose = c( 100, 50 )
  )
  arm1 = Arm(
    name = "BrasTest1",
    size = 32,
    administrations = list( administrationRespPK ),
    samplingTimes = list(
      SamplingTimes( outcome = "Cc", samplings = c( 0.5, 1, 2, 6, 12, 24, 48, 120 ) ),
      SamplingTimes( outcome = "E",  samplings = c( 0, 24, 48, 72, 120 ) )
    ),
    initialConditions = list( Cc = "dose_Cc/V", E = "Rin/kout" )
  )
  design1 = Design( name = "design1", arms = list( arm1 ) )
  ev = Evaluation(
    name = "",
    modelEquations = modelEquations,
    modelParameters = modelParameters,
    modelError = modelError,
    designs = list( design1 ),
    fimType = "population",
    outputs = list( RespPK = "Cc", RespPD = "E" ),
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  )
  ev = run( ev )

  outdir  = file.path( tempdir(), "pfim-ode-report-test" )
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  outfile = "ode_pkpd.html"

  expect_error(
    Report( ev, outputPath = outdir, outputFile = outfile,
            plotOptions = list( unitTime = "h", unitOutcomes = c( "ng/mL", "unit" ) ) ),
    NA
  )
  expect_true( file.exists( file.path( outdir, outfile ) ) )
} )
