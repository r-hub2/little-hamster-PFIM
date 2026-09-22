test_that("updateSamplingTimes keeps outcome-specific sampling grids", {
  arm = Arm(
    name = "a1", size = 10,
    administrations = list(
      Administration( outcome = "RespPK", timeDose = 0, dose = 100 )
    ),
    samplingTimes = list(
      SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 1, 2, 6 ) ),
      SamplingTimes( outcome = "RespPD", samplings = c( 0, 12, 24, 48 ) )
    )
  )
  samplingData = getSamplingData( arm )
  arm2 = updateSamplingTimes( arm, samplingData )

  pk = prop( arm2, "samplingTimes" )[[ 1L ]]
  pd = prop( arm2, "samplingTimes" )[[ 2L ]]
  expect_true( all( c( 12, 24, 48 ) %in% prop( pd, "samplings" ) ) )
  expect_false( any( c( 12, 24, 48 ) %in% prop( pk, "samplings" ) ) )
  expect_true( all( c( 0.5, 1, 2, 6 ) %in% prop( pk, "samplings" ) ) )
})

test_that( "response and SI plots densify while keeping design sampling times", {
  source( system.file( "examples", "evaluation-minimal.R", package = "PFIM" ) )
  plots = plotEvaluation( ev, list() )
  p = plots[[ 1L ]][[ 1L ]][[ 1L ]]
  expect_s3_class( p, "ggplot" )
  expect_gt( nrow( p$data ), 4L )
  expect_true( all( c( 0.33, 1.5, 5, 12 ) %in% round( p$data$time, 2 ) ) )
} )

test_that( "response and SI plots are skipped with a single sampling time", {
  source( system.file( "examples", "evaluation-minimal.R", package = "PFIM" ) )
  arm1 = prop( prop( ev, "designs" )[[ 1L ]], "arms" )[[ 1L ]]
  singleArm = Arm(
    name              = prop( arm1, "name" ),
    size              = prop( arm1, "size" ),
    administrations   = prop( arm1, "administrations" ),
    samplingTimes     = list( SamplingTimes( outcome = "RespPK", samplings = c( 10 ) ) )
  )
  model       = rebuildEvalModel( ev, finiteDifference = TRUE )
  outputNames = as.list( prop( model, "outputNames" ) )
  fim         = defineFim( ev )
  expect_false( .pfimArmPlotsEnabled( singleArm, model, outputNames ) )
  expect_length( processArmEvaluationResults( singleArm, model, fim, "design1", list() )$design1$arm1, 0L )
  expect_length( processArmEvaluationSI( singleArm, model, fim, "design1", list() )$design1$arm1, 0L )
  expect_true( .pfimArmPlotsEnabled( arm1, model, outputNames ) )
})
