# Shared fixtures for cas10-style covariate + IOV tests (PK 1-cpt, ka/V/Cl).

cas10_model_from_library = function() {
  list( PKModel = "Linear1FirstOrderSingleDose_kaClV" )
}

cas10_model_parameters = function() {
  list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1,   omega = sqrt( 0.09 ) ),
                    gamma = sqrt( 0.0225 ) ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ),
                    gamma = sqrt( 0.0225 ) ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = 2,   omega = sqrt( 0.09 ) ),
                    gamma = sqrt( 0.0225 ) )
  )
}

cas10_covariates = function() {
  sex = Covariate(
    name = "Sex", categories = c( "M", "F" ), categoriesProportions = c( 0.5, 0.5 ),
    effects = list( "F" = c( "V" = log( 1.2 ) ) )
  )
  treatment4 = Covariate(
    name = "Treatment", categories = c( "A", "B" ),
    sequences            = list( c( "A", "B", "A", "B" ), c( "B", "A", "B", "A" ) ),
    sequencesProportions = c( 0.5, 0.5 ),
    effects = list( "B" = c( "Cl" = log( 1.1 ) ) )
  )
  list( sex, treatment4 )
}

cas10_design = function() {
  errorModelRespPK = Constant( output = "RespPK", sigmaInter = 0.1 )
  administration   = Administration( outcome = "RespPK", timeDose = 0, dose = 30 )
  samplingTimes    = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6, 8 ) )
  arm1 = Arm(
    name = "arm1", size = 40,
    administrations = list( administration ),
    samplingTimes   = list( samplingTimes )
  )
  list(
    design1    = Design( name = "design1", arms = list( arm1 ) ),
    modelError = list( errorModelRespPK )
  )
}

#' Two identical cas10 arms (tests independent \code{evaluationFim} per arm).
cas10_design_two_arms = function() {
  fx    = cas10_design()
  arm1  = prop( fx$design1, "arms" )[[ 1L ]]
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 30 )
  st    = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6, 8 ) )
  arm2  = Arm(
    name = "arm2", size = 40,
    administrations = list( admin ),
    samplingTimes   = list( st )
  )
  Design( name = "design_two_arms", arms = list( arm1, arm2 ) )
}

cas10_evaluation = function( name, fimType ) {
  fx = cas10_design()
  Evaluation(
    name = name,
    modelFromLibrary        = cas10_model_from_library(),
    modelParameters         = cas10_model_parameters(),
    modelCovariates         = cas10_covariates(),
    modelCovariatesEquation = "exponential",
    modelError              = fx$modelError,
    designs                 = list( fx$design1 ),
    fimType                 = fimType,
    outputs                 = list( "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
}
