# Poster cov/IOV fixtures (1-cpt PK, population FIM).

cov_iov_model_equations = function() {
  list(
    RespPK = "dose_RespPK/V * ka/(ka - Cl/V) * (exp(-Cl/V * t) - exp(-ka * t))"
  )
}

cov_iov_model_error = function() {
  list( Constant( output = "RespPK", sigmaInter = 0.1 ) )
}

cov_iov_parameters = function( cl_mu = 2, with_iov = FALSE ) {
  gamma = if ( with_iov ) sqrt( 0.0225 ) else 0
  list(
    ModelParameter( name = "ka", distribution = LogNormal( mu = 1,   omega = sqrt( 0.09 ) ), gamma = gamma ),
    ModelParameter( name = "V",  distribution = LogNormal( mu = 3.5, omega = sqrt( 0.09 ) ), gamma = gamma ),
    ModelParameter( name = "Cl", distribution = LogNormal( mu = cl_mu, omega = sqrt( 0.09 ) ), gamma = gamma )
  )
}

cov_iov_sex = function() {
  Covariate(
    name = "Sex", categories = c( "M", "F" ), categoriesProportions = c( 0.5, 0.5 ),
    effects = list( F = c( V = log( 1.2 ) ) )
  )
}

cov_iov_treatment = function() {
  Covariate(
    name = "Treatment", categories = c( "A", "B" ),
    sequences = list( c( "A", "B" ), c( "B", "A" ) ),
    sequencesProportions = c( 0.5, 0.5 ),
    effects = list( B = c( Cl = log( 1.1 ) ) )
  )
}

cov_iov_evaluation_design = function() {
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 30 )
  sts   = SamplingTimes( outcome = "RespPK", samplings = c( 0.5, 2, 4, 6, 8 ) )
  arm   = Arm(
    name = "arm1", size = 40,
    administrations = list( admin ),
    samplingTimes   = list( sts )
  )
  Design( name = "design1", arms = list( arm ) )
}

cov_iov_simplex_design = function() {
  samps = c( 0, 2, 3, 8, 12, 24, 36, 50, 72, 120 )
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 30 )
  sts   = SamplingTimes( outcome = "RespPK", samplings = samps )
  stc   = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = samps,
    samplingsWindows = list( c( 0, 24 ), c( 35, 130 ) ),
    numberOfTimesByWindows = c( 6, 4 ),
    minSampling = c( 1, 2 )
  )
  arm = Arm(
    name = "arm1", size = 100,
    administrations = list( admin ),
    samplingTimes   = list( sts ),
    samplingTimesConstraints = list( stc )
  )
  Design( name = "design1", arms = list( arm ) )
}

cov_iov_mult_design = function() {
  samps = c( 0.5, 2, 4, 6, 8 )
  admin = Administration( outcome = "RespPK", timeDose = 0, dose = 30 )
  sts   = SamplingTimes( outcome = "RespPK", samplings = samps )
  stc   = SamplingTimeConstraints(
    outcome = "RespPK",
    initialSamplings = samps,
    numberOfsamplingsOptimisable = 3
  )
  adc = AdministrationConstraints( outcome = "RespPK", doses = list( 30 ) )
  arm = Arm(
    name = "arm1", size = 40,
    administrations = list( admin ),
    samplingTimes   = list( sts ),
    administrationsConstraints = list( adc ),
    samplingTimesConstraints   = list( stc )
  )
  Design( name = "design1", arms = list( arm ) )
}

cov_iov_evaluation_cases = function() {
  list(
    list(
      id = "cas1_noCov_noIOV", file = "cas1_noCov_noIOV.txt",
      cl_mu = 2, iov = FALSE, covariates = list(), cov_test = FALSE, ncol = 7L
    ),
    list(
      id = "cas2_noCov_withIOV", file = "cas2_noCov_withIOV.txt",
      cl_mu = 2, iov = TRUE, covariates = list(), cov_test = FALSE, ncol = 10L
    ),
    list(
      id = "cas3_withCov_noIOV", file = "cas3_withCov_noIOV.txt",
      cl_mu = 2, iov = FALSE, covariates = list( cov_iov_sex() ), cov_test = TRUE, ncol = 8L
    ),
    list(
      id = "cas4_withCov_withIOV", file = "cas4_withCov_withIOV.txt",
      cl_mu = 2, iov = TRUE, covariates = list( cov_iov_sex() ), cov_test = TRUE, ncol = 11L
    ),
    list(
      id = "cas5_withCov_mixed", file = "cas5_withCov_mixed.txt",
      cl_mu = 2, iov = TRUE,
      covariates = list( cov_iov_sex(), cov_iov_treatment() ),
      cov_test = TRUE, ncol = 12L, numberOfOccasions = 2
    ),
    list(
      id = "cas6_noCovFixed_withCovOccasion_withIOV",
      file = "cas6_noCovFixed_withCovOccasion_withIOV.txt",
      cl_mu = 2, iov = TRUE, covariates = list( cov_iov_treatment() ),
      cov_test = TRUE, ncol = 11L, numberOfOccasions = 2
    )
  )
}

build_cov_iov_evaluation = function( case ) {
  args = list(
    name                    = case$id,
    modelParameters         = cov_iov_parameters( cl_mu = case$cl_mu, with_iov = case$iov ),
    modelCovariates         = case$covariates,
    modelCovariatesEquation = "exponential",
    modelEquations          = cov_iov_model_equations(),
    modelError              = cov_iov_model_error(),
    designs                 = list( cov_iov_evaluation_design() ),
    fimType                 = "population",
    outputs                 = list( RespPK = "RespPK" ),
    odeSolverParameters     = list( atol = 1e-8, rtol = 1e-8 )
  )
  if ( !is.null( case$numberOfOccasions ) )
    args$numberOfOccasions = case$numberOfOccasions
  do.call( Evaluation, args )
}
