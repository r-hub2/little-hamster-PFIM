#' @title SamplingTimeConstraints
#' @description
#' Constraints on sampling times for design optimization (windows, fixed times,
#' number of points to optimize, minimum spacing).
#' @param outcome Character string: outcome name.
#' @param initialSamplings Initial sampling-time vector.
#' @param fixedTimes Times that must remain fixed during optimization.
#' @param numberOfsamplingsOptimisable Number of sampling times to optimize.
#' @param samplingsWindows List of time windows for each optimizable sample.
#' @param numberOfTimesByWindows Number of samples per window.
#' @param minSampling Minimum spacing between consecutive samples.
#' @export

SamplingTimeConstraints = new_class( "SamplingTimeConstraints", package = "PFIM",

                                     properties = list( outcome = new_property(class_character, default = character(0)),
                                                        initialSamplings = new_property(class_vector, default = numeric(0)),
                                                        fixedTimes = new_property(class_vector, default = numeric(0)),
                                                        numberOfsamplingsOptimisable = new_property(class_double, default = 0.0),
                                                        samplingsWindows = new_property(class_list, default = list()),
                                                        numberOfTimesByWindows = new_property(class_vector, default = numeric(0)),
                                                        minSampling = new_property(class_vector, default = numeric(0))),
                                     constructor = function( outcome = character(0),
                                                             initialSamplings = numeric(0),
                                                             fixedTimes = numeric(0),
                                                             numberOfsamplingsOptimisable = 0.0,
                                                             samplingsWindows = list(),
                                                             numberOfTimesByWindows = numeric(0),
                                                             minSampling = numeric(0) ) {
                                       new_object(
                                         SamplingTimeConstraints,
                                         outcome                      = outcome,
                                         initialSamplings             = initialSamplings,
                                         fixedTimes                   = fixedTimes,
                                         numberOfsamplingsOptimisable = numberOfsamplingsOptimisable,
                                         samplingsWindows             = samplingsWindows,
                                         numberOfTimesByWindows       = numberOfTimesByWindows,
                                         minSampling                  = minSampling
                                       )
                                     } )

generateSamplingsFromSamplingConstraints = new_generic( "generateSamplingsFromSamplingConstraints", c( "samplingTimeConstraints" ) )
checkSamplingTimeConstraintsForMetaheuristic = new_generic( "checkSamplingTimeConstraintsForMetaheuristic", c( "samplingTimesConstraints", "arm"  ) )

#' Draw feasible sampling times under constraints
#' @name generateSamplingsFromSamplingConstraints
#' @export

method( generateSamplingsFromSamplingConstraints, SamplingTimeConstraints ) = function( samplingTimeConstraints ) {

  # get min sampling constraint
  minSampling = prop( samplingTimeConstraints, "minSampling" )
  # get samplings window constraints
  samplingsWindow = prop( samplingTimeConstraints, "samplingsWindows" )
  # get numberOfTimesByWindows constraints
  numberOfTimesByWindows = prop( samplingTimeConstraints, "numberOfTimesByWindows" )
  # generate samplings from sampling constraints
  intervalsConstraints = list()

  minSamplingAndNumberOfTimesByWindows = as.data.frame( list( minSampling, numberOfTimesByWindows ) )
  data = t( as.data.frame( samplingsWindow ) )
  inputRandomSpaced = as.data.frame( do.call( "cbind", list( data, minSamplingAndNumberOfTimesByWindows ) ) )

  colnames( inputRandomSpaced ) = c("min","max","delta","n")
  rownames( inputRandomSpaced ) = NULL

  for( iter in seq_len( length( inputRandomSpaced$n ) ) )
  {
    min = inputRandomSpaced$min[iter]
    max = inputRandomSpaced$max[iter]
    delta = inputRandomSpaced$delta[iter]
    n = inputRandomSpaced$n[iter]

    distance = max-min-(n-1)*delta

    ind = runif(n,0,1)
    tmp = distance*sort( ind )
    interval = min + tmp + delta * seq(0,n-1,1)

    intervalsConstraints[[iter]] = interval
  }

  return( unlist( intervalsConstraints ) )
}

#' Check sampling feasibility for PSO, PGBO or simplex
#' @name checkSamplingTimeConstraintsForMetaheuristic
#' @export

method( checkSamplingTimeConstraintsForMetaheuristic, list( SamplingTimeConstraints, Arm ) ) = function( samplingTimesConstraints, arm, newSamplings, outcome ) {

  armName = prop(arm, 'name')

  # get min sampling constraint
  minSampling = prop( samplingTimesConstraints, "minSampling")

  # get samplings window constraints
  samplingsWindow =  prop( samplingTimesConstraints, "samplingsWindows")
  # get numberOfTimesByWindows constraints
  numberOfTimesByWindows =  prop( samplingTimesConstraints, "numberOfTimesByWindows")

  # get the constraints min, max, delta and n
  minSamplingAndNumberOfTimesByWindows = as.data.frame( list( minSampling, numberOfTimesByWindows ) )
  data = t( as.data.frame(samplingsWindow ) )
  inputRandomSpaced = as.data.frame( do.call( "cbind", list( data, minSamplingAndNumberOfTimesByWindows ) ) )
  colnames( inputRandomSpaced ) = c("min","max","delta","n")
  rownames( inputRandomSpaced ) = NULL

  testForConstraintsWindowsLength = list()
  testForConstraintsMinimalSampling = list()

  for( iter in seq_len( length( inputRandomSpaced$n ) ) )
  {
    min = inputRandomSpaced$min[iter]
    max = inputRandomSpaced$max[iter]

    # check the constraint numberOfTimesByWindows
    testForConstraintsWindowsLength[[ armName ]][[ outcome ]][[iter]] = newSamplings[ newSamplings >= min & newSamplings <= max ]

    # check the constraint minSampling
    testForConstraintsMinimalSampling[[ armName ]][[ outcome ]][[iter]] = diff( testForConstraintsWindowsLength[[ armName ]][[ outcome ]][[ iter ]] )

    # case for one sampling
    if ( length( testForConstraintsMinimalSampling[[ armName ]][[ outcome ]][[iter]] ) == 0 )
    {
      testForConstraintsMinimalSampling[[ armName ]][[ outcome ]][[iter]] = 0
    }
  }

  # tests for the constraint
  testForConstraintsWindowsLengthtmp = map_int( testForConstraintsWindowsLength[[armName]][[outcome]], ~ length(.x) )

  if ( all( testForConstraintsWindowsLengthtmp == numberOfTimesByWindows ) )
  {
    constraintWindowsLength = TRUE
  }else{
    constraintWindowsLength = FALSE
  }

  testForConstraintsMinimalSamplingTmp = map_dbl(testForConstraintsMinimalSampling[[armName]][[outcome]], ~ min(.x) )

  if ( all( testForConstraintsMinimalSamplingTmp >= minSampling ) )
  {
    constraintMinimalSampling = TRUE
  }else{
    constraintMinimalSampling = FALSE
  }
  return( list( constraintWindowsLength = constraintWindowsLength, constraintMinimalSampling = constraintMinimalSampling ) )
}
