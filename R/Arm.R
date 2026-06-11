#' @title Arm
#' @description
#' One experimental arm: subjects, dosing, sampling times, and slots filled during
#' design evaluation (model predictions, gradients, variance, FIM).
#' @param name Character string: arm identifier.
#' @param size Numeric: number of subjects in the arm.
#' @param administrations List of \code{Administration} objects.
#' @param initialConditions List of ODE initial conditions (name = variable, value = numeric).
#' @param samplingTimes List of \code{SamplingTimes} objects.
#' @param administrationsConstraints List of \code{AdministrationConstraints} objects.
#' @param samplingTimesConstraints List of \code{SamplingTimeConstraints} objects.
#' @param evaluationModel Model predictions at sampling times (filled by \code{evaluateArm}).
#' @param evaluationGradients Gradients of responses w.r.t. parameters (nested if covariates/IOV).
#' @param evaluationVariance Residual variance structure for the arm.
#' @param evaluationFim FIM object after \code{evaluateFim}.
#' @param covariateSampling Optional \code{CovariateSamplingMethods} object.
#' @include Model.R
#' @include Fim.R
#' @include MultiplicativeAlgorithm.R
#' @include FedorovWynnAlgorithm.R
#' @include SimplexAlgorithm.R
#' @include PSOAlgorithm.R
#' @include PGBOAlgorithm.R
#' @include CovariateSamplingMethods.R
#' @export

Arm = new_class("Arm", package = "PFIM",

                properties = list(
                  name = new_property(class_character, default = character(0)),
                  size = new_property(class_double, default = numeric(0)),
                  administrations = new_property(class_list, default = list()),
                  initialConditions = new_property(class_list, default = list()),
                  samplingTimes = new_property(class_list, default = list()),
                  administrationsConstraints = new_property(class_list, default = list()),
                  samplingTimesConstraints = new_property(class_list, default = list()),
                  evaluationModel = new_property(class_list, default = list()),
                  evaluationGradients = new_property(class_list, default = list()),
                  evaluationVariance = new_property(class_list, default = list()),
                  evaluationFim = new_property(Fim, default = NULL),
                  covariateSampling = new_property(CovariateSamplingMethods, default = NULL)
                ))



#' Evaluate model, gradients, variance and FIM for one arm
#' @param arm An \code{Arm} object.
#' @param ... \code{model} and \code{fim} for the \code{Arm} method.
#' @name evaluateArm
#' @export
evaluateArm = new_generic( "evaluateArm", c( "arm" ) )

#' Sampling times and grid extent for evaluation plots
#' @param arm An \code{Arm} object.
#' @param ... Not used.
#' @name getSamplingData
#' @export
getSamplingData = new_generic( "getSamplingData", c( "arm" ) )

#' Refine the time grid used in response and SI plots
#' @param arm An \code{Arm} object.
#' @param ... \code{samplingData} from \code{\link{getSamplingData}}.
#' @name updateSamplingTimes
#' @export
updateSamplingTimes = new_generic( "updateSamplingTimes", c( "arm" ) )

#' Model response plots for one arm
#' @param arm An \code{Arm} object.
#' @param model A \code{Model} object.
#' @param fim A \code{Fim} object.
#' @param ... \code{designName} and \code{plotOptions}.
#' @name processArmEvaluationResults
#' @export
processArmEvaluationResults = new_generic( "processArmEvaluationResults", c( "arm", "model", "fim" ) )

#' Sensitivity-index plots for one arm
#' @param arm An \code{Arm} object.
#' @param model A \code{Model} object.
#' @param fim A \code{Fim} object.
#' @param ... \code{designName} and \code{plotOptions}.
#' @name processArmEvaluationSI
#' @export
processArmEvaluationSI = new_generic( "processArmEvaluationSI", c( "arm", "model", "fim" ) )

#' ggplot of predicted responses with sampling points
#' @param arm An \code{Arm} object.
#' @param ... Evaluation outputs and plot options for the \code{Arm} method.
#' @name plotEvaluationResults
#' @export
plotEvaluationResults = new_generic( "plotEvaluationResults", c( "arm" ) )

#' ggplot of parameter sensitivities over time
#' @param arm An \code{Arm} object.
#' @param ... Gradient data and plot options for the \code{Arm} method.
#' @name plotEvaluationSI
#' @export
plotEvaluationSI = new_generic( "plotEvaluationSI", c( "arm" ) )

#' Dose and sampling summary for reports
#' @param arm An \code{Arm} object.
#' @param ... Not used.
#' @name getArmData
#' @export
getArmData = new_generic( "getArmData", c( "arm" ) )

#' Administration and sampling constraints for an optimizer
#' @param arm An \code{Arm} object.
#' @param optimizationAlgorithm An optimization algorithm object.
#' @param ... Not used.
#' @name getArmConstraints
#' @export
getArmConstraints = new_generic( "getArmConstraints", c( "arm", "optimizationAlgorithm" ) )

#' Tabular administration settings for reports
#' @param arm An \code{Arm} object.
#' @param ... Not used.
#' @name armAdministration
#' @export
armAdministration = new_generic( "armAdministration", c( "arm" ) )

method( evaluateArm, Arm ) = function( arm, model, fim ) {
  model = .pfimPrepareModelForEvaluation( model, arm )
  prop( arm, "evaluationModel" )     = evaluateModel( model, arm )
  prop( arm, "evaluationGradients" ) = evaluateModelGradient( model, arm )
  prop( arm, "evaluationVariance" )  = evaluateModelVariance( model, arm )
  armFim = .duplicateFim( fim )
  prop( arm, "evaluationFim" ) = evaluateFim( armFim, model, arm )
  arm
}

method( armAdministration, Arm ) = function( arm ) {
  armName   = prop( arm, "name" )
  armSize   = round( prop( arm, "size" ), 2 )
  adminList = prop( arm, "administrations" )
  map( adminList, function( adm ) {
    list(
      "Design name"        = "Design optimized",
      "Arms name"          = armName,
      "Number of subjects" = as.character( armSize ),
      "Outcome"            = prop( adm, "outcome" ),
      "Dose"               = as.character( prop( adm, "dose" ) ),
      "Time of dose"       = if ( length( prop( adm, "timeDose" ) ) > 0L )
        as.character( prop( adm, "timeDose" ) ) else ".",
      "tau"                = as.character( prop( adm, "tau" ) ),
      "Tinf"               = if ( length( prop( adm, "Tinf" ) ) > 0L )
        as.character( prop( adm, "Tinf" ) ) else "."
    )
  } )
}

method( getArmData, Arm ) = function( arm ) {
  armName         = prop( arm, "name" )
  armSize         = round( prop( arm, "size" ), 2 )
  administrations = prop( arm, "administrations" )
  doseList = map( administrations, function( adm ) {
    list( outcome = prop( adm, "outcome" ), dose = prop( adm, "dose" ) )
  } )
  doseDict = set_names(
    map( doseList, ~ paste( .x$dose, collapse = ", " ) ),
    map_chr( doseList, ~ .x$outcome )
  )
  samplingList     = prop( arm, "samplingTimes" )
  samplingOutcomes = map_chr( samplingList, ~ prop( .x, "outcome" ) )
  samplingTimes    = map( samplingList, ~ prop( .x, "samplings" ) )
  map2( samplingOutcomes, samplingTimes, function( outc, samps ) {
    doseVal = if ( outc %in% names( doseDict ) ) doseDict[[ outc ]] else "."
    list(
      "Arms name"          = armName,
      "Number of subjects" = as.character( armSize ),
      "Outcome"            = outc,
      "Dose"               = doseVal,
      "Sampling times"     = paste0( "(", paste( round( samps, 2 ), collapse = ", " ), ")" )
    )
  } )
}

method( getSamplingData, Arm ) = function( arm ) {
  samplingTimes = prop( arm, "samplingTimes" )
  samplings = map( samplingTimes, ~ prop( .x, "samplings" ) ) |>
    set_names( map_chr( samplingTimes, ~ prop( .x, "outcome" ) ) )
  list(
    samplingTimes = samplingTimes,
    samplings     = samplings,
    samplingMax   = samplings |> flatten_dbl() |> max()
  )
}

method( updateSamplingTimes, Arm ) = function( arm, samplingData ) {
  prop( arm, "samplingTimes" ) = map( samplingData$samplingTimes, function( st ) {
    prop( st, "samplings" ) = sort( unique( c(
      samplingData$samplings |> flatten_dbl(),
      seq( 0.0, samplingData$samplingMax, 0.1 )
    ) ) )
    st
  } )
  arm
}

.armConstraintsDiscrete = function( arm ) {
  armName = prop( arm, "name" )
  armSize = prop( arm, "size" )
  admins  = map( prop( arm, "administrationsConstraints" ), function( ac ) {
    outcome = prop( ac, "outcome" )
    doses   = paste0( "(", paste( unlist( prop( ac, "doses" ) ), collapse = ", " ), ")" )
    set_names( list( doses ), outcome )
  } ) |> flatten()
  map( prop( arm, "samplingTimesConstraints" ), function( sc ) {
    outcome = prop( sc, "outcome" )
    doseConstraints = admins[[ outcome ]]
    if ( is.null( doseConstraints ) ) doseConstraints = "."
    list(
      "Arms name"                       = armName,
      "Number of subjects"              = armSize,
      "Outcome"                         = outcome,
      "Initial samplings"               = paste0( "(", paste( prop( sc, "initialSamplings" ), collapse = ", " ), ")" ),
      "Fixed times"                     = paste0( "(", paste( prop( sc, "fixedTimes" ), collapse = ", " ), ")" ),
      "Number of samplings optimisable" = as.character( prop( sc, "numberOfsamplingsOptimisable" ) ),
      "Dose constraints"                = doseConstraints
    )
  } )
}

method( getArmConstraints, list( Arm, MultiplicativeAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsDiscrete( arm )
}

method( getArmConstraints, list( Arm, FedorovWynnAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsDiscrete( arm )
}

method( getArmConstraints, list( Arm, SimplexAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsContinuous( arm )
}

method( getArmConstraints, list( Arm, PSOAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsContinuous( arm )
}

method( getArmConstraints, list( Arm, PGBOAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsContinuous( arm )
}

#' @keywords internal
.normalizePlotOptions = function( plotOptions, outputNames ) {
  unitTime = plotOptions$unitTime
  if ( is.null( unitTime ) || length( unitTime ) == 0L )
    unitTime = " "
  else if ( length( unitTime ) > 1L )
    unitTime = unitTime[[ 1L ]]
  out = unlist( outputNames, use.names = FALSE )
  if ( length( out ) == 0L ) out = names( outputNames )
  raw = plotOptions$unitOutcomes
  if ( is.null( raw ) || length( raw ) == 0L ) {
    unitOutcomes = stats::setNames( rep( " ", length( out ) ), out )
  } else {
    raw = as.character( unlist( raw, use.names = FALSE ) )
    if ( length( raw ) == 1L && length( out ) > 1L )
      raw = rep( raw, length( out ) )
    unitOutcomes = stats::setNames( raw[ seq_along( out ) ], out )
  }
  list( unitTime = unitTime, unitOutcomes = unitOutcomes )
}

method( processArmEvaluationResults, list( Arm, Model, Fim ) ) = function(
    arm, model, fim, designName, plotOptions ) {
  outputNames = as.list( prop( model, "outputNames" ) )
  samplingData = getSamplingData( arm )
  arm = updateSamplingTimes( arm, samplingData )
  model = defineModelAdministration( model, arm )
  evaluationModel = evaluateModel( model, arm )
  plotEvaluationResults( arm, evaluationModel, outputNames, samplingData, designName, plotOptions )
}

method( processArmEvaluationSI, list( Arm, Model, Fim ) ) = function(
    arm, model, fim, designName, plotOptions ) {
  outputNames = as.list( prop( model, "outputNames" ) )
  samplingData = getSamplingData( arm )
  arm = updateSamplingTimes( arm, samplingData )
  model = defineModelAdministration( model, arm )
  parametersNames = prop( model, "modelParameters" ) |> map_chr( ~ prop( .x, "name" ) )
  evaluationModelGradient = evaluateModelGradient( model, arm )
  timeSeq = seq( from = 0, by = 0.1, length.out = nrow( pluck( evaluationModelGradient, 1 ) ) )
  evaluationModelGradient = map2(
    outputNames, evaluationModelGradient,
    function( outputName, gradient ) data.frame( time = timeSeq, gradient )
  ) |> set_names( outputNames )
  plotEvaluationSI(
    arm, evaluationModelGradient, parametersNames,
    outputNames, samplingData, designName, plotOptions
  )
}

method( plotEvaluationResults, Arm ) = function(
    arm, evaluationModel, outputNames, samplingData, designName, plotOptions ) {
  units     = .normalizePlotOptions( plotOptions, outputNames )
  unitXAxis = units$unitTime
  unitYAxis = units$unitOutcomes
  armName   = prop( arm, "name" )
  plotList  = list()
  plotList[[ designName ]] = list()
  plotList[[ designName ]][[ armName ]] = map2(
    outputNames, samplingData$samplings,
    function( outputName, sampling ) {
      data = evaluationModel[[ outputName ]]
      samplingPoints = data[ data$time %in% sampling, ]
      ggplot( data, aes( x = time, y = .data[[ outputName ]] ) ) +
        geom_line() +
        geom_point( data = samplingPoints, aes( x = time, y = .data[[ outputName ]] ), color = "red" ) +
        labs(
          x = paste0( "Time (", unitXAxis, ")\n\nDesign: ", sub( "_", " ", designName ), "      Arm: ", armName ),
          y = paste0( outputName, " (", unitYAxis[[ outputName ]], ")\n" )
        ) +
        scale_x_continuous(
          breaks = pretty_breaks( n = 10 ),
          sec.axis = sec_axis( ~ . * 1, breaks = round( sampling, 2 ), name = "Sampling times" )
        ) +
        scale_y_continuous( breaks = pretty_breaks( n = 10 ) ) +
        theme(
          legend.position = "none",
          axis.title.x.top = element_text( color = "red", vjust = 2.0 ),
          axis.text.x.top  = element_text( angle = 90, hjust = 0, color = "red" ),
          plot.title = element_text( size = 16, hjust = 0.5 ),
          axis.title.x = element_text( size = 16 ),
          axis.title.y = element_text( size = 16 ),
          axis.text.x  = element_text( size = 16, angle = 90, vjust = 0.5 ),
          axis.text.y  = element_text( size = 16 ),
          strip.text.x = element_text( size = 16 )
        )
    }
  ) |> set_names( outputNames )
  plotList
}

method( plotEvaluationSI, Arm ) = function(
    arm, evaluationModelGradient, parametersNames, outputNames,
    samplingData, designName, plotOptions ) {
  unitXAxis = .normalizePlotOptions( plotOptions, outputNames )$unitTime
  armName   = prop( arm, "name" )
  plotList  = list()
  plotList[[ designName ]] = list()
  plotList[[ designName ]][[ armName ]] = map2(
    outputNames, samplingData$samplings,
    function( outputName, sampling ) {
      gradientData = evaluationModelGradient[[ outputName ]]
      minYAxis = min( gradientData[ , parametersNames ], na.rm = TRUE )
      maxYAxis = max( gradientData[ , parametersNames ], na.rm = TRUE )
      map( parametersNames, function( parameterName ) {
        data = as_tibble( gradientData[ , c( "time", parameterName ) ] )
        names( data )[ 2L ] = "parameterValue"
        samplingPoints = data[ data$time %in% sampling, ]
        ggplot( data, aes( x = time, y = parameterValue ) ) +
          geom_line() +
          geom_point( data = samplingPoints, color = "red" ) +
          labs(
            y = paste0( "df/d", parameterName ),
            x = paste0(
              "Time (", unitXAxis, ")\n\n",
              "Design: ", gsub( "_", " ", designName ), "   ",
              "Arm: ", armName, "   ",
              "Output: ", outputName, "   ",
              "Parameter: ", parameterName
            )
          ) +
          scale_x_continuous(
            breaks = pretty_breaks( n = 10 ),
            sec.axis = sec_axis( ~ ., breaks = round( sampling, 2 ), name = "Sampling times" )
          ) +
          scale_y_continuous( breaks = pretty_breaks( n = 10 ), limits = c( minYAxis, maxYAxis ) ) +
          theme(
            legend.position = "none",
            axis.title.x.top = element_text( color = "red", vjust = 2.0 ),
            axis.text.x.top  = element_text( angle = 90, hjust = 0, color = "red" ),
            plot.title = element_text( size = 16, hjust = 0.5 ),
            axis.title.x = element_text( size = 16 ),
            axis.title.y = element_text( size = 16 ),
            axis.text.x  = element_text( size = 16, angle = 90, vjust = 0.5 ),
            axis.text.y  = element_text( size = 16 ),
            strip.text.x = element_text( size = 16 )
          )
      } ) |> set_names( parametersNames )
    }
  ) |> set_names( outputNames )
  plotList
}
