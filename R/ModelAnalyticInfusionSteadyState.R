#' @title ModelAnalyticInfusionSteadyState
#' @description Analytic infusion model at steady state.
#' @inheritParams ModelInfusion
#' @param wrapperModelAnalyticInfusion                   Wrapper for the analytic solver.
#' @param functionArgumentsModelAnalyticInfusion         A list with the function arguments.
#' @param functionArgumentsSymbolModelAnalyticInfusion   A list with the function argument symbols.
#' @param solverInputs                                   A list with the solver inputs.
#' @inheritParams Model
#' @include ModelInfusion.R
#' @include ModelAnalytic.R
#' @export

ModelAnalyticInfusionSteadyState = new_class( "ModelAnalyticInfusionSteadyState",
                                              package = "PFIM",
                                              parent  = ModelInfusion,
                                              properties = list(
                                                wrapperModelAnalyticInfusion                 = new_property(class_list, default = list()),
                                                functionArgumentsModelAnalyticInfusion       = new_property(class_list, default = list()),
                                                functionArgumentsSymbolModelAnalyticInfusion = new_property(class_list, default = list()),
                                                solverInputs                                 = new_property(class_list, default = list())
                                              ))


method( defineModelWrapper, ModelAnalyticInfusionSteadyState ) = function( model, evaluation ) {

  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )
  libraryOutcomesWithAdmin   = .administeredLibraryOutcomeNames( evaluation )

  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNames      = paste0( "t_",    outcomesWithAdministration )
  TinfNames      = paste0( "Tinf_", outcomesWithAdministration )

  equations                  = prop( evaluation, "modelEquations" )
  equationsDuringWithAdmin   = equations$duringInfusion[  names( equations$duringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithAdmin    = equations$afterInfusion[   names( equations$afterInfusion  ) %in% libraryOutcomesWithAdmin ]
  equationsDuringWithNoAdmin = equations$duringInfusion[ !names( equations$duringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithNoAdmin  = equations$afterInfusion[  !names( equations$afterInfusion  ) %in% libraryOutcomesWithAdmin ]

  outputsForEvaluation = prop( evaluation, "outputs" )
  outputAdmin   = unlist( outputsForEvaluation[[1L]] )
  outputNoAdmin = if ( length( outputsForEvaluation ) >= 2L ) unlist( outputsForEvaluation[[2L]] ) else character(0L)

  # Steady-state adds "tau" to the full shared argument list.
  functionArguments = unique( c( doseNames, TinfNames, outcomesWithAdministration,
                                 parameterNames, timeNames, "tau" ) )

  prop( model, "wrapperModelAnalyticInfusion" ) = list(
    functionDefinitionDuringInfusionWithAdmin   = .buildAnalyticWrapper(
      equationsDuringWithAdmin, functionArguments,
      .libraryEquationTimeNames( evaluation, names( equationsDuringWithAdmin ) ), outputAdmin ),
    functionDefinitionDuringInfusionWithNoAdmin = .buildAnalyticWrapper(
      equationsDuringWithNoAdmin, functionArguments,
      .libraryEquationTimeNames( evaluation, names( equationsDuringWithNoAdmin ) ), outputNoAdmin ),
    functionDefinitionAfterInfusionWithAdmin    = .buildAnalyticWrapper(
      equationsAfterWithAdmin, functionArguments,
      .libraryEquationTimeNames( evaluation, names( equationsAfterWithAdmin ) ), outputAdmin ),
    functionDefinitionAfterInfusionWithNoAdmin  = .buildAnalyticWrapper(
      equationsAfterWithNoAdmin, functionArguments,
      .libraryEquationTimeNames( evaluation, names( equationsAfterWithNoAdmin ) ), outputNoAdmin )
  )
  prop( model, "functionArgumentsModelAnalyticInfusion" ) = list( functionArguments = functionArguments )
  prop( model, "functionArgumentsSymbolModelAnalyticInfusion" ) = list(
    functionArgumentsSymbol = map( functionArguments, as.symbol )
  )
  prop( model, "outputNames" ) = unlist( names( equations$duringInfusion ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


method( defineModelAdministration, ModelAnalyticInfusionSteadyState ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  samplingTimes              = prop( arm,   "samplingTimes" )
  samplings                  = map( samplingTimes, ~ prop( .x, "samplings" ) ) |>
    unlist() |> sort() |> unique()
  maxSampling = max( samplings )

  solverInputs = map( administrations, function( adm ) {
    tau      = prop( adm, "tau" )
    dosing   = .alignAdministrationDosing( adm )
    timeDose = dosing$timeDose
    dose     = dosing$dose
    Tinf     = dosing$Tinf

    if ( tau != 0 ) {
      timeDose = seq( 0, maxSampling, tau )
      dose     = rep( dose, length( timeDose ) )
      Tinf     = rep( Tinf, length( timeDose ) )
    }

    # outcomes do not accumulate each other's infusion-window labels.
    duringAndAfter = rep( "afterInfusion", length( samplings ) )

    Tinfs = map2( timeDose, timeDose + Tinf, c )
    samplingsDuringInfusion = map( Tinfs, function( iv ) {
      samplings |> keep( ~ .x >= min(iv) & .x < max(iv) )
    }) |> unlist() |> unique()
    duringAndAfter[ samplings %in% samplingsDuringInfusion ] = "duringInfusion"

    samplingTimeDoses = timeDose |> map( ~ ifelse( samplings - .x > 0, samplings - .x, 0 ) )
    indicesDoses = map_int( samplings, function( s ) {
      idx = which( s >= timeDose ); idx[ length(idx) ]
    })

    data = data.frame( duringAndAfter, indicesDoses, samplings, samplingTimeDoses )
    colnames( data ) = c( "duringAndAfter", "indicesDoses", "samplings",
                          paste0( "samplingTimeDoses", seq_along( dose ) ) )
    list( data = data, dose = dose, Tinf = Tinf, tau = tau )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}


method( evaluateModel, ModelAnalyticInfusionSteadyState ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = map_chr( administrations, ~ prop( .x, "outcome" ) )
  outputNames                = prop( model, "outputNames" ) |> unlist()
  solverInputs               = prop( model, "solverInputs" )
  samplings                  = prop( model, "samplings" )

  wrappers        = prop( model, "wrapperModelAnalyticInfusion" )
  fnDuringAdmin   = wrappers$functionDefinitionDuringInfusionWithAdmin
  fnDuringNoAdmin = wrappers$functionDefinitionDuringInfusionWithNoAdmin
  fnAfterAdmin    = wrappers$functionDefinitionAfterInfusionWithAdmin
  fnAfterNoAdmin  = wrappers$functionDefinitionAfterInfusionWithNoAdmin

  functionArguments = prop( model, "functionArgumentsModelAnalyticInfusion" )$functionArguments
  functionArgumentsSymbol = prop( model, "functionArgumentsSymbolModelAnalyticInfusion" )$functionArgumentsSymbol
  argsTemplate = stats::setNames( functionArgumentsSymbol, functionArguments )

  mu = .extractMu( prop( model, "modelParameters" ) )
  list2env( as.list( mu ), envir = environment() )

  evaluationModelTmpList = map( seq_along( samplings ), function( iterTime ) {

    outcomesResults = map( outcomesWithAdministration, function( outcome ) {
      data         = solverInputs[[ outcome ]]$data
      tau          = solverInputs[[ outcome ]]$tau
      duringAfter  = data$duringAndAfter[ iterTime ]
      indicesDoses = data$indicesDoses[ iterTime ]
      sampCols     = colnames( data )[ str_detect( colnames( data ), "samplingTimeDoses" ) ]
      sampTimes    = as.numeric( data[ iterTime, sampCols ] ) |> unlist() |> unname()

      evalAdmin = if ( duringAfter == "duringInfusion" ) {
        if ( indicesDoses == 1L ) {
          assign( paste0( "t_", outcome ), sampTimes[ indicesDoses ], envir = environment() )
          assign( paste0( "dose_", outcome ), solverInputs[[ outcome ]]$dose[ indicesDoses ], envir = environment() )
          assign( paste0( "Tinf_", outcome ), solverInputs[[ outcome ]]$Tinf[ indicesDoses ], envir = environment() )
          do.call( fnDuringAdmin, argsTemplate ) |> unlist()
        } else {
          sampTimesUsed = sampTimes[ seq_len( indicesDoses ) ]
          samplingDuring = tail( sampTimesUsed, 1L )
          samplingAfter  = sampTimesUsed[ seq_len( indicesDoses - 1L ) ]

          assign( paste0( "t_", outcome ), samplingDuring, envir = environment() )
          assign( paste0( "dose_", outcome ), solverInputs[[ outcome ]]$dose[ indicesDoses ], envir = environment() )
          assign( paste0( "Tinf_", outcome ), solverInputs[[ outcome ]]$Tinf[ indicesDoses ], envir = environment() )

          evalDuring = do.call( fnDuringAdmin, argsTemplate ) |> unlist()
          evalDuring + sum( map_dbl( seq_len( indicesDoses - 1L ), function( idx ) {
            assign( paste0( "t_", outcome ), samplingAfter[ idx ], envir = environment() )
            assign( paste0( "dose_", outcome ), solverInputs[[ outcome ]]$dose[ idx ], envir = environment() )
            assign( paste0( "Tinf_", outcome ), solverInputs[[ outcome ]]$Tinf[ idx ], envir = environment() )
            do.call( fnAfterAdmin, argsTemplate ) |> unlist()
          } ) )
        }
      } else if ( indicesDoses == 1L ) {
        assign( paste0( "t_", outcome ), sampTimes[ indicesDoses ], envir = environment() )
        assign( paste0( "dose_", outcome ), solverInputs[[ outcome ]]$dose[ indicesDoses ], envir = environment() )
        assign( paste0( "Tinf_", outcome ), solverInputs[[ outcome ]]$Tinf[ indicesDoses ], envir = environment() )
        do.call( fnAfterAdmin, argsTemplate ) |> unlist()
      } else {
        sampTimesUsed = sampTimes[ seq_len( indicesDoses ) ]
        samplingDuring = tail( sampTimesUsed, 1L )
        samplingAfter  = sampTimesUsed[ seq_len( indicesDoses - 1L ) ]

        assign( paste0( "t_", outcome ), samplingDuring, envir = environment() )
        assign( paste0( "dose_", outcome ), solverInputs[[ outcome ]]$dose[ indicesDoses ], envir = environment() )
        assign( paste0( "Tinf_", outcome ), solverInputs[[ outcome ]]$Tinf[ indicesDoses ], envir = environment() )

        evalAfter = do.call( fnAfterAdmin, argsTemplate ) |> unlist()
        evalAfter + sum( map_dbl( seq_len( indicesDoses - 1L ), function( idx ) {
          assign( paste0( "t_", outcome ), samplingAfter[ idx ], envir = environment() )
          assign( paste0( "dose_", outcome ), solverInputs[[ outcome ]]$dose[ idx ], envir = environment() )
          assign( paste0( "Tinf_", outcome ), solverInputs[[ outcome ]]$Tinf[ idx ], envir = environment() )
          do.call( fnAfterAdmin, argsTemplate ) |> unlist()
        } ) )
      }

      assign( outcome, evalAdmin, envir = environment() )
      evalNoAdmin = do.call( fnAfterNoAdmin, argsTemplate ) |> unlist()

      if ( is.null( evalNoAdmin ) || length( evalNoAdmin ) == 0L )
        data.frame( evalAdmin )
      else
        data.frame( evalAdmin, evalNoAdmin )
    })

    data.frame( time = samplings[[ iterTime ]], do.call( cbind, outcomesResults ) )
  })

  evaluationModelTmp = list_rbind( evaluationModelTmpList )
  colnames( evaluationModelTmp ) = c( "time", outputNames )

  samplings_by_output = set_names(
    map( prop( arm, "samplingTimes" ), ~ prop( .x, "samplings" ) ),
    outputNames
  )
  set_names(
    map( outputNames, ~ evaluationModelTmp[ evaluationModelTmp$time %in% samplings_by_output[[ .x ]], c( "time", .x ) ] ),
    outputNames
  )
}


method( definePKModel, list( ModelAnalyticInfusionSteadyState, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}
