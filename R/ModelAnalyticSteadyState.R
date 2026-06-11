#' ModelAnalyticSteadyState
#' @title ModelAnalyticSteadyState
#' @description Analytic model at steady state (tau in equations).
#' @inheritParams ModelAnalytic
#' @param wrapperModelAnalytic                   Wrapper for the analytic solver.
#' @param functionArgumentsModelAnalytic         A list with the function arguments.
#' @param functionArgumentsSymbolModelAnalytic   A list with the function argument symbols.
#' @param solverInputs                           A list with the solver inputs.
#' @inheritParams Model
#' @include ModelAnalytic.R
#' @export

ModelAnalyticSteadyState = new_class( "ModelAnalyticSteadyState",
                                      package = "PFIM",
                                      parent  = ModelAnalytic,
                                      properties = list(
                                        wrapperModelAnalytic                 = new_property(class_list, default = list()),
                                        functionArgumentsModelAnalytic       = new_property(class_list, default = list()),
                                        functionArgumentsSymbolModelAnalytic = new_property(class_list, default = list()),
                                        solverInputs                         = new_property(class_list, default = list())
                                      ))


method( defineModelWrapper, ModelAnalyticSteadyState ) = function( model, evaluation ) {

  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )

  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNames      = paste0( "t_",    outcomesWithAdministration )

  equations            = prop( evaluation, "modelEquations" )
  equationsWithAdmin   = equations[  names( equations ) %in% outcomesWithAdministration ]
  equationsWithNoAdmin = equations[ !names( equations ) %in% outcomesWithAdministration ]

  outputsForEvaluation = prop( evaluation, "outputs" )
  outputAdmin   = unlist( outputsForEvaluation[[1L]] )
  outputNoAdmin = if ( length( outputsForEvaluation ) >= 2L ) unlist( outputsForEvaluation[[2L]] ) else character(0L)

  # Steady-state adds "tau" to the function arguments
  functionArgumentsWithAdmin   = unique( c( doseNames, parameterNames, timeNames, "tau" ) )
  functionArgumentsWithNoAdmin = unique( c( outcomesWithAdministration, parameterNames, timeNames, "tau" ) )

  prop( model, "wrapperModelAnalytic" ) = list(
    functionDefinitionWithAdmin   = .buildAnalyticWrapper( equationsWithAdmin,   functionArgumentsWithAdmin,   timeNames, outputAdmin   ),
    functionDefinitionWithNoAdmin = .buildAnalyticWrapper( equationsWithNoAdmin, functionArgumentsWithNoAdmin, timeNames, outputNoAdmin )
  )
  prop( model, "functionArgumentsModelAnalytic" ) = list(
    functionArgumentsWithAdmin   = functionArgumentsWithAdmin,
    functionArgumentsWithNoAdmin = functionArgumentsWithNoAdmin
  )
  prop( model, "outputNames" ) = unlist( names( equations ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


method( defineModelAdministration, ModelAnalyticSteadyState ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  samplingTimes              = prop( arm,   "samplingTimes" )
  samplings                  = map( samplingTimes, ~ prop( .x, "samplings" ) ) |>
    unlist() |> sort() |> unique()

  solverInputs = map( administrations, function( adm ) {
    tau         = prop( adm, "tau" )
    dosing      = .alignAdministrationDosing( adm )
    timeDose    = dosing$timeDose
    dose        = dosing$dose
    maxSampling = max( samplings )

    if ( tau != 0 ) {
      timeDose = seq( 0, maxSampling, tau )
      dose     = rep( dose, length( timeDose ) )
    }

    # Relative sampling times from each dose event: positive delay or full time.
    timeDose = timeDose |>
      map( ~ ifelse( samplings - .x > 0, samplings - .x, samplings ) ) |>
      reduce( cbind )

    indicesDoses = map_int( seq_len( nrow( timeDose ) ),
                            ~ length( unique( timeDose[.x, ] ) ) )

    list( data = data.frame( timeDose, indicesDoses ), dose = dose, tau = tau )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}


method( evaluateModel, ModelAnalyticSteadyState ) = function( model, arm ) {

  parameters                 = prop( model, "modelParameters" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  outputNames                = prop( model, "outputNames" )
  samplings                  = prop( model, "samplings" )
  solverInputs               = prop( model, "solverInputs" )

  wrapperModelAnalytic          = prop( model, "wrapperModelAnalytic" )
  functionDefinitionWithAdmin   = wrapperModelAnalytic$functionDefinitionWithAdmin
  functionDefinitionWithNoAdmin = wrapperModelAnalytic$functionDefinitionWithNoAdmin

  mu = .extractMu( parameters )

  evaluationModelTmpList = map( seq_along( samplings ), function( iterTime ) {

    # Initialise argument list with parameter values at each time step
    currentArgsNoAdmin = as.list( mu )

    outcomesResults = vector( "list", length( outcomesWithAdministration ) )
    for ( i in seq_along( outcomesWithAdministration ) ) {
      outcome = outcomesWithAdministration[[ i ]]
      data         = solverInputs[[ outcome ]]$data
      dose         = solverInputs[[ outcome ]]$dose
      indicesDoses = data$indicesDoses[ iterTime ]

      timesVec = as.numeric( data[ iterTime, seq_len( indicesDoses ) ] )
      dosesVec = as.numeric( dose[ seq_len( indicesDoses ) ] )

      argsAdmin = c( currentArgsNoAdmin,
                     set_names( list( timesVec, dosesVec ),
                                c( paste0( "t_", outcome ), paste0( "dose_", outcome ) ) ) )
      argsAdmin[[ "tau" ]] = solverInputs[[ outcome ]]$tau

      evaluationOutcomeWithAdmin = sum( do.call( functionDefinitionWithAdmin, argsAdmin )[[ 1L ]] )
      currentArgsNoAdmin[[ outcome ]] = evaluationOutcomeWithAdmin
      currentArgsNoAdmin[[ "tau" ]]   = solverInputs[[ outcome ]]$tau

      evaluationOutcomeWithNoAdmin = do.call( functionDefinitionWithNoAdmin, currentArgsNoAdmin )[[ 1L ]]

      outcomesResults[[ i ]] = if ( is.null( evaluationOutcomeWithNoAdmin ) ||
                                    length( evaluationOutcomeWithNoAdmin ) == 0L ) {
        data.frame( Admin = evaluationOutcomeWithAdmin )
      } else {
        data.frame( Admin = evaluationOutcomeWithAdmin, NoAdmin = evaluationOutcomeWithNoAdmin )
      }
    }
    data.frame( time = samplings[[iterTime]], do.call( cbind, outcomesResults ) )
  })

  evaluationModelTmp = list_rbind( evaluationModelTmpList )
  colnames( evaluationModelTmp ) = c( "time", outputNames )

  samplings_by_output = set_names( map( prop( arm, "samplingTimes" ), ~ prop( .x, "samplings" ) ), outputNames )
  set_names(
    map( outputNames, ~ evaluationModelTmp[ evaluationModelTmp$time %in% samplings_by_output[[.x]], c("time", .x) ] ),
    outputNames
  )
}


method( definePKModel, list( ModelAnalyticSteadyState, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}


method( definePKPDModel, list( ModelAnalyticSteadyState, ModelAnalytic, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {
    c( prop( pkModel, "modelEquations" ), prop( pdModel, "modelEquations" ) )
  }


method( definePKPDModel, list( ModelAnalyticSteadyState, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    equations = c(
      convertPKModelAnalyticToPKModelODE( pkModel ),
      prop( pdModel, "modelEquations" )
    )
    eq = remapPkpdLibraryEquations( equations, pfimproject )
    set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ) ) )
  }
