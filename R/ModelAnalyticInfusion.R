#' @title ModelAnalyticInfusion
#' @description Closed-form model with infusion (dose in equations).
#' @inheritParams ModelInfusion
#' @param wrapperModelAnalyticInfusion                   Wrapper for the analytic solver.
#' @param functionArgumentsModelAnalyticInfusion         A list with the function arguments.
#' @param functionArgumentsSymbolModelAnalyticInfusion   A list with the function argument symbols.
#' @param solverInputs                                   A list with the solver inputs.
#' @inheritParams Model
#' @include ModelInfusion.R
#' @include ModelAnalytic.R
#' @export

ModelAnalyticInfusion = new_class(
  "ModelAnalyticInfusion",
  package = "PFIM",
  parent  = ModelInfusion,
  properties = list(
    wrapperModelAnalyticInfusion                 = new_property(class_list, default = list()),
    functionArgumentsModelAnalyticInfusion       = new_property(class_list, default = list()),
    functionArgumentsSymbolModelAnalyticInfusion = new_property(class_list, default = list()),
    solverInputs                                 = new_property(class_list, default = list())
  ))


method( defineModelWrapper, ModelAnalyticInfusion ) = function( model, evaluation ) {

  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )
  libraryOutcomesWithAdmin   = .administeredLibraryOutcomeNames( evaluation )

  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNames      = paste0( "t_",    outcomesWithAdministration )
  TinfNames      = paste0( "Tinf_", outcomesWithAdministration )

  equations               = prop( evaluation, "modelEquations" )
  equationsDuringInfusion = equations$duringInfusion
  equationsAfterInfusion  = equations$afterInfusion

  equationsDuringWithAdmin   = equationsDuringInfusion[  names( equationsDuringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithAdmin    = equationsAfterInfusion[   names( equationsAfterInfusion  ) %in% libraryOutcomesWithAdmin ]
  equationsDuringWithNoAdmin = equationsDuringInfusion[ !names( equationsDuringInfusion ) %in% libraryOutcomesWithAdmin ]
  equationsAfterWithNoAdmin  = equationsAfterInfusion[  !names( equationsAfterInfusion  ) %in% libraryOutcomesWithAdmin ]

  outputsForEvaluation = prop( evaluation, "outputs" )
  outputAdmin   = unlist( outputsForEvaluation[[1L]] )
  outputNoAdmin = if ( length( outputsForEvaluation ) >= 2L ) unlist( outputsForEvaluation[[2L]] ) else character(0L)

  # All four wrappers share the same full argument signature.
  functionArguments = unique( c( doseNames, TinfNames, outcomesWithAdministration, parameterNames, timeNames ) )

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
  prop( model, "outputNames" ) = unlist( names( equationsDuringInfusion ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


method( defineModelAdministration, ModelAnalyticInfusion ) = function( model, arm ) {

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
    list( data = data, dose = dose, Tinf = Tinf )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}


method( evaluateModel, ModelAnalyticInfusion ) = function( model, arm ) {

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

  mu = .extractMu( prop( model, "modelParameters" ) )

  evaluationModelTmpList = map( seq_along( samplings ), function( iterTime ) {

    # Initialise argument list with parameter values; will be extended with PK values
    argsBase = as.list( mu )

    outcomesResults = vector( "list", length( outcomesWithAdministration ) )
    for ( i in seq_along( outcomesWithAdministration ) ) {
      outcome = outcomesWithAdministration[[ i ]]
      data         = solverInputs[[ outcome ]]$data
      duringAfter  = data$duringAndAfter[ iterTime ]
      indicesDoses = data$indicesDoses[ iterTime ]
      sampCols     = colnames( data )[ str_detect( colnames( data ), "samplingTimeDoses" ) ]
      sampTimes    = as.numeric( data[ iterTime, sampCols ] )

      tName    = paste0( "t_", outcome )
      doseName = paste0( "dose_", outcome )
      TinfName = paste0( "Tinf_", outcome )

      argsCall = argsBase

      evalAdmin = if ( duringAfter == "duringInfusion" ) {
        if ( indicesDoses == 1L ) {
          argsCall[[ tName ]]    = sampTimes[1L]
          argsCall[[ doseName ]] = solverInputs[[ outcome ]]$dose[1L]
          argsCall[[ TinfName ]] = solverInputs[[ outcome ]]$Tinf[1L]
          do.call( fnDuringAdmin, argsCall )[[ 1L ]]
        } else {
          idxBefore = seq_len( indicesDoses - 1L )
          argsCall[[ tName ]]    = sampTimes[ idxBefore ]
          argsCall[[ doseName ]] = solverInputs[[ outcome ]]$dose[ idxBefore ]
          argsCall[[ TinfName ]] = solverInputs[[ outcome ]]$Tinf[ idxBefore ]
          sumAfter = sum( do.call( fnAfterAdmin, argsCall )[[ 1L ]] )

          argsCall[[ tName ]]    = sampTimes[ indicesDoses ]
          argsCall[[ doseName ]] = solverInputs[[ outcome ]]$dose[ indicesDoses ]
          argsCall[[ TinfName ]] = solverInputs[[ outcome ]]$Tinf[ indicesDoses ]
          sumAfter + do.call( fnDuringAdmin, argsCall )[[ 1L ]]
        }
      } else {
        idxAll = seq_len( indicesDoses )
        argsCall[[ tName ]]    = sampTimes[ idxAll ]
        argsCall[[ doseName ]] = solverInputs[[ outcome ]]$dose[ idxAll ]
        argsCall[[ TinfName ]] = solverInputs[[ outcome ]]$Tinf[ idxAll ]
        sum( do.call( fnAfterAdmin, argsCall )[[ 1L ]] )
      }

      argsBase[[ outcome ]] = evalAdmin
      evalNoAdmin = do.call( fnAfterNoAdmin, argsBase )[[ 1L ]]

      outcomesResults[[ i ]] = if ( is.null( evalNoAdmin ) || length( evalNoAdmin ) == 0L ) {
        data.frame( evalAdmin )
      } else {
        data.frame( evalAdmin, evalNoAdmin )
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


method( definePKModel, list( ModelAnalyticInfusion, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}


method( definePKPDModel, list( ModelAnalyticInfusion, ModelAnalytic, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {
    pkEq = prop( pkModel, "modelEquations" )
    pdEq = prop( pdModel, "modelEquations" )
    list(
      duringInfusion = c( pkEq$duringInfusion, pdEq ),
      afterInfusion  = c( pkEq$afterInfusion,  pdEq )
    )
  }


method( definePKPDModel, list( ModelAnalyticInfusion, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    pkEqODE = remapPkpdLibraryEquations(
      convertPKModelAnalyticToPKModelODE( pkModel ),
      pfimproject
    )
    pdEq = remapPkpdLibraryEquations(
      prop( pdModel, "modelEquations" ),
      pfimproject
    )
    derivNames = .derivativeNamesFromCompartments( pfimproject, 2L )

    list(
      duringInfusion = .combineInfusionPkPdEquations(
        pkEqODE$duringInfusion, pdEq, derivNames
      ),
      afterInfusion = .combineInfusionPkPdEquations(
        pkEqODE$afterInfusion, pdEq, derivNames
      )
    )
  }
