#' @title ModelAnalytic
#' @description Closed-form (non-ODE) PK/PD model with bolus or explicit dosing.
#' @param wrapperModelAnalytic                   Wrapper for the analytic solver.
#' @inheritParams Model
#' @param functionArgumentsModelAnalytic         A list with the function arguments of the wrapper.
#' @param functionArgumentsSymbolModelAnalytic   A list with the function argument symbols.
#' @param solverInputs                           A list with the solver inputs.
#' @inheritParams Model
#' @include Model.R
#' @include ModelODE.R
#' @export

ModelAnalytic = new_class(
  "ModelAnalytic",
  package = "PFIM",
  parent  = Model,
  properties = list(
    wrapperModelAnalytic                 = new_property(class_list, default = list()),
    functionArgumentsModelAnalytic       = new_property(class_list, default = list()),
    functionArgumentsSymbolModelAnalytic = new_property(class_list, default = list()),
    solverInputs                         = new_property(class_list, default = list())
  ))

convertPKModelAnalyticToPKModelODE = new_generic( "convertPKModelAnalyticToPKModelODE", c("pkModel") )

# Package-internal helper (not exported).
# Compiles an analytic equation set into a callable R function.
#
# equations        : named list/vector of character equation strings.
# functionArguments: character vector of all function arguments (signature).
# timeNames        : character vector of per-outcome time-variable names
#                    (e.g. "t_C1").  "\\bt\\b" occurrences in each equation
#                    are replaced by the corresponding timeName.
# returnNames      : character vector of names to include in the return value.

.buildAnalyticWrapper = function( equations, functionArguments, timeNames, returnNames ) {
  if ( length( equations ) == 0L ) return( function(...) NULL )

  # Map bare "t" to the per-equation time variable (t_<admin outcome>).
  eqNames = names( equations )
  body = map_chr( seq_along( eqNames ), function( i ) {
    nm = eqNames[[ i ]]
    line = sprintf( "%s = %s", nm, equations[[ nm ]] )
    t_var = if ( length( timeNames ) >= i && nzchar( timeNames[[ i ]] ) ) {
      timeNames[[ i ]]
    } else {
      paste0( "t_", nm )
    }
    line = str_replace_all( line, "\\bt\\b", t_var )
    if ( grepl( "\\bt\\b", line, perl = TRUE ) )
      stop( "Bare time variable 't' remains in analytic equation for '", nm, "'.", call. = FALSE )
    line
  })

  body = paste( body, collapse = "\n" )
  body = sprintf( "%s\nreturn(list(c(%s)))", body, paste( returnNames, collapse = ", " ) )
  eval( parse( text = sprintf( "function(%s) { %s }",
                               paste( functionArguments, collapse = ", " ), body ) ) )
}


method( defineModelWrapper, ModelAnalytic ) = function( model, evaluation ) {

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

  functionArgumentsWithAdmin    = unique( c( doseNames, parameterNames, timeNames ) )
  functionArgumentsWithNoAdmin  = unique( c( outcomesWithAdministration, parameterNames, timeNames ) )

  functionArgumentsSymbolWithAdmin   = map( functionArgumentsWithAdmin,   as.symbol )
  functionArgumentsSymbolWithNoAdmin = map( functionArgumentsWithNoAdmin, as.symbol )

  prop( model, "wrapperModelAnalytic" ) = list(
    functionDefinitionWithAdmin   = .buildAnalyticWrapper( equationsWithAdmin,   functionArgumentsWithAdmin,   timeNames, outputAdmin   ),
    functionDefinitionWithNoAdmin = .buildAnalyticWrapper( equationsWithNoAdmin, functionArgumentsWithNoAdmin, timeNames, outputNoAdmin )
  )
  prop( model, "functionArgumentsModelAnalytic" ) = list(
    functionArgumentsWithAdmin   = functionArgumentsWithAdmin,
    functionArgumentsWithNoAdmin = functionArgumentsWithNoAdmin
  )
  prop( model, "functionArgumentsSymbolModelAnalytic" ) = list(
    functionArgumentsSymbolWithAdmin   = functionArgumentsSymbolWithAdmin,
    functionArgumentsSymbolWithNoAdmin = functionArgumentsSymbolWithNoAdmin
  )
  prop( model, "outputNames" ) = unlist( names( equations ) )
  prop( model, "outcomesWithAdministration" ) = outcomesWithAdministration

  return( model )
}


method( defineModelAdministration, ModelAnalytic ) = function( model, arm ) {

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

    # Relative sampling times from each dose event: positive delay or 0.
    timeDose = timeDose |>
      map( ~ ifelse( samplings - .x > 0, samplings - .x, samplings ) ) |>
      reduce( cbind )

    # Number of unique doses active at each sampling time.
    indicesDoses = if ( is.null( dim( timeDose ) ) ) {
      1L
    } else {
      map_int( seq_len( nrow( timeDose ) ), ~ length( unique( timeDose[.x, ] ) ) )
    }
    list( data = data.frame( timeDose, indicesDoses ), dose = dose )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}

#' evaluateAnalyticCore: core analytic evaluation (shared by covariate path).
#' @name evaluateAnalyticCore
#' @param model A \code{ModelAnalytic} object.
#' @param arm   An \code{Arm} object.
#' @keywords internal

evaluateAnalyticCore = function( model, arm ) {

  parameters                 = prop( model, "modelParameters" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  outputNames                = prop( model, "outputNames" )
  samplings                  = prop( model, "samplings" )
  solverInputs               = prop( model, "solverInputs" )

  wrapperModelAnalytic          = prop( model, "wrapperModelAnalytic" )
  functionDefinitionWithAdmin   = wrapperModelAnalytic$functionDefinitionWithAdmin
  functionDefinitionWithNoAdmin = wrapperModelAnalytic$functionDefinitionWithNoAdmin

  # Extract mu values as a named numeric list (shared across all time steps)
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
                                c( paste0( "t_", outcome ),
                                   paste0( "dose_", outcome ) ) ) )

      evaluationOutcomeWithAdmin = sum( do.call( functionDefinitionWithAdmin, argsAdmin )[[ 1L ]] )
      currentArgsNoAdmin[[ outcome ]] = evaluationOutcomeWithAdmin

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


method( evaluateModel, ModelAnalytic ) = function( model, arm ) {
  if ( usesCovariateOccasionStructure( model ) )
    evaluateModelWithCovariates( model, arm, evaluateAnalyticCore )
  else
    evaluateAnalyticCore( model, arm )
}

#' Convert an analytic infusion model to ODE form
#' @name convertPKModelAnalyticToPKModelODE
#' @export

method( convertPKModelAnalyticToPKModelODE, ModelAnalytic ) = function( pkModel ) {

  pkModelEquations = prop( pkModel, "modelEquations" )
  eq = pluck( pkModelEquations, 1 )

  dtEquationPKsubstitute = D( parse( text = eq ), "t" ) |> deparse() |> str_c( collapse = "" )

  pkModelEquations = if ( str_detect( eq, "Cl" ) )
    str_c( dtEquationPKsubstitute, "+(Cl/V)*", eq, "- (Cl/V)*RespPK" )
  else
    str_c( dtEquationPKsubstitute, "+k*", eq, "- k*RespPK" )

  pkModelEquations |> str_replace_all( " ", "" ) |> (\(x) paste( Simplify(x) ))()
}


method( definePKModel, list( ModelAnalytic, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}


method( definePKPDModel, list( ModelAnalytic, ModelAnalytic, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {
    c( prop( pkModel, "modelEquations" ), prop( pdModel, "modelEquations" ) )
  }


method( definePKPDModel, list( ModelAnalytic, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    equations = c(
      convertPKModelAnalyticToPKModelODE( pkModel ),
      prop( pdModel, "modelEquations" )
    )
    eq = remapPkpdLibraryEquations( equations, pfimproject )
    set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ) ) )
  }
