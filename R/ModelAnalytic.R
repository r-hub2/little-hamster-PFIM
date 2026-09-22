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
#' @return An S7 object of class \code{ModelAnalytic}.
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

#' Convert analytic PK equations to ODE-compatible form
#' @param pkModel A PK model object.
#' @param ... Optional method arguments.
#' @name convertPKModelAnalyticToPKModelODE
#' @keywords internal
convertPKModelAnalyticToPKModelODE = new_generic( "convertPKModelAnalyticToPKModelODE", c("pkModel") )

# Compile analytic equations into an R function (equations, args, timeNames, returnNames).
#' Build executable wrapper for analytic equations.
#' @param equations Named character vector of analytic equations.
#' @param functionArguments Character vector of function argument names.
#' @param timeNames Character vector of time variable names.
#' @param returnNames Character vector of returned output names.
#' @return Function evaluating analytic equations.
#' @noRd
#' @keywords internal
.buildAnalyticWrapper = function( equations, functionArguments, timeNames, returnNames ) {
  if ( length( equations ) == 0L ) return( function(...) NULL )

  # Map bare "t" to the per-equation time variable (t_<this outcome>).
  eqNames = names( equations )
  body = map_chr( seq_along( eqNames ), function( i ) {
    nm = eqNames[[ i ]]
    line = sprintf( "%s = %s", nm, equations[[ nm ]] )
    t_var = if ( length( timeNames ) >= i &&
                 .pfimIsNonEmptyScalar( timeNames[[ i ]] ) ) {
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
  fn = eval( parse( text = sprintf( "function(%s) { %s }",
                                    paste( functionArguments, collapse = ", " ), body ) ) )
  # Free symbols (e.g. global covariates logtWT, SEX) resolve in knit/global env.
  environment( fn ) = .pfimUserSymbolEnv()
  fn
}


#' Compile analytic wrappers for administered and passive outcomes.
#'
#' Splits library equations into outcomes that receive dosing versus those that
#' do not (e.g. PD linked to PK), then builds two R functions via
#' \code{.buildAnalyticWrapper}. Administered wrappers take \code{dose_*} /
#' \code{t_*}; passive wrappers take the administered outcome values as inputs.
#' @return Updated model with compiled analytic wrappers.
#' @name defineModelWrapper
#' @keywords internal
method( defineModelWrapper, ModelAnalytic ) = function( model, evaluation ) {

  # Outcomes that receive dosing (from design administrations).
  outcomesWithAdministration = .getOutcomesFromEvaluation( evaluation )

  # Formal names injected into each analytic wrapper: dose_*, t_*, params.
  parameters     = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames      = paste0( "dose_", outcomesWithAdministration )
  timeNamesAdmin = paste0( "t_",    outcomesWithAdministration )

  # Split equations into administered vs passive (no-admin) outcomes.
  equations            = prop( evaluation, "modelEquations" )
  equationsWithAdmin   = equations[  names( equations ) %in% outcomesWithAdministration ]
  equationsWithNoAdmin = equations[ !names( equations ) %in% outcomesWithAdministration ]

  outputAdmin   = names( equationsWithAdmin )
  outputNoAdmin = names( equationsWithNoAdmin )
  timeNamesNoAdmin = .libraryEquationTimeNames( evaluation, outputNoAdmin )

  # Passive wrappers also receive administered outcome values (e.g. Conc for PD)
  # and their own t_<outcome>, not t_<administered>.
  functionArgumentsWithAdmin    = unique( c( doseNames, parameterNames, timeNamesAdmin ) )
  functionArgumentsWithNoAdmin  = unique( c( outcomesWithAdministration, parameterNames, timeNamesNoAdmin ) )

  functionArgumentsSymbolWithAdmin   = map( functionArgumentsWithAdmin,   as.symbol )
  functionArgumentsSymbolWithNoAdmin = map( functionArgumentsWithNoAdmin, as.symbol )

  # Two evaluators: administered outcomes and passive/linked outcomes.
  prop( model, "wrapperModelAnalytic" ) = list(
    functionDefinitionWithAdmin   = .buildAnalyticWrapper( equationsWithAdmin,   functionArgumentsWithAdmin,   timeNamesAdmin,   outputAdmin   ),
    functionDefinitionWithNoAdmin = .buildAnalyticWrapper( equationsWithNoAdmin, functionArgumentsWithNoAdmin, timeNamesNoAdmin, outputNoAdmin )
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


#' Build per-outcome dose and relative-time tables for the analytic solver.
#'
#' For \code{tau != 0}, expands a single dose into a regular grid
#' \code{0, tau, 2*tau, ...} up to the last sampling time. At each sampling,
#' records relative times since each past dose and how many distinct dose
#' contributions are active (for linear superposition).
#' @return Updated model with \code{samplings} and \code{solverInputs}.
#' @name defineModelAdministration
#' @keywords internal
method( defineModelAdministration, ModelAnalytic ) = function( model, arm ) {

  administrations            = prop( arm,   "administrations" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  samplings                  = .analyticSamplingGrid( arm )

  solverInputs = map( administrations, function( adm ) {
    dosing = .alignAdministrationDosing( adm )
    tbl    = .analyticRelativeDoseTable(
      samplings, dosing$timeDose, dosing$dose, prop( adm, "tau" )
    )
    list( data = tbl$data, dose = tbl$dose )
  }) |> set_names( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}

#' Evaluate analytic bolus concentrations at all sampling times.
#'
#' For each observation time and administered outcome, evaluates the closed-form
#' formula over all active doses (superposition via \code{sum}), then evaluates
#' passive (non-admin) equations with the administered prediction injected into
#' the shared argument list.
#' @param model A \code{ModelAnalytic} object.
#' @param arm   An \code{Arm} object.
#' @return Named list of output data frames at requested sampling times.
#' @name evaluateAnalyticCore
#' @keywords internal
evaluateAnalyticCore = function( model, arm ) {

  parameters                 = prop( model, "modelParameters" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  outputNames                = prop( model, "outputNames" )
  samplings                  = prop( model, "samplings" )
  solverInputs               = prop( model, "solverInputs" )

  # Compiled wrappers from defineModelWrapper().
  wrappers = .analyticBindUserEnv(
    prop( model, "wrapperModelAnalytic" )$functionDefinitionWithAdmin,
    prop( model, "wrapperModelAnalytic" )$functionDefinitionWithNoAdmin
  )
  fnAdmin   = wrappers[[ 1L ]]
  fnNoAdmin = wrappers[[ 2L ]]
  mu        = .extractMu( parameters )

  tmp = .analyticEvalGrid(
    samplings, outcomesWithAdministration, as.list( mu ),
    function( iterTime, outcome, args ) {
      data         = solverInputs[[ outcome ]]$data
      dose         = solverInputs[[ outcome ]]$dose
      indicesDoses = data$indicesDoses[ iterTime ]
      argsAdmin    = c(
        args,
        set_names(
          list(
            as.numeric( data[ iterTime, seq_len( indicesDoses ) ] ),
            as.numeric( dose[ seq_len( indicesDoses ) ] )
          ),
          c( paste0( "t_", outcome ), paste0( "dose_", outcome ) )
        )
      )
      admin = sum( do.call( fnAdmin, argsAdmin )[[ 1L ]] )
      args[[ outcome ]] = admin
      argsNoAdmin = .pfimFillMissingTimeFormals( fnNoAdmin, args, samplings[[ iterTime ]] )
      list( admin = admin, noAdmin = do.call( fnNoAdmin, argsNoAdmin )[[ 1L ]], args = args )
    }
  )
  .analyticFinishEvaluation( tmp, outputNames, arm )
}


#' Dispatch: covariate/occasion structure -> specialised path; else core evaluator.
#' @return Named list of output data frames at requested sampling times.
#' @name evaluateModel
#' @keywords internal
method( evaluateModel, ModelAnalytic ) = function( model, arm ) {
  .analyticDispatchEvaluate( model, arm, evaluateAnalyticCore )
}

#' Convert one analytic PK concentration formula to an ODE derivative.
#'
#' Differentiates the closed form w.r.t. time, then adds elimination in
#' clearance (\code{Cl/V}) or rate-constant (\code{k}) form so the RHS matches
#' the library ODE convention.
#' @param equation Character analytic PK formula (library notation).
#' @param stateName State variable replaced in the elimination term (default \code{RespPK}).
#' @return Character scalar ODE right-hand side.
#' @noRd
#' @keywords internal
.convertAnalyticPkExprToOde = function( equation, stateName = "RespPK" ) {
  # d(analytic)/dt, then cancel elimination of the closed form and re-introduce
  # elimination of the ODE state (Cl/V or k depending on parameterisation).
  dt = D( parse( text = equation ), "t" ) |> deparse() |> str_c( collapse = "" )
  out = if ( str_detect( equation, "Cl" ) ) {
    str_c( dt, "+(Cl/V)*", equation, "-(Cl/V)*", stateName )
  } else {
    str_c( dt, "+k*", equation, "-k*", stateName )
  }
  out |> str_replace_all( " ", "" ) |> (\(x) paste( Simplify( x ) ))()
}

#' Convert the first analytic PK equation to ODE-compatible form.
#' @param pkModel First argument of generic.
#' @return Character vector of ODE-compatible PK equations.
#' @name convertPKModelAnalyticToPKModelODE
#' @keywords internal
method( convertPKModelAnalyticToPKModelODE, ModelAnalytic ) = function( pkModel ) {
  .convertAnalyticPkExprToOde( pluck( prop( pkModel, "modelEquations" ), 1 ) )
}


#' Return PK equations stored on the analytic model (library / user).
#' @param pkModel First argument of generic.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return List of PK equations from \code{pkModel}.
#' @name definePKModel
#' @keywords internal
method( definePKModel, list( ModelAnalytic, PFIMProject ) ) = function( pkModel, pfimproject ) {
  prop( pkModel, "modelEquations" )
}


#' Concatenate analytic PK and analytic PD equation lists.
#' @param pkModel First argument of generic.
#' @param pdModel Second model combined with PK equations.
#' @param pfimproject \code{PFIMProject} object (unused).
#' @return Concatenated analytic PK and PD equation list.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelAnalytic, ModelAnalytic, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {
    c( prop( pkModel, "modelEquations" ), prop( pdModel, "modelEquations" ) )
  }


#' Convert analytic PK to ODE, append ODE PD, and remap library compartments.
#' @param pkModel First argument of generic.
#' @param pdModel ODE PD model whose equations are appended.
#' @param pfimproject \code{PFIMProject} used for equation remapping.
#' @return Named list of remapped PK/PD ODE equations.
#' @name definePKPDModel
#' @keywords internal
method( definePKPDModel, list( ModelAnalytic, ModelODE, PFIMProject ) ) =
  function( pkModel, pdModel, pfimproject ) {

    # Analytic PK -> ODE derivative, then join with PD and remap RespPK/E tokens.
    equations = c(
      convertPKModelAnalyticToPKModelODE( pkModel ),
      prop( pdModel, "modelEquations" )
    )
    eq = remapPkpdLibraryEquations( equations, pfimproject )
    set_names( eq, .derivativeNamesFromCompartments( pfimproject, length( eq ), names( eq ) ) )
  }
