#' @title ModelODE
#' @description Base class for models integrated with \code{deSolve}.
#' @inheritParams Model
#' @include Model.R
#' @return An S7 object of class \code{ModelODE}.
#' @export

ModelODE = new_class( "ModelODE", package = "PFIM", parent = Model )

#' Check whether a value is a finite scalar.
#' @param x Object to validate.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.scalar_finite = function( x ) {
  length( x ) == 1L && !is.na( x ) && is.finite( x )
}

#' Extract representative mu values from model parameters.
#'
#' Preference order: distribution mu when omega is present; else mu when value
#' is missing; else the reserved \code{value} slot; else mu again. Used to seed
#' ODE/analytic evaluation environments.
#' @param parameters List of \code{ModelParameter} objects.
#' @return Named numeric vector of parameter values used by ODE solvers.
#' @noRd
#' @keywords internal
.extractMu = function( parameters ) {
  set_names(
    map_dbl( parameters, function( parameter ) {
      distribution = prop( parameter, "distribution" )
      value        = prop( parameter, "value" )

      if ( !is.null( distribution ) ) {
        omega = prop( distribution, "omega" )
        mu    = prop( distribution, "mu" )
        # Random-effect parameter: evaluate at typical value mu.
        if ( .scalar_finite( omega ) && omega != 0 )
          return( mu )
        if ( .scalar_finite( mu ) && mu != 0 && !.scalar_finite( value ) )
          return( mu )
      }
      if ( .scalar_finite( value ) )
        return( value )
      if ( !is.null( distribution ) && .scalar_finite( prop( distribution, "mu" ) ) )
        return( prop( distribution, "mu" ) )
      NA_real_
    }),
    map_chr( parameters, ~ prop( .x, "name" ) )
  )
}

#' Build sorted unique sampling grid including time zero.
#'
#' Always prepends \code{0} so deSolve starts at the initial condition even when
#' the design's first sample is strictly positive.
#' @param samplingTimes List of \code{SamplingTimes} objects.
#' @return Numeric vector of sampling times.
#' @noRd
#' @keywords internal
.buildSamplings = function( samplingTimes ) {
  s = map( samplingTimes, ~ prop( .x, "samplings" ) ) |> unlist() |> sort() |> unique()
  c( 0.0, s[ s != 0 ] )
}

#' Extract administered outcomes from project design arms.
#' @param evaluation \code{PFIMProject} or \code{Evaluation} object.
#' @return Character vector of outcome names with administrations.
#' @noRd
#' @keywords internal
.getOutcomesFromEvaluation = function( evaluation ) {
  prop( evaluation, "designs" ) |>
    map( \(d) prop( d, "arms" ) ) |>
    list_flatten() |>
    map( \(arm) prop( arm, "administrations" ) ) |>
    list_flatten() |>
    map_chr( \(adm) prop( adm, "outcome" ) ) |>
    unique()
}

#' Flatten equation list names (handles during/after infusion nests).
#' @noRd
#' @keywords internal
.pfimModelEquationNames = function( equations ) {
  if ( !is.list( equations ) || !length( equations ) )
    return( character( 0 ) )
  nms = names( equations )
  nested = vapply( equations, is.list, logical( 1L ) )
  if ( !is.null( nms ) && any( is.na( nms ) ) )
    .pfimStop( "model equation names must be non-empty (got NA)." )
  out = if ( is.null( nms ) ) character( 0 ) else
    nms[ !nested & !is.na( nms ) & nzchar( nms ) ]
  if ( any( nested ) )
    out = c( out, unlist( lapply( equations[ nested ], .pfimModelEquationNames ),
                          use.names = FALSE ) )
  unique( out )
}

#' Compartment / catalogue names from model equations (\code{Deriv_} stripped).
#' @noRd
#' @keywords internal
.pfimOdeCompartmentNames = function( equations ) {
  sub( "^Deriv_", "", .pfimModelEquationNames( equations ) )
}

#' State names in deSolve / wrapper order (\code{variableNames} or equations).
#' @noRd
#' @keywords internal
.pfimOdeStateNames = function( model ) {
  vn = prop( model, "variableNames" )
  if ( length( vn ) ) {
    vn = as.character( vn )
    if ( any( is.na( vn ) | !nzchar( vn ) ) )
      .pfimStop( "ODE state names must be non-empty (got NA or \"\")." )
    return( vn )
  }
  .pfimOdeCompartmentNames( prop( model, "modelEquations" ) )
}

#' Map output aliases (e.g. \code{RespPK}) to state names from model formulas.
#' @noRd
#' @keywords internal
.pfimOutputStateAliases = function( model ) {
  of = .getModelOutputFormulas( model )
  if ( !length( of ) ) return( character( 0 ) )
  vals = vapply( names( of ), function( nm ) {
    x = of[[ nm ]]
    if ( is.character( x ) && length( x ) == 1L ) return( x )
    if ( is.symbol( x ) || is.name( x ) ) return( as.character( x ) )
    NA_character_
  }, character( 1L ) )
  stats::setNames( vals, names( of ) )
}

#' Resolve an administration outcome to an ODE state (alias via \code{outputs}).
#' @noRd
#' @keywords internal
.pfimResolveOutcomeToState = function( outcome, states, aliases ) {
  if ( outcome %in% states ) return( outcome )
  if ( length( aliases ) && outcome %in% names( aliases ) ) {
    target = unname( aliases[[ outcome ]] )
    if ( .pfimIsNonEmptyScalar( target ) && target %in% states )
      return( target )
  }
  outcome
}

#' Reorder a named IC vector to match ODE state order (missing states filled with 0).
#' @noRd
#' @keywords internal
.pfimAlignOdeStates = function( ic, states ) {
  if ( !length( states ) ) return( unlist( ic ) )
  ic = unlist( ic )
  if ( length( ic ) && !is.null( names( ic ) ) ) {
    unknown = setdiff( names( ic ), states )
    unknown = unknown[ !is.na( unknown ) & nzchar( unknown ) ]
    if ( length( unknown ) )
      .pfimStop(
        "initialConditions name(s) not in the ODE states: ",
        paste( unknown, collapse = ", " ),
        ". Declared states: ", paste( states, collapse = ", " ), "."
      )
  }
  out = stats::setNames( rep( 0, length( states ) ), states )
  if ( length( ic ) && !is.null( names( ic ) ) ) {
    hit = intersect( names( ic ), states )
    out[ hit ] = ic[ hit ]
  }
  out
}

#' Suffixes of \code{dose_*} tokens in equation text.
#' @noRd
#' @keywords internal
.pfimDoseSuffixes = function( equations ) {
  .pfimPrefixedTokenSuffixes( equations, "dose_" )
}

#' Suffixes of \code{Tinf_*} tokens in equation text.
#' @noRd
#' @keywords internal
.pfimTinfSuffixes = function( equations ) {
  .pfimPrefixedTokenSuffixes( equations, "Tinf_" )
}

#' Suffixes of \code{prefix*} tokens in equation text.
#' @noRd
#' @keywords internal
.pfimPrefixedTokenSuffixes = function( equations, prefix ) {
  txt = unlist( equations, use.names = FALSE )
  if ( !length( txt ) ) return( character( 0 ) )
  pat = paste0( prefix, "([A-Za-z][A-Za-z0-9]*)" )
  hits = unlist( regmatches( txt, gregexpr( pat, txt, perl = TRUE ) ), use.names = FALSE )
  unique( sub( paste0( "^", prefix ), "", hits ) )
}

#' Allowed administration and sampling outcome names for a project.
#' @noRd
#' @keywords internal
.pfimAllowedOutcomes = function( project ) {
  equations = projectProp( project, "modelEquations" )
  compartments = .pfimOdeCompartmentNames( equations )
  dose_names = .pfimDoseSuffixes( equations )
  outputs = projectProp( project, "outputs" )
  out_names = names( outputs )
  out_vals = unname( unlist( outputs, use.names = FALSE ) )
  alias_admin = character( 0 )
  if ( length( outputs ) ) {
    mapped = vapply( out_names, function( nm ) {
      val = unlist( outputs[[ nm ]], use.names = FALSE )
      if ( .pfimIsNonEmptyScalar( val ) &&
           val %in% c( compartments, dose_names ) )
        return( nm )
      NA_character_
    }, character( 1L ) )
    alias_admin = mapped[ !is.na( mapped ) ]
  }
  list(
    admin = unique( c( compartments, dose_names, alias_admin ) ),
    sampling = unique( c( compartments, dose_names, out_names, out_vals ) )
  )
}

#' Stop when arm administration / sampling outcomes are not model states or aliases.
#' @noRd
#' @keywords internal
.pfimValidateProjectOutcomes = function( project ) {
  equations = projectProp( project, "modelEquations" )
  if ( !length( equations ) ) return( invisible( NULL ) )
  allowed = .pfimAllowedOutcomes( project )
  if ( !length( allowed$admin ) ) return( invisible( NULL ) )
  designs = projectProp( project, "designs" )
  for ( design in designs ) {
    for ( arm in prop( design, "arms" ) ) {
      for ( adm in prop( arm, "administrations" ) ) {
        o = prop( adm, "outcome" )
        if ( .pfimIsNonEmptyScalar( o ) && !o %in% allowed$admin )
          .pfimStop(
            "Administration(outcome = \"", o,
            "\") is not a model state, dose_* target, or outputs alias ",
            "mapping to a state. Allowed: ",
            paste( allowed$admin, collapse = ", " ), "."
          )
      }
      for ( st in prop( arm, "samplingTimes" ) ) {
        o = prop( st, "outcome" )
        if ( .pfimIsNonEmptyScalar( o ) && !o %in% allowed$sampling )
          .pfimStop(
            "SamplingTimes(outcome = \"", o,
            "\") is not a model state or output alias. Allowed: ",
            paste( allowed$sampling, collapse = ", " ), "."
          )
      }
    }
  }
  .pfimValidateDoseTokenAdministrations( project, equations )
  invisible( NULL )
}

#' Stop when \code{dose_*} / \code{Tinf_*} tokens have no matching administration.
#' @noRd
#' @keywords internal
.pfimValidateDoseTokenAdministrations = function( project, equations ) {
  tokens = unique( c( .pfimDoseSuffixes( equations ), .pfimTinfSuffixes( equations ) ) )
  if ( !length( tokens ) ) return( invisible( NULL ) )
  outputs = projectProp( project, "outputs" )
  aliases = if ( length( outputs ) ) {
    vapply( names( outputs ), function( nm ) {
      val = unlist( outputs[[ nm ]], use.names = FALSE )
      if ( .pfimIsNonEmptyScalar( val ) ) val else NA_character_
    }, character( 1L ) )
  } else {
    character( 0 )
  }
  aliases = aliases[ !is.na( aliases ) ]
  resolved = character( 0 )
  designs = projectProp( project, "designs" )
  for ( design in designs ) {
    for ( arm in prop( design, "arms" ) ) {
      for ( adm in prop( arm, "administrations" ) ) {
        o = prop( adm, "outcome" )
        if ( !.pfimIsNonEmptyScalar( o ) ) next
        resolved = c( resolved, o )
        if ( length( aliases ) && o %in% names( aliases ) )
          resolved = c( resolved, unname( aliases[[ o ]] ) )
        if ( length( aliases ) && o %in% aliases )
          resolved = c( resolved, names( aliases )[ aliases == o ] )
      }
    }
  }
  resolved = unique( resolved )
  missing = setdiff( tokens, resolved )
  if ( length( missing ) )
    .pfimStop(
      "equations reference dose_/Tinf_ \"", paste( missing, collapse = ", " ),
      "\" but no Administration(outcome) resolves to ",
      if ( length( missing ) == 1L ) "that name" else "those names",
      ". Administer the matching state (or an outputs alias). Allowed administrations: ",
      paste( unique( c( resolved, tokens ) ), collapse = ", " ), "."
    )
  invisible( NULL )
}

#' Operators and math functions allowed in initial-condition expressions.
#' @noRd
#' @keywords internal
.pfimIcSafeNames = c(
  "+", "-", "*", "/", "^", "(", "{", "c",
  "exp", "log", "log10", "log1p", "expm1", "sqrt", "abs", "pmax", "pmin"
)

#' Build an evaluation environment for initial-condition expressions.
#'
#' Uses \code{baseenv()} as parent so arithmetic operators resolve without
#' capturing the caller environment.
#' @param mu Named numeric vector of typical parameter values.
#' @param extra Optional named list of extra bindings (e.g. \code{dose_*}).
#' @return Environment ready for \code{eval()}.
#' @noRd
#' @keywords internal
.odeIcEvalEnv = function( mu, extra = list() ) {
  env = new.env( parent = baseenv() )
  list2env( c( as.list( mu ), extra ), envir = env )
  env
}

#' Evaluate one IC expression in a sealed env after checking free names.
#' @noRd
#' @keywords internal
.pfimEvalIcExpr = function( ic, env, label ) {
  parsed = parse( text = ic )
  if ( !length( parsed ) )
    .pfimStop( "initial condition '", label, "' is empty." )
  used = unique( all.names( parsed[[ 1L ]] ) )
  allowed = c( .pfimIcSafeNames, ls( envir = env, all.names = TRUE ) )
  bad = setdiff( used, allowed )
  if ( length( bad ) )
    .pfimStop(
      "initial condition '", label, "' uses unknown name(s): ",
      paste( bad, collapse = ", " ), "."
    )
  eval( parsed, envir = env )
}

#' Resolve library equation keys that correspond to administered outcomes.
#'
#' Design administrations name user outcomes (\code{Cc}); library equations use
#' keys like \code{RespPK}. Map through \code{outputs} so during/after splits
#' match the catalogue, not the design labels.
#' @param evaluation \code{PFIMProject} or \code{Evaluation} object.
#' @return Character vector of output names to treat as administered.
#' @noRd
#' @keywords internal
.administeredLibraryOutcomeNames = function( evaluation ) {
  adminOutcomes = .getOutcomesFromEvaluation( evaluation )
  outputs       = prop( evaluation, "outputs" )
  if ( length( outputs ) == 0L ) return( adminOutcomes )
  outputValues = unlist( outputs, use.names = TRUE )
  outputNames  = names( outputValues )
  if ( is.null( outputNames ) ) outputNames = as.character( outputValues )
  kept = outputValues %in% adminOutcomes
  if ( any( kept ) ) outputNames[ kept ] else adminOutcomes
}

#' Build time-variable name for one library outcome.
#' @param evaluation \code{PFIMProject} carrying output mapping.
#' @param libraryOutcomeName Character outcome name from library equations.
#' @return Character scalar time-variable name.
#' @noRd
#' @keywords internal
.libraryEquationTimeName = function( evaluation, libraryOutcomeName ) {
  outputs = prop( evaluation, "outputs" )
  admin   = outputs[[ libraryOutcomeName ]]
  paste0( "t_", if ( .pfimIsNonEmptyScalar( admin ) ) admin else libraryOutcomeName )
}

#' Build time-variable names for multiple library outcomes.
#' @param evaluation \code{PFIMProject} carrying output mapping.
#' @param libraryOutcomeNames Character vector of library outcome names.
#' @return Character vector of time-variable names.
#' @noRd
#' @keywords internal
.libraryEquationTimeNames = function( evaluation, libraryOutcomeNames ) {
  libraryOutcomeNames = as.character( libraryOutcomeNames )
  if ( !length( libraryOutcomeNames ) ) return( character( 0 ) )
  map_chr( libraryOutcomeNames, ~ .libraryEquationTimeName( evaluation, .x ) )
}

#' Build an executable ODE wrapper from equation strings.
#' @param equations Named character vector of equation definitions.
#' @param functionArguments Character vector of wrapper argument names.
#' @return Function returning derivative vectors for ODE integration.
#' @noRd
#' @keywords internal
.buildODEWrapper = function( equations, functionArguments ) {
  body = paste(
    paste( map_chr( names( equations ),
                    ~ sprintf( "%s = %s", .x, equations[[.x]] ) ),
           collapse = "\n" ),
    sprintf( "return(list(c(%s)))", paste( names( equations ), collapse = ", " ) ),
    sep = "\n"
  )
  fn = eval( parse( text = sprintf( "function(%s) { %s }",
                                    paste( functionArguments, collapse = ", " ), body ) ) )
  # Free symbols (global covariates) resolve in knit/global env, not the package ns.
  environment( fn ) = .pfimUserSymbolEnv()
  fn
}

# Memoize output-formula maps / parsed expressions across repeated ODE evals.
.pfimOutputFormulaCache = new.env( parent = emptyenv() )
.pfimParsedOutputFormulaCache = new.env( parent = emptyenv() )

#' Parse one output-formula entry when needed.
#' @param x Character formula expression or already parsed object.
#' @return Parsed expression or original object.
#' @noRd
#' @keywords internal
.parseOutputFormulaEntry = function( x ) {
  if ( is.character( x ) ) parse( text = x ) else x
}

#' Resolve output column name in ODE simulation data.
#'
#' deSolve may emit \code{RespPK.Cc}-style dotted names when an output formula
#' aliases a state; try the alias before falling back to \code{out_name}.
#' @param out_name Character requested output name.
#' @param evaluationModelTmp Data frame returned by ODE simulation.
#' @param outputFormula Optional output formula mapping.
#' @return Character column name present in \code{evaluationModelTmp}.
#' @noRd
#' @keywords internal
.odeColumnForOutput = function( out_name, evaluationModelTmp, outputFormula ) {
  cols = names( evaluationModelTmp )
  if ( out_name %in% cols ) return( out_name )
  if ( length( outputFormula ) && out_name %in% names( outputFormula ) ) {
    target = outputFormula[[ out_name ]]
    if ( is.character( target ) && length( target ) == 1L ) {
      dotted = paste0( out_name, ".", target )
      if ( dotted %in% cols ) return( dotted )
      if ( target %in% cols ) return( target )
    }
  }
  out_name
}

#' Build cache lookup key for model output formulas.
#' @param model \code{Model} object containing output names.
#' @return Character scalar lookup key.
#' @noRd
#' @keywords internal
.pfimOutputFormulaLookupKey = function( model ) {
  outputs = prop( model, "outputFormula" ) %||% list()
  formula_sig = if ( length( outputs ) )
    paste( names( outputs ), as.character( outputs ), sep = "=", collapse = "|" )
  else ""
  paste(
    class( model )[[ 1L ]],
    paste( prop( model, "outputNames" ), collapse = "," ),
    formula_sig,
    sep = "::"
  )
}

#' Store output formula mapping in cache.
#'
#' Clears the parsed-expression memo for this key so string updates are not
#' served as stale \code{parse()} results.
#' @param model \code{Model} object serving as cache key source.
#' @param outputs Named list of output formulas.
#' @return Invisibly returns updated \code{model}.
#' @noRd
#' @keywords internal
.setModelOutputFormulas = function( model, outputs ) {
  lookupKey = .pfimOutputFormulaLookupKey( model )
  .pfimCacheStore( .pfimOutputFormulaCache, lookupKey, outputs )
  .pfimParsedOutputFormulaCache[[ lookupKey ]] = NULL
  prop( model, "outputFormula" ) = outputs
  invisible( model )
}

#' Retrieve output formulas from cache or model slot.
#' @param model \code{Model} object serving as cache key source.
#' @return Named list of output formulas.
#' @noRd
#' @keywords internal
.getModelOutputFormulas = function( model ) {
  cached = .pfimOutputFormulaCache[[ .pfimOutputFormulaLookupKey( model ) ]]
  if ( !is.null( cached ) ) return( cached )
  prop( model, "outputFormula" )
}

#' Retrieve parsed output formulas with memoization.
#' @param model \code{Model} object serving as cache key source.
#' @return Named list of parsed output expressions.
#' @noRd
#' @keywords internal
.getOutputFormulaParsed = function( model ) {
  key = .pfimOutputFormulaLookupKey( model )
  cached = .pfimParsedOutputFormulaCache[[ key ]]
  if ( !is.null( cached ) ) return( cached )
  outputs = .getModelOutputFormulas( model )
  if ( !length( outputs ) ) return( outputs )
  parsed = map( outputs, .parseOutputFormulaEntry )
  .pfimCacheStore( .pfimParsedOutputFormulaCache, key, parsed )
  parsed
}

#' Evaluate initial-condition expressions against parameter environment.
#'
#' Numeric ICs pass through; character ICs are \code{eval}'d with typical
#' values bound locally so expressions like \code{"dose/V"} resolve.
#' @param model \code{ModelODE} object providing parameter values.
#' @param arm \code{Arm} object containing initial-condition definitions.
#' @return Named numeric vector of initial conditions.
#' @noRd
#' @keywords internal
.evalInitialConditionsImpl = function( model, arm ) {
  mu = .extractMu( prop( model, "modelParameters" ) )
  env = .odeIcEvalEnv( mu )
  imap( prop( arm, "initialConditions" ), function( ic, compartment ) {
    if ( is.numeric( ic ) )
      return( ic )
    .pfimEvalIcExpr( ic, env, compartment )
  } ) |> unlist()
}

#' Evaluate initial conditions for ODE models (numeric or expression IC).
#' @return Named numeric vector of evaluated initial conditions.
#' @name evaluateInitialConditions
#' @keywords internal
method( evaluateInitialConditions, ModelODE ) = function( model, arm ) {
  .evalInitialConditionsImpl( model, arm )
}
