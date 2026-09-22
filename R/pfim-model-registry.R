# Model-class registry: detect() tried in registration order.
#
# `pfim_register_model_class()` stores factory + optional detect() handlers.
# Resolution: explicit `modelClass` on the project, else first matching detect(),
# else legacy fallback in model-type-dispatch.R.
#
# Built-in catch-alls (ModelAnalytic / ModelODEDoseNotInEquations) have no
# detect() so custom classes registered after package load can still match.
# Specific builtins stay before custom entries when registered at load time;
# call pfim_register_model_class() early (or before a catch-all) for overrides.

#' @include model-type-dispatch.R
#' @include ModelAnalytic.R
#' @include ModelAnalyticSteadyState.R
#' @include ModelAnalyticInfusion.R
#' @include ModelAnalyticInfusionSteadyState.R
#' @include ModelODE.R
#' @include model-ode-bolus.R
#' @include ModelODEInfusion.R
#' @include ModelODEInfusionDoseInEquation.R
#' @keywords internal
NULL

.pfimModelClassRegistry = new.env( parent = emptyenv() )

#' Ordered list of registered model class names (registration order).
#' @noRd
#' @keywords internal
.pfimModelClassOrder = function() {
  if ( exists( "__order__", envir = .pfimModelClassRegistry, inherits = FALSE ) )
    get( "__order__", envir = .pfimModelClassRegistry )
  else
    character( 0L )
}

#' Append a class name to the registry order (no duplicates).
#' @noRd
#' @keywords internal
.pfimModelClassOrderAppend = function( className ) {
  ord = .pfimModelClassOrder()
  if ( !className %in% ord )
    assign( "__order__", c( ord, className ), envir = .pfimModelClassRegistry )
}

#' Invoke a registered detect() handler (pfimproject or legacy eq/ic).
#'
#' Prefers the project object whenever the signature looks project-oriented
#' (\code{function(project)}, \code{function(project, ...)}, \code{function(...)}).
#' Only plain two-argument \code{function(equations, ic)} stays on the legacy path -
#' so \code{function(project, ...)} is never fed equation lists.
#' @noRd
#' @keywords internal
.pfimCallModelDetect = function( detectFn, pfimproject ) {
  fmls = formals( detectFn )
  nms  = names( fmls )
  if ( is.null( nms ) ) nms = rep( "", length( fmls ) )
  n = length( fmls )

  projectish = c( "pfimproject", "project", "evaluation", "optimization", "x", "object" )
  legacy_eq  = c( "equations", "eq", "modelEquations" )
  legacy_ic  = c( "initialConditions", "ic", "init" )

  # Explicit legacy: two+ named args that are not a project (and no ...).
  if ( n >= 2L && !( "..." %in% nms ) &&
       !nms[[ 1L ]] %in% projectish &&
       ( nms[[ 1L ]] %in% legacy_eq || nms[[ 1L ]] == "" ) ) {
    eqs = projectProp( pfimproject, "modelEquations" )
    ic  = .initialConditionsFromProject( pfimproject )
    return( isTRUE( detectFn( eqs, ic ) ) )
  }

  # Default / documented API: pass the project.
  isTRUE( detectFn( pfimproject ) )
}

#' Register a model class
#' @param className Character S7 class name (e.g. \code{"ModelAnalytic"}).
#' @param factory Zero-argument function returning a \code{Model} object.
#' @param detect Optional \code{function(pfimproject)} returning logical.
#'   Legacy \code{function(equations, initialConditions)} is still accepted.
#' @return Invisibly, \code{NULL}.
#' @export
pfim_register_model_class = function( className, factory, detect = NULL ) {
  if ( !is.character( className ) || !.pfimIsNonEmptyScalar( className ) )
    stop( "className must be a non-empty character string.", call. = FALSE )
  if ( !is.function( factory ) )
    stop( "factory must be a function.", call. = FALSE )
  if ( !is.null( detect ) && !is.function( detect ) )
    stop( "detect must be NULL or a function.", call. = FALSE )
  assign(
    className,
    list( factory = factory, detect = detect ),
    envir = .pfimModelClassRegistry
  )
  .pfimModelClassOrderAppend( className )
  invisible( NULL )
}

#' Resolve model class name for a project
#' @param pfimproject A \code{PFIMProject}, \code{Evaluation}, or \code{Optimization}.
#' @return Character class name with attribute \code{source} (\code{"registry"} or \code{"legacy"}).
#' @export
pfim_resolve_model_class = function( pfimproject ) {
  info = .pfimResolveModelClassInfo( pfimproject )
  structure( info$class, source = info$source )
}

#' Resolve class name and provenance (explicit / registry / legacy).
#' @noRd
#' @keywords internal
.pfimResolveModelClassInfo = function( pfimproject ) {
  explicit = projectProp( pfimproject, "modelClass" )
  if ( .pfimIsNonEmptyScalar( explicit ) ) {
    if ( !exists( explicit, envir = .pfimModelClassRegistry, inherits = FALSE ) )
      stop( "Unknown model class: ", explicit, call. = FALSE )
    return( list( class = explicit, source = "explicit" ) )
  }

  # First matching detect() wins (specific builtins first; catch-alls have none).
  order = .pfimModelClassOrder()
  idx = detect_index( order, function( nm ) {
    entry = get( nm, envir = .pfimModelClassRegistry )
    !is.null( entry$detect ) && .pfimCallModelDetect( entry$detect, pfimproject )
  } )
  if ( idx > 0L )
    return( list( class = order[[ idx ]], source = "registry" ) )

  list(
    class  = .selectModelClass( .detectModelFeaturesFromProject( pfimproject ) ),
    source = "legacy"
  )
}

#' Instantiate a registered model class via its factory.
#' @noRd
#' @keywords internal
.pfimInstantiateModelClass = function( className ) {
  if ( !exists( className, envir = .pfimModelClassRegistry, inherits = FALSE ) )
    stop( "Unknown model class: ", className, call. = FALSE )
  get( className, envir = .pfimModelClassRegistry )$factory()
}

#' Persist resolved class name onto the project when modelClass was empty.
#' @noRd
#' @keywords internal
.pfimSyncModelClass = function( pfimproject, className ) {
  current = projectProp( pfimproject, "modelClass" )
  if ( !.pfimIsNonEmptyScalar( current ) )
    projectProp( pfimproject, "modelClass" ) = className
  invisible( className )
}

#' Register built-in analytic and ODE model classes (most specific first).
#'
#' Catch-all classes are registered without \code{detect} so legacy selection
#' (and user-registered detectors) can still run.
#' @noRd
#' @keywords internal
.pfimInitBuiltinModelRegistry = function() {
  feats = function( pfimproject ) .detectModelFeaturesFromProject( pfimproject )

  pfim_register_model_class(
    "ModelODEBolus",
    function() ModelODEBolus(),
    detect = function( pfimproject ) {
      f = feats( pfimproject ); f$isODE && f$doseInInitialConditions
    }
  )
  pfim_register_model_class(
    "ModelODEInfusionDoseInEquation",
    function() ModelODEInfusionDoseInEquation(),
    detect = function( pfimproject ) {
      f = feats( pfimproject ); f$isODE && f$hasInfusion && f$doseInEquation
    }
  )
  pfim_register_model_class(
    "ModelODEDoseInEquations",
    function() ModelODEDoseInEquations(),
    detect = function( pfimproject ) {
      f = feats( pfimproject ); f$isODE && f$doseInEquation
    }
  )
  # Catch-all ODE: no detect - legacy .selectModelClass picks it when needed.
  pfim_register_model_class(
    "ModelODEDoseNotInEquations",
    function() ModelODEDoseNotInEquations(),
    detect = NULL
  )
  pfim_register_model_class(
    "ModelAnalyticInfusionSteadyState",
    function() ModelAnalyticInfusionSteadyState(),
    detect = function( pfimproject ) {
      f = feats( pfimproject )
      !f$isODE && f$hasInfusion && f$doseInEquation && f$hasTau
    }
  )
  pfim_register_model_class(
    "ModelAnalyticInfusion",
    function() ModelAnalyticInfusion(),
    detect = function( pfimproject ) {
      f = feats( pfimproject )
      !f$isODE && f$hasInfusion && f$doseInEquation
    }
  )
  pfim_register_model_class(
    "ModelAnalyticSteadyState",
    function() ModelAnalyticSteadyState(),
    detect = function( pfimproject ) {
      f = feats( pfimproject ); !f$isODE && f$hasTau
    }
  )
  # Catch-all analytic: no detect - leaves room for custom detect() handlers.
  pfim_register_model_class(
    "ModelAnalytic",
    function() ModelAnalytic(),
    detect = NULL
  )
}

.pfimInitBuiltinModelRegistry()
