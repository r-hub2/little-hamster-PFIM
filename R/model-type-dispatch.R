#' Match a whole-token pattern in model equation strings (Perl regex).
#'
#' Recurses into nested equation lists (e.g. during/after infusion).
#' @param equations Character vector or nested list of equation strings.
#' @param pattern Perl regular expression tested with \code{grepl}.
#' @return Logical scalar: \code{TRUE} if any equation matches \code{pattern}.
#' @noRd
#' @keywords internal
.equationsContainPattern = function( equations, pattern ) {
  map_lgl( equations, function( eq ) {
    if ( is.list( eq ) )
      any( map_lgl( eq, ~ grepl( pattern, .x, perl = TRUE ) ) )
    else
      grepl( pattern, eq, perl = TRUE )
  }) |> any()
}

#' Whether any administration in the project has a positive dosing-interval tau.
#'
#' Distinguishes steady-state dosing (\code{Administration(tau > 0)}) from a
#' model parameter coincidentally named \code{tau} that appears in equations.
#' @noRd
#' @keywords internal
.adminHasPositiveTau = function( pfimproject ) {
  designs = tryCatch( projectProp( pfimproject, "designs" ), error = function( e ) list() )
  if ( !length( designs ) ) return( FALSE )
  any( map_lgl( designs, function( design ) {
    any( map_lgl( prop( design, "arms" ), function( arm ) {
      any( map_lgl( prop( arm, "administrations" ), function( adm ) {
        tau = prop( adm, "tau" )
        length( tau ) == 1L && is.finite( tau ) && tau > 0
      } ) )
    } ) )
  } ) )
}

#' Detect structural features used to select the model S7 class.
#'
#' Inspects equation leaf names and tokens for ODE markers (\code{Deriv_}),
#' infusion duration (\code{Tinf_}), dose terms, and steady-state dosing.
#' Steady-state (\code{hasTau}) requires both a \code{\\btau\\b} token in the
#' equations \emph{and} \code{hasAdminTau = TRUE} (positive administration tau).
#' @param equations Named equation list from the project.
#' @param initialConditions Character/numeric initial-condition expressions.
#' @param hasAdminTau Logical; \code{TRUE} when an administration has \code{tau > 0}.
#' @return Named list of logical feature flags.
#' @noRd
#' @keywords internal
.detectModelFeatures = function( equations, initialConditions, hasAdminTau = FALSE ) {
  # Leaf names carry Deriv_* for ODE RHS keys.
  leafNames = .getListLastName( equations )
  hasTauToken = .equationsContainPattern( equations, "\\btau\\b" )
  list(
    isODE                   = any( grepl( "Deriv_", leafNames, fixed = TRUE ) ),
    # Parameter named tau in equations alone must not force SteadyState.
    hasTau                  = isTRUE( hasAdminTau ) && hasTauToken,
    hasTauToken             = hasTauToken,
    hasAdminTau             = isTRUE( hasAdminTau ),
    hasInfusion             = .equationsContainPattern( equations, "\\bTinf_" ),
    doseInEquation          = .equationsContainPattern( equations, "\\bdose_" ),
    doseInInitialConditions = any( grepl( "\\bdose_", initialConditions, perl = TRUE ) )
  )
}

#' Features for a project (equations + administration tau).
#' @noRd
#' @keywords internal
.detectModelFeaturesFromProject = function( pfimproject ) {
  .detectModelFeatures(
    projectProp( pfimproject, "modelEquations" ),
    .initialConditionsFromProject( pfimproject ),
    hasAdminTau = .adminHasPositiveTau( pfimproject )
  )
}

#' Flatten initial conditions from every arm in every design.
#' @param pfimproject A \code{PFIMProject} (or Evaluation/Optimization) object.
#' @return Unnamed vector of initial-condition expressions across arms.
#' @noRd
#' @keywords internal
.initialConditionsFromProject = function( pfimproject ) {
  map( projectProp( pfimproject, "designs" ), function( design ) {
    arms = prop( design, "arms" )
    names( arms ) = map_chr( arms, ~ prop( .x, "name" ) )
    map( arms, ~ prop( .x, "initialConditions" ) )
  }) |> unlist()
}

#' Map detected equation features to a concrete \code{Model*} class name.
#'
#' Priority: ODE bolus IC -> ODE infusion dose-in-eq -> ODE dose-in-eq ->
#' ODE dose-not-in-eq -> analytic infusion SS -> analytic infusion ->
#' analytic SS -> analytic bolus.
#' @param feats Feature list from \code{.detectModelFeatures()}.
#' @return Character scalar class name (without package prefix).
#' @noRd
#' @keywords internal
.selectModelClass = function( feats ) {
  if ( feats$isODE ) {
    # Dose in IC ⇒ bolus-at-t0 with deSolve events for later doses.
    if ( feats$doseInInitialConditions )
      return( "ModelODEBolus" )
    if ( feats$hasInfusion && feats$doseInEquation )
      return( "ModelODEInfusionDoseInEquation" )
    if ( feats$doseInEquation )
      return( "ModelODEDoseInEquations" )
    return( "ModelODEDoseNotInEquations" )
  }
  # Analytic branch: infusion + tau ⇒ steady state; infusion alone; tau alone; else bolus.
  if ( feats$hasInfusion && feats$doseInEquation && feats$hasTau )
    return( "ModelAnalyticInfusionSteadyState" )
  if ( feats$hasInfusion && feats$doseInEquation )
    return( "ModelAnalyticInfusion" )
  if ( feats$hasTau )
    return( "ModelAnalyticSteadyState" )
  "ModelAnalytic"
}
