#' @keywords internal
.equationsContainPattern = function( equations, pattern ) {
  map_lgl( equations, function( eq ) {
    if ( is.list( eq ) )
      any( map_lgl( eq, ~ grepl( pattern, .x, fixed = TRUE ) ) )
    else
      grepl( pattern, eq, fixed = TRUE )
  }) |> any()
}

#' @keywords internal
.detectModelFeatures = function( equations, initialConditions ) {
  leafNames = .getListLastName( equations )
  list(
    isODE                   = any( grepl( "Deriv_", leafNames, fixed = TRUE ) ),
    hasTau                  = .equationsContainPattern( equations, "tau" ),
    hasInfusion             = .equationsContainPattern( equations, "Tinf_" ),
    doseInEquation          = .equationsContainPattern( equations, "dose_" ),
    doseInInitialConditions = any( grepl( "dose_", initialConditions, fixed = TRUE ) )
  )
}

#' @keywords internal
.initialConditionsFromProject = function( pfimproject ) {
  map( projectProp( pfimproject, "designs" ), function( design ) {
    arms = prop( design, "arms" )
    names( arms ) = map_chr( arms, ~ prop( .x, "name" ) )
    map( arms, ~ prop( .x, "initialConditions" ) )
  }) |> unlist()
}

#' @keywords internal
.selectModelClass = function( feats ) {
  if ( feats$isODE ) {
    if ( feats$doseInInitialConditions )
      return( "ModelODEBolus" )
    if ( feats$hasInfusion && feats$doseInEquation )
      return( "ModelODEInfusionDoseInEquation" )
    if ( feats$doseInEquation )
      return( "ModelODEDoseInEquations" )
    return( "ModelODEDoseNotInEquations" )
  }
  if ( feats$hasInfusion && feats$doseInEquation && feats$hasTau )
    return( "ModelAnalyticInfusionSteadyState" )
  if ( feats$hasInfusion && feats$doseInEquation )
    return( "ModelAnalyticInfusion" )
  if ( feats$hasTau )
    return( "ModelAnalyticSteadyState" )
  "ModelAnalytic"
}

#' @keywords internal
.instantiateModelClass = function( className ) {
  switch(
    className,
    ModelODEBolus                    = ModelODEBolus(),
    ModelODEInfusionDoseInEquation   = ModelODEInfusionDoseInEquation(),
    ModelODEDoseInEquations          = ModelODEDoseInEquations(),
    ModelODEDoseNotInEquations       = ModelODEDoseNotInEquations(),
    ModelAnalyticInfusionSteadyState = ModelAnalyticInfusionSteadyState(),
    ModelAnalyticInfusion            = ModelAnalyticInfusion(),
    ModelAnalyticSteadyState         = ModelAnalyticSteadyState(),
    ModelAnalytic                    = ModelAnalytic(),
    stop( "Unknown model class: ", className, call. = FALSE )
  )
}
