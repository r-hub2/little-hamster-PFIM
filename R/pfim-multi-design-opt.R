#' @include Optimization.R
#' @include pfim-fim-cache.R
#' @keywords internal
NULL

# Multi-design optimization driver (continuous and discrete).
#
# When `length(designs) > 1`, run one independent search per design. The primary
# `optimisationDesign` / algorithm outputs stay on design 1; other designs are
# stored under `optimisationAlgorithmOutputs$perDesign`.

#' Run an optimizer once per design, then reassemble the Optimization.
#'
#' @param optimizationObject An \code{Optimization} with one or more designs.
#' @param optimizationAlgorithm Algorithm object (Multiplicative, FW, PSO, ...).
#' @param one_design_fn Function(\code{optimizationObject}, \code{algorithm}) that
#'   optimizes the single design currently in \code{projectProp(..., "designs")}.
#' @return Updated \code{Optimization} with all designs optimized.
#' @noRd
#' @keywords internal
.pfimOptimizeDesigns = function(
    optimizationObject,
    optimizationAlgorithm,
    one_design_fn,
    begin_cache = TRUE
) {
  designs = projectProp( optimizationObject, "designs" )
  # Single design: optionally open a FIM cache scope and delegate immediately.
  if ( length( designs ) <= 1L ) {
    if ( isTRUE( begin_cache ) )
      .pfimFimCacheBegin( optimizationObject )
    return( one_design_fn( optimizationObject, optimizationAlgorithm ) )
  }

  if ( isTRUE( begin_cache ) )
    .pfimFimCacheBegin( optimizationObject )
  state = reduce(
    seq_along( designs ),
    function( state, i ) {
      # Restrict the project to one design for this search.
      projectProp( state$obj, "designs" ) = list( designs[[ i ]] )
      obj = one_design_fn( state$obj, optimizationAlgorithm )
      opt_design = projectProp( obj, "designs" )[[ 1L ]]
      per = list(
        name               = prop( opt_design, "name" ),
        optimisationDesign = prop( obj, "optimisationDesign" ),
        optimalArms        = pluck(
          prop( obj, "optimisationAlgorithmOutputs" ),
          "optimalArms"
        )
      )
      primary = if ( i == 1L ) {
        list(
          optimisationDesign           = prop( obj, "optimisationDesign" ),
          optimisationAlgorithmOutputs = prop( obj, "optimisationAlgorithmOutputs" )
        )
      } else {
        state$primary
      }
      list(
        obj        = obj,
        optimized  = c( state$optimized, list( opt_design ) ),
        per_design = c( state$per_design, list( per ) ),
        primary    = primary
      )
    },
    .init = list(
      obj = optimizationObject,
      optimized = list(),
      per_design = list(),
      primary = NULL
    )
  )

  optimizationObject = state$obj
  projectProp( optimizationObject, "designs" ) = state$optimized
  prop( optimizationObject, "optimisationDesign" ) = state$primary$optimisationDesign
  outputs = state$primary$optimisationAlgorithmOutputs
  outputs$perDesign = state$per_design
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = outputs
  optimizationObject
}

#' Alias kept for continuous call sites (opens FIM cache in the driver).
#' @noRd
#' @keywords internal
.pfimContinuousOptimizeDesigns = function(
    optimizationObject,
    optimizationAlgorithm,
    one_design_fn
) {
  .pfimOptimizeDesigns(
    optimizationObject, optimizationAlgorithm, one_design_fn, begin_cache = TRUE
  )
}

#' Discrete call sites: \code{generateFimsFromConstraints} opens its own cache.
#' @noRd
#' @keywords internal
.pfimDiscreteOptimizeDesigns = function(
    optimizationObject,
    optimizationAlgorithm,
    one_design_fn
) {
  .pfimOptimizeDesigns(
    optimizationObject, optimizationAlgorithm, one_design_fn, begin_cache = FALSE
  )
}
