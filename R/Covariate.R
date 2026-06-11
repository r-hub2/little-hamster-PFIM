#' @description
#' The class \code{Covariate} is the base class for all covariate types.
#'
#' The constructor doubles as a smart factory: when \code{categories} is
#' supplied it dispatches to the correct concrete subclass:
#' \itemize{
#'   \item \code{sequences} + \code{sequencesProportions} â†’ \code{CategoricalCovariateWithIOV}
#'   \item \code{categoriesProportions}                   â†’ \code{CategoricalCovariate}
#'   \item otherwise                                      â†’ plain \code{Covariate}
#' }
#'
#' @title Covariate
#' @param name                  Character â€” covariate identifier.
#' @param effects               Named list of covariate effects per category.
#' @param categories            (Optional) Character vector of category labels.
#' @param categoriesProportions (Optional) Numeric vector summing to 1.
#' @param sequences             (Optional) List of category sequences (IOV).
#' @param sequencesProportions  (Optional) Numeric vector summing to 1.
#' @export

Covariate = new_class( "Covariate", package = "PFIM",
  properties = list(
    name    = class_character,
    effects = class_list
  ),
  constructor = function( name,
                          effects               = list(),
                          categories            = NULL,
                          categoriesProportions = NULL,
                          sequences             = NULL,
                          sequencesProportions  = NULL ) {

    if ( !is.null( categories ) ) {
      if ( !is.null( sequences ) && !is.null( sequencesProportions ) )
        return( CategoricalCovariateWithIOV(
          name                 = name,
          categories           = categories,
          sequences            = sequences,
          sequencesProportions = sequencesProportions,
          effects              = effects ) )

      if ( !is.null( categoriesProportions ) )
        return( CategoricalCovariate(
          name                  = name,
          categories            = categories,
          categoriesProportions = categoriesProportions,
          effects               = effects ) )
    }

    new_object( Covariate, name = name, effects = effects )
  }
)

getCovariateEffects    = new_generic( "getCovariateEffects",    "covariate" )
createEffectVector     = new_generic( "createEffectVector",     "covariate" )
getCategoryOfReference = new_generic( "getCategoryOfReference", "covariate" )

#' Reference category for a covariate
#' @name getCategoryOfReference
#' @export

method( getCategoryOfReference, Covariate ) = function( covariate )
  prop( covariate, "categories" )[[ 1L ]]

#' Covariate effect indicator vector
#' @name createEffectVector
#' @export

method( createEffectVector, Covariate ) = function( covariate, category, effectVector ) {

  if ( category == getCategoryOfReference( covariate ) ||
       !category %in% names( prop( covariate, "effects" ) ) )
    return( effectVector )

  covEffects = prop( covariate, "effects" )
  effectVector[ names( covEffects[[ category ]] ) ] = covEffects[[ category ]]
  effectVector
}
