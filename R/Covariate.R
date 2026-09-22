#' Check categorical labels, effects, and non-negative proportions.
#' @noRd
#' @keywords internal
.validateCategoricalCategories = function( categories, effects,
                                           proportions = NULL,
                                           sequences = NULL,
                                           proportionLabel = "categoriesProportions",
                                           covName = NULL ) {
  if ( !length( categories ) || any( is.na( categories ) | !nzchar( categories ) ) )
    .pfimStop( "CategoricalCovariate: categories must be non-empty strings." )
  if ( anyDuplicated( categories ) )
    .pfimStop( "CategoricalCovariate: categories must be unique." )
  if ( !is.null( proportions ) ) {
    if ( length( proportions ) != length( categories ) )
      .pfimStop( proportionLabel, " must have one value per category." )
    if ( any( !is.finite( proportions ) ) || any( proportions < 0 ) )
      .pfimStop( proportionLabel, " must be finite and non-negative." )
  }
  extra = setdiff( names( effects ), categories )
  if ( length( extra ) )
    .pfimStop(
      "CategoricalCovariate: effects name(s) not in categories: ",
      paste( extra, collapse = ", " ), "."
    )
  .pfimWarnZeroCovariateEffects( effects, categories, covName )
  if ( !is.null( sequences ) ) {
    used = unique( unlist( sequences, use.names = FALSE ) )
    bad = setdiff( used, categories )
    if ( length( bad ) )
      .pfimStop(
        "CategoricalCovariateWithIOV: sequence label(s) not in categories: ",
        paste( bad, collapse = ", " ), "."
      )
  }
  invisible( NULL )
}

#' Warn when a declared non-reference effect is exactly zero (beta stays in the FIM).
#' @noRd
#' @keywords internal
.pfimWarnZeroCovariateEffects = function( effects, categories, covName = NULL ) {
  if ( !length( effects ) || !length( categories ) )
    return( invisible( NULL ) )
  ref = categories[[ 1L ]]
  walk( names( effects ), function( cat ) {
    if ( identical( cat, ref ) ) return( NULL )
    eff = effects[[ cat ]]
    nms = names( eff )
    if ( is.null( nms ) || !length( eff ) ) return( NULL )
    zeros = nms[ !is.na( nms ) & nzchar( nms ) & is.finite( eff ) & eff == 0 ]
    if ( !length( zeros ) ) return( NULL )
    who = if ( .pfimIsNonEmptyScalar( covName ) )
      paste0( "covariate '", covName, "' category '", cat, "'" )
    else
      paste0( "category '", cat, "'" )
    .pfimWarn(
      who, ": effect(s) ", paste( zeros, collapse = ", " ),
      " are 0; the corresponding \u03b2 parameter(s) stay in the FIM ",
      "so SE under the null remains defined."
    )
  } )
  invisible( NULL )
}

#' Check that a proportion vector sums to 1 (within floating-point tolerance).
#' @param x Numeric vector.
#' @param label Name used in error messages.
#' @noRd
#' @keywords internal
.validateUnitProportions = function( x, label ) {
  if ( length( x ) == 0L ) return( invisible( NULL ) )
  s = sum( x )
  tol = sqrt( .Machine$double.eps ) * max( 1, length( x ) )
  if ( abs( s - 1 ) > tol )
    stop( sprintf( "%s must sum to 1 (got %.8g).", label, s ), call. = FALSE )
  invisible( NULL )
}

#' @description
#' The class \code{Covariate} is the base class for all covariate types.
#'
#' The constructor doubles as a smart factory: when \code{categories} is
#' supplied it dispatches to the correct concrete subclass:
#' \itemize{
#'   \item \code{sequences} + \code{sequencesProportions} -> \code{CategoricalCovariateWithIOV}
#'   \item \code{categoriesProportions}                   -> \code{CategoricalCovariate}
#'   \item otherwise                                      -> plain \code{Covariate}
#' }
#'
#' @title Covariate
#' @param name                  Character: covariate identifier.
#' @param effects               Named list of covariate effects per category.
#'   A declared value of 0 still creates a \eqn{\beta} column in the population
#'   FIM (SE under the null); a warning is issued. Omit the parameter to drop it.
#' @param categories            (Optional) Character vector of category labels.
#' @param categoriesProportions (Optional) Numeric vector summing to 1.
#' @param sequences             (Optional) List of category sequences (IOV).
#' @param sequencesProportions  (Optional) Numeric vector summing to 1.
#' @return An S7 object of class \code{Covariate}.
#' @export
Covariate = new_class( "Covariate", package = "PFIM",
  properties = list(
    name    = class_character,
    effects = class_list
  ),
  validator = function( self ) {
    n = prop( self, "name" )
    if ( !.pfimIsNonEmptyScalar( n ) )
      return( "Covariate: name must be a non-empty string." )
    NULL
  },
  constructor = function( name,
                          effects               = list(),
                          categories            = NULL,
                          categoriesProportions = NULL,
                          sequences             = NULL,
                          sequencesProportions  = NULL ) {

    # Factory branch: concrete subclasses validate proportions themselves.
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

    new_object( S7_object(), name = name, effects = effects )
  }
)

#' Return covariate effect vectors by category
#' @param covariate A \code{Covariate} object.
#' @param ... Optional method arguments.
#' @name getCovariateEffects
#' @return A numeric vector or list of covariate effects.
#' @keywords internal
getCovariateEffects    = new_generic( "getCovariateEffects",    "covariate" )

#' Build one effect vector for a selected category
#' @param covariate A \code{Covariate} object.
#' @param ... Optional method arguments.
#' @name createEffectVector
#' @return A numeric vector of covariate effects.
#' @keywords internal
createEffectVector     = new_generic( "createEffectVector",     "covariate" )

#' Return the reference category for a covariate
#' @param covariate A \code{Covariate} object.
#' @param ... Optional method arguments.
#' @name getCategoryOfReference
#' @return A character value naming the reference category.
#' @keywords internal
getCategoryOfReference = new_generic( "getCategoryOfReference", "covariate" )

#' Reference category for a covariate (first element of \code{categories}).
#' @name getCategoryOfReference
#' @keywords internal

method( getCategoryOfReference, Covariate ) = function( covariate )
  prop( covariate, "categories" )[[ 1L ]]

#' Overlay category-specific effects onto a zero (or base) effect vector.
#'
#' The reference category leaves \code{effectVector} unchanged; other categories
#' write their named \code{effects} entries into matching parameter slots.
#' @name createEffectVector
#' @keywords internal

method( createEffectVector, Covariate ) = function( covariate, category, effectVector ) {

  if ( category == getCategoryOfReference( covariate ) ||
       !category %in% names( prop( covariate, "effects" ) ) )
    return( effectVector )

  covEffects = prop( covariate, "effects" )
  effectVector[ names( covEffects[[ category ]] ) ] = covEffects[[ category ]]
  effectVector
}
