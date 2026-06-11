#' @description
#' The class \code{CategoricalCovariate} represents a fixed categorical
#' covariate (constant across occasions): binary or multi-level group
#' membership (sex, genotype, centre, â€¦).
#'
#' The constructor validates that \code{categoriesProportions} sums to 1
#' (within floating-point tolerance) before creating the object.
#'
#' @title CategoricalCovariate
#' @param name                  Covariate identifier string.
#' @param categories            Character vector of category labels.
#'                              The first element is the reference level.
#' @param categoriesProportions Numeric vector of proportions; must sum to 1.
#' @param effects               Named list of covariate effects per category.
#' @include Covariate.R
#' @export

CategoricalCovariate = new_class( "CategoricalCovariate",
  parent = Covariate,
  properties = list(
    categories            = class_character,
    categoriesProportions = class_double
  ),
  constructor = function( name, categories, categoriesProportions,
                          effects = list() ) {
    new_object( CategoricalCovariate,
                name                  = name,
                effects               = effects,
                categories            = categories,
                categoriesProportions = categoriesProportions )
  }
)

#' Estimated covariate effects on parameters
#' @name getCovariateEffects
#' @export

method( getCovariateEffects, CategoricalCovariate ) = function( covariate, nullVector ) {
  cats = prop( covariate, "categories" )
  cats |>
    map( ~ createEffectVector( covariate, .x, nullVector ) ) |>
    set_names( cats )
}
