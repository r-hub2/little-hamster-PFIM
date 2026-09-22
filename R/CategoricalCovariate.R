#' @description
#' The class \code{CategoricalCovariate} represents a fixed categorical
#' covariate (constant across occasions): binary or multi-level group
#' membership (sex, genotype, centre, ...).
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
#' @return An S7 object of class \code{CategoricalCovariate}.
#' @export

CategoricalCovariate = new_class( "CategoricalCovariate",
  package = "PFIM",
  parent = Covariate,
  properties = list(
    categories            = class_character,
    categoriesProportions = class_double
  ),
  constructor = function( name, categories, categoriesProportions,
                          effects = list() ) {
    .validateUnitProportions( categoriesProportions, "categoriesProportions" )
    .validateCategoricalCategories(
      categories, effects, proportions = categoriesProportions, covName = name
    )
    new_object(
      Covariate( name = name, effects = effects ),
      categories            = categories,
      categoriesProportions = categoriesProportions
    )
  }
)

#' Effect vectors for every category of a fixed categorical covariate.
#'
#' Returns a named list (one entry per category) built by
#' \code{createEffectVector()} from a shared null/base vector.
#' @name getCovariateEffects
#' @keywords internal

method( getCovariateEffects, CategoricalCovariate ) = function( covariate, nullVector ) {
  cats = prop( covariate, "categories" )
  cats |>
    map( ~ createEffectVector( covariate, .x, nullVector ) ) |>
    set_names( cats )
}
