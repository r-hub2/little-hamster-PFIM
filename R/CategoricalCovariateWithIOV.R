#' @description
#' The class \code{CategoricalCovariateWithIOV} represents a categorical
#' covariate that varies across occasions (inter-occasion variability design).
#'
#' Typical use: crossover treatment sequences (AB/BA), period effects, etc.
#'
#' Each element of \code{sequences} is a character vector of category labels,
#' one per occasion.  \code{sequencesProportions} gives the fraction of subjects
#' following each sequence and must sum to 1.
#'
#' @title CategoricalCovariateWithIOV
#' @param name                Character: covariate identifier.
#' @param categories          Character vector of category labels;
#'                            first element is the reference level.
#' @param sequences           Named list of character vectors (one per sequence
#'                            group); names are auto-generated as
#'                            \code{"sequence_1"}, \code{"sequence_2"}, ...
#' @param sequencesProportions Numeric vector summing to 1.
#' @param effects             Named list of covariate effects per category.
#' @include Covariate.R
#' @return An S7 object of class \code{CategoricalCovariateWithIOV}.
#' @export

CategoricalCovariateWithIOV = new_class( "CategoricalCovariateWithIOV",
  package = "PFIM",
  parent = Covariate,
  properties = list(
    categories           = class_character,
    sequences            = class_list,
    sequencesProportions = class_double
  ),
  constructor = function( name, categories, sequences, sequencesProportions,
                          effects = list() ) {

    .validateUnitProportions( sequencesProportions, "sequencesProportions" )
    .validateCategoricalCategories(
      categories, effects, sequences = sequences,
      proportionLabel = "sequencesProportions", covName = name
    )
    # Stable names for occasion nesting even when the user list was unnamed.
    names( sequences ) = paste0( "sequence_", seq_along( sequences ) )

    new_object(
      Covariate( name = name, effects = effects ),
      categories           = categories,
      sequences            = sequences,
      sequencesProportions = sequencesProportions
    )
  }
)

#' Effect vectors nested by sequence and occasion for IOV covariates.
#'
#' Each sequence yields a list of occasion-level effect vectors (same length as
#' the sequence), used when averaging the FIM over crossover designs.
#' @name getCovariateEffects
#' @keywords internal

method( getCovariateEffects, CategoricalCovariateWithIOV ) = function( covariate, effectVector ) {
  sequences = prop( covariate, "sequences" )
  nOcc      = length( pluck( sequences, 1L ) )
  occNames  = paste0( "occasion_", seq_len( nOcc ) )

  imap( sequences, function( seqVals, seqName )
    seq_len( nOcc ) |>
      map( ~ createEffectVector( covariate, seqVals[[ .x ]], effectVector ) ) |>
      set_names( occNames )
  )
}
