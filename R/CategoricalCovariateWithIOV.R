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
#' @param name                Character â€” covariate identifier.
#' @param categories          Character vector of category labels;
#'                            first element is the reference level.
#' @param sequences           Named list of character vectors (one per sequence
#'                            group); names are auto-generated as
#'                            \code{"sequence_1"}, \code{"sequence_2"}, â€¦
#' @param sequencesProportions Numeric vector summing to 1.
#' @param effects             Named list of covariate effects per category.
#' @include Covariate.R
#' @export

CategoricalCovariateWithIOV = new_class( "CategoricalCovariateWithIOV",
  parent = Covariate,
  properties = list(
    categories           = class_character,
    sequences            = class_list,
    sequencesProportions = class_double
  ),
  constructor = function( name, categories, sequences, sequencesProportions,
                          effects = list() ) {

    names( sequences ) = paste0( "sequence_", seq_along( sequences ) )

    new_object( CategoricalCovariateWithIOV,
                name                 = name,
                effects              = effects,
                categories           = categories,
                sequences            = sequences,
                sequencesProportions = sequencesProportions )
  }
)

#' Estimated covariate effects on parameters
#' @name getCovariateEffects
#' @export

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
