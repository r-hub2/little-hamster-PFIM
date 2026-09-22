#' Arm constraints for optimization HTML reports.
#'
#' S7 generic \code{getArmConstraints(arm, optimizationAlgorithm)} dispatches on
#' the optimizer class: discrete algorithms (Multiplicative, Fedorov-Wynn) use
#' fixed-time / dose tables; continuous ones (Simplex, PSO, PGBO) use sampling
#' windows. Shared row builders live in \code{.armConstraintsDiscrete} here and
#' \code{.armConstraintsContinuous} in \code{pfim-utils.R}.
#' @include Arm.R
#' @include MultiplicativeAlgorithm.R
#' @include FedorovWynnAlgorithm.R
#' @include SimplexAlgorithm.R
#' @include PSOAlgorithm.R
#' @include PGBOAlgorithm.R
#' @include pfim-utils.R
#' @title Arm constraints (reports)
#' @name pfim-arm-constraints
#' @keywords internal
NULL

#' Administration and sampling constraints for an optimizer
#' @param arm An \code{Arm} object.
#' @param optimizationAlgorithm An optimization algorithm object.
#' @param ... Not used.
#' @name getArmConstraints
#' @return A list of arm-level optimization constraints.
#' @keywords internal
getArmConstraints = new_generic( "getArmConstraints", c( "arm", "optimizationAlgorithm" ) )

#' Build discrete arm constraints table for reporting.
#'
#' One row per sampling-times constraint, with dose bounds looked up by outcome.
#' @param arm \code{Arm} object containing administration and sampling constraints.
#' @return List of named rows describing discrete optimization constraints.
#' @noRd
#' @keywords internal
.armConstraintsDiscrete = function( arm ) {
  armName = prop( arm, "name" )
  armSize = prop( arm, "size" )
  # Map outcome -> formatted dose list for the dose-constraints column.
  admins  = map( prop( arm, "administrationsConstraints" ), function( ac ) {
    outcome = prop( ac, "outcome" )
    doses   = paste0( "(", paste( unlist( prop( ac, "doses" ) ), collapse = ", " ), ")" )
    set_names( list( doses ), outcome )
  } ) |> flatten()
  map( prop( arm, "samplingTimesConstraints" ), function( sc ) {
    outcome = prop( sc, "outcome" )
    doseConstraints = admins[[ outcome ]]
    if ( is.null( doseConstraints ) ) doseConstraints = "."
    list(
      "Arms name"                       = armName,
      "Number of subjects"              = armSize,
      "Outcome"                         = outcome,
      "Initial samplings"               = paste0( "(", paste( prop( sc, "initialSamplings" ), collapse = ", " ), ")" ),
      "Fixed times"                     = paste0( "(", paste( prop( sc, "fixedTimes" ), collapse = ", " ), ")" ),
      "Number of samplings optimisable" = as.character( prop( sc, "numberOfsamplingsOptimisable" ) ),
      "Dose constraints"                = doseConstraints
    )
  } )
}

# Discrete optimizers share the same constraint table layout.
method( getArmConstraints, list( Arm, MultiplicativeAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsDiscrete( arm )
}

method( getArmConstraints, list( Arm, FedorovWynnAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsDiscrete( arm )
}

# Continuous optimizers report sampling windows (see .armConstraintsContinuous).
method( getArmConstraints, list( Arm, SimplexAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsContinuous( arm )
}

method( getArmConstraints, list( Arm, PSOAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsContinuous( arm )
}

method( getArmConstraints, list( Arm, PGBOAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  .armConstraintsContinuous( arm )
}
