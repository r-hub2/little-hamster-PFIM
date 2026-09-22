#' Constraint tables for optimization reports.
#'
#' Builds \code{kableExtra} tables of arm-level dose/sampling constraints for
#' discrete (Fedorov-Wynn / Multiplicative) and continuous (PSO / PGBO / Simplex)
#' optimizers. Uses \code{getArmConstraints()} from \code{pfim-arm-constraints.R}.
#' @include pfim-arm-constraints.R
#' @include pfim-utils.R
#' @name pfim-constraints-helpers
NULL

#' Build a kableExtra table from arm constraints (discrete optimizers).
#' @noRd
#' @keywords internal
.constraintsKbl = function( arms, optimizationAlgorithm, col_names ) {
  armsConstraints = map( pluck( arms, 1L ), ~ getArmConstraints( .x, optimizationAlgorithm ) )
  df = map( armsConstraints, ~ map( .x, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) ) |>
    list_flatten() |>
    list_rbind()
  colnames( df ) = col_names
  kbl( df, align = c( "l", rep( "c", ncol( df ) - 1L ) ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
}

#' kableExtra table of discrete-optimizer arm constraints for reports.
#' @param optimizationAlgorithm An optimization algorithm object (Simplex, Fedorov-Wynn, etc.).
#' @param arms List of \code{Arm} objects from the design.
#' @return A \code{kableExtra} table.
#' @noRd
#' @keywords internal
.pfimConstraintsTableDiscrete = function( optimizationAlgorithm, arms ) {
  .constraintsKbl(
    arms, optimizationAlgorithm,
    c( "Arms name", "Number of subjects", "Outcome",
       "Initial samplings", "Fixed times",
       "Number of samplings optimisable", "Dose constraints" )
  )
}

#' kableExtra table of continuous-optimizer arm constraints for reports.
#' @param optimizationAlgorithm A continuous optimizer object (PSO, PGBO, Simplex on windows).
#' @param arms List of \code{Arm} objects from the design.
#' @return A \code{kableExtra} table.
#' @noRd
#' @keywords internal
.pfimConstraintsTableContinuous = function( optimizationAlgorithm, arms ) {
  armsConstraints = map( pluck( arms, 1L ), ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = .constraintsArmsTable( armsConstraints )
  colnames( armsConstraints ) = c(
    "Arms name", "Number of subjects", "Outcome",
    "Initial samplings", "Samplings windows", "Number of times by windows", "Min sampling"
  )
  kbl( armsConstraints, align = c( "l", rep( "c", ncol( armsConstraints ) - 1L ) ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
}
