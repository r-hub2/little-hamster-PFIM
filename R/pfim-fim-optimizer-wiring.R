# Wiring: setOptimalArms + generateReportOptimization methods
# Loaded after FIM classes and algorithm classes (see DESCRIPTION Collate).

#' Load-time \code{setOptimalArms} and report methods.
#'
#' Registers S7 methods that map each (FIM type x optimizer) pair to the helper
#' that builds optimal arms after \code{optimizeDesign()}, and attaches HTML
#' report templates per FIM type. Continuous optimizers (Simplex / PSO / PGBO)
#' share one path; discrete Multiplicative / Fedorov-Wynn use population-specific
#' helpers that allocate subjects from weights / frequencies.
#' @name pfim-fim-optimizer-wiring
#' @include Fim.R
#' @include pfim-fim-report-render.R
#' @include pfim-fim-optimal-arms.R
#' @include MultiplicativeAlgorithm.R
#' @include FedorovWynnAlgorithm.R
#' @include SimplexAlgorithm.R
#' @include PSOAlgorithm.R
#' @include PGBOAlgorithm.R
#' @include PopulationFim.R
#' @include IndividualFim.R
#' @include BayesianFim.R
#' @keywords internal
NULL

#' Optimal arms already stored on continuous optimizer outputs.
#' @param optimizationAlgorithm Simplex / PSO / PGBO after \code{optimizeDesign()}.
#' @return List of \code{Arm} objects.
#' @noRd
#' @keywords internal
.setOptimalArmsContinuous = function( optimizationAlgorithm, ... ) {
  prop( optimizationAlgorithm, "optimizerOutputs" )$optimalArms
}

#' Register one \code{setOptimalArms} method: (fimClass, algoClass) -> helper.
#' @noRd
#' @keywords internal
.pfimRegisterSetOptimalArms = function( fimClass, algoClass, helper ) {
  method( setOptimalArms, list( fimClass, algoClass ) ) =
    function( fim, optimizationAlgorithm, ... ) helper( optimizationAlgorithm, ... )
}

#' Register continuous \code{setOptimalArms} for all three FIM types.
#' @noRd
#' @keywords internal
.pfimRegisterSetOptimalArmsContinuous = function( algoClass ) {
  walk(
    list( PopulationFim, IndividualFim, BayesianFim ),
    ~ .pfimRegisterSetOptimalArms( .x, algoClass, .setOptimalArmsContinuous )
  )
}

# Discrete: all FIM types allocate subjects from weights / frequencies.
walk( list( IndividualFim, BayesianFim, PopulationFim ), function( .fim ) {
  .pfimRegisterSetOptimalArms( .fim, MultiplicativeAlgorithm, .setOptimalArmsMultiplicativePopulation )
  .pfimRegisterSetOptimalArms( .fim, FedorovWynnAlgorithm,    .setOptimalArmsFedorovWynnPopulation )
} )

walk(
  list( SimplexAlgorithm, PSOAlgorithm, PGBOAlgorithm ),
  .pfimRegisterSetOptimalArmsContinuous
)

.pfimRegisterOptimizationReports( PopulationFim, "Population" )
.pfimRegisterOptimizationReports( IndividualFim, "Individual" )
.pfimRegisterOptimizationReports( BayesianFim,   "Bayesian" )
