#' @title ModelODEInfusion
#' @description ODE model with infusion inputs.
#' @inheritParams ModelInfusion
#' @include ModelInfusion.R
#' @export

ModelODEInfusion = new_class( "ModelODEInfusion", package = "PFIM", parent = ModelInfusion )

method( evaluateInitialConditions, ModelODEInfusion ) = function( model, arm ) {
  .evalInitialConditionsImpl( model, arm )
}
