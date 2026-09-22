#' @title ModelInfusion
#' @description Base class for infusion administration models.
#' @inheritParams Model
#' @include Model.R
#' @include ModelAnalytic.R
#' @return An S7 object of class \code{ModelInfusion}.
#' @export

ModelInfusion = new_class( "ModelInfusion", package = "PFIM", parent = Model )

#' Convert analytic during/after infusion formulas to ODE form.
#'
#' Applies \code{.convertAnalyticPkExprToOde} to every equation in the
#' \code{duringInfusion} and \code{afterInfusion} lists so infusion analytic
#' models can be remapped onto ODE library compartments.
#' @param pkModel An infusion analytic model with during/after equation lists.
#' @return List with \code{duringInfusion} and \code{afterInfusion} ODE strings.
#' @name convertPKModelAnalyticToPKModelODE
#' @keywords internal
method( convertPKModelAnalyticToPKModelODE, ModelInfusion ) = function( pkModel ) {
  pkEq = prop( pkModel, "modelEquations" )
  list(
    duringInfusion = map( pkEq$duringInfusion, .convertAnalyticPkExprToOde ),
    afterInfusion  = map( pkEq$afterInfusion,  .convertAnalyticPkExprToOde )
  )
}
