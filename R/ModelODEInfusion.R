#' @title ModelODEInfusion
#' @description ODE model with infusion inputs.
#'
#' Inherits from \code{ModelODE} (numerical nature) rather than
#' \code{ModelInfusion}, so \code{S7_inherits(x, ModelODE)} is TRUE for
#' infusion ODEs. Use \code{.pfimIsInfusionModel()} for the infusion capability.
#' @inheritParams ModelODE
#' @include ModelODE.R
#' @include ModelInfusion.R
#' @return An S7 object of class \code{ModelODEInfusion}.
#' @export

ModelODEInfusion = new_class( "ModelODEInfusion", package = "PFIM", parent = ModelODE )

#' Evaluate initial conditions for ODE infusion models.
#'
#' Delegates to \code{.evalInitialConditionsImpl}: substitutes typical values
#' into arm initial-condition expressions (numeric values pass through).
#' @return Numeric vector of evaluated initial conditions.
#' @name evaluateInitialConditions
#' @keywords internal
method( evaluateInitialConditions, ModelODEInfusion ) = function( model, arm ) {
  .evalInitialConditionsImpl( model, arm )
}

#' TRUE when \code{model} is integrated with \code{deSolve} (any \code{ModelODE}).
#' @noRd
#' @keywords internal
.pfimIsOdeModel = function( model ) {
  S7::S7_inherits( model, ModelODE )
}

#' TRUE when \code{model} uses infusion (analytic \code{ModelInfusion} or ODE infusion).
#'
#' S7 has single inheritance: ODE-infusion classes inherit \code{ModelODE}, not
#' \code{ModelInfusion}. This predicate is the polymorphic infusion check.
#' @noRd
#' @keywords internal
.pfimIsInfusionModel = function( model ) {
  S7::S7_inherits( model, ModelInfusion ) ||
    S7::S7_inherits( model, ModelODEInfusion )
}
