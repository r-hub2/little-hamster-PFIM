#' @title ModelInfusion
#' @description Base class for infusion administration models.
#' @inheritParams Model
#' @include Model.R
#' @export

ModelInfusion = new_class( "ModelInfusion", package = "PFIM", parent = Model )

#' Convert an analytic infusion model to ODE form
#' @name convertPKModelAnalyticToPKModelODE
#' @export

method( convertPKModelAnalyticToPKModelODE, ModelInfusion ) = function( pkModel ) {

  pkModelEquations = prop( pkModel, "modelEquations" )
  pkModelEquations = list( duringInfusion = pkModelEquations$duringInfusion,
                           afterInfusion  = pkModelEquations$afterInfusion )

  convertEquation = function( equation ) {
    dtEquationPKsubstitute = D( parse( text = equation ), "t" )
    dtEquationPKsubstitute = str_c( deparse( dtEquationPKsubstitute ), collapse = "" )

    if ( str_detect( equation, "Cl" ) ) {
      str_c( dtEquationPKsubstitute, "+(Cl/V)*", equation, "- (Cl/V)*RespPK" )
    } else {
      str_c( dtEquationPKsubstitute, "+k*", equation, "- k*RespPK" )
    }
  }

  pkModelEquations = map( pkModelEquations, ~ map( .x, convertEquation ) )

  return( pkModelEquations )
}
