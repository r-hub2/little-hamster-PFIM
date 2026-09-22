#' @title LibraryOfPKModels
#' @description
#' Built-in pharmacokinetic (PK) model library (bolus, infusion, absorption,
#' one- and two-compartment and Michaelis-Menten structures).
#'
#' Exported \code{Linear*} helpers return equation *strings* (not S7 objects):
#' two-compartment forms expand hybrid rates (alpha/beta) and partial-fraction
#' coefficients (A/B) via string substitution. Prefer looking up names in
#' \code{pkModelLibrary} for new code.
#' @param models Named list of built-in PK model equation definitions.
#' @return An S7 object of class \code{LibraryOfPKModels}.
#' @examples
#' \donttest{
#' names(prop(pkModelLibrary, "models"))[1:5]
#' nzchar(Linear2BolusSingleDose_ClQV1V2())
#' }
#' @include LibraryOfModels.R
#' @export

LibraryOfPKModels = new_class( "LibraryOfPKModels", package = "PFIM", parent = LibraryOfModels )

# Two-compartment closed forms share hybrid rates (alpha/beta) and A/B
# coefficients. Substitutions are wrapped so B/beta is never B/0.5 * (...).
.pfimEqParen = function( s ) str_c( "(", s, ")" )

.pfimEqFill = function( template, repl, wrap = TRUE ) {
  keys = names( repl )[ order( nchar( names( repl ) ), decreasing = TRUE ) ]
  out = template
  for ( k in keys ) {
    val = repl[[ k ]]
    if ( isTRUE( wrap ) ) val = .pfimEqParen( val )
    out = str_replace_all( out, paste0( "\\b", k, "\\b" ), val )
  }
  str_replace_all( out, "\\s+", " " )
}

.pfimTwoCptRates = function( param = c( "ClQ", "micro" ), wrap = TRUE ) {
  param = match.arg( param )
  if ( identical( param, "ClQ" ) ) {
    p = list(
      beta = "0.5*(Q/V1+Q/V2+Cl/V1-sqrt((Q/V1+Q/V2+Cl/V1)^2-4*Q/V2*Cl/V1))",
      k21  = "Q/V2",
      V    = "V1"
    )
    p$alpha = .pfimEqFill( "(Q/V2*Cl/V1)/beta", p, wrap = wrap )
  } else {
    p = list(
      beta = "0.5*(k12+k21+k-sqrt((k12+k21+k)^2-4*k21*k))",
      k21  = "k21",
      V    = "V"
    )
    p$alpha = .pfimEqFill( "(k21*k)/beta", p, wrap = wrap )
  }
  p
}

.pfimTwoCptAB = function( p, first_order = FALSE, wrap = TRUE ) {
  tpl = if ( isTRUE( first_order ) ) {
    list(
      A = "ka/V*(k21-alpha)/(beta-alpha)/(ka-alpha)",
      B = "ka/V*(k21-beta)/(alpha-beta)/(ka-beta)"
    )
  } else {
    list(
      A = "1/V*(alpha-k21)/(alpha-beta)",
      B = "1/V*(beta-k21)/(beta-alpha)"
    )
  }
  list(
    A = .pfimEqFill( tpl$A, p, wrap = wrap ),
    B = .pfimEqFill( tpl$B, p, wrap = wrap )
  )
}

.pfimLinear2Bolus = function( param, ss = FALSE ) {
  p = .pfimTwoCptRates( param )
  tpl = if ( isTRUE( ss ) ) {
    "dose_RespPK*(A*exp(-alpha*t)/(1-exp(-alpha*tau))+B*exp(-beta*t)/(1-exp(-beta*tau)))"
  } else {
    "dose_RespPK*(A*exp(-alpha*t)+B*exp(-beta*t))"
  }
  .pfimEqFill( tpl, c( .pfimTwoCptAB( p ), p[ c( "alpha", "beta" ) ] ) )
}

.pfimLinear2FirstOrder = function( param, ss = FALSE ) {
  p = .pfimTwoCptRates( param )
  tpl = if ( isTRUE( ss ) ) {
    "dose_RespPK*(A*exp(-alpha*t)/(1-exp(-alpha*tau))+B*exp(-beta*t)/(1-exp(-beta*tau))-(A+B)*exp(-ka*t)/(1-exp(-ka*tau)))"
  } else {
    "dose_RespPK*(A*exp(-alpha*t)+B*exp(-beta*t)-(A+B)*exp(-ka*t))"
  }
  .pfimEqFill( tpl, c( .pfimTwoCptAB( p, first_order = TRUE ), p[ c( "alpha", "beta" ) ] ) )
}

.pfimLinear2Infusion = function( param, ss = FALSE ) {
  p = .pfimTwoCptRates( param )
  ab = c( .pfimTwoCptAB( p ), p[ c( "alpha", "beta" ) ] )
  if ( isTRUE( ss ) ) {
    during = "dose_RespPK/Tinf_RespPK*(A/alpha*(1-exp(-alpha*t)+exp(-alpha*tau)*(1-exp(-alpha*Tinf_RespPK))*exp(-alpha*(t-Tinf_RespPK))/(1-exp(-alpha*tau)))+B/beta*(1-exp(-beta*t)+exp(-beta*tau)*(1-exp(-beta*Tinf_RespPK))*exp(-beta*(t-Tinf_RespPK))/(1-exp(-beta*tau))))"
    after  = "dose_RespPK/Tinf_RespPK*(A/alpha*(1-exp(-alpha*Tinf_RespPK))*exp(-alpha*(t-Tinf_RespPK))/(1-exp(-alpha*tau))+B/beta*(1-exp(-beta*Tinf_RespPK))*exp(-beta*(t-Tinf_RespPK))/(1-exp(-beta*tau)))"
  } else {
    during = "dose_RespPK/Tinf_RespPK*(A/alpha*(1-exp(-alpha*t))+B/beta*(1-exp(-beta*t)))"
    after  = "dose_RespPK/Tinf_RespPK*(A/alpha*(1-exp(-alpha*Tinf_RespPK))*exp(-alpha*(t-Tinf_RespPK))+B/beta*(1-exp(-beta*Tinf_RespPK))*exp(-beta*(t-Tinf_RespPK)))"
  }
  list(
    duringInfusion = list( "RespPK" = .pfimEqFill( during, ab ) ),
    afterInfusion  = list( "RespPK" = .pfimEqFill( after, ab ) )
  )
}

.pfimMm1BolusVmKm = list( "Deriv_C1" = "-Vm*C1/(Km+C1)" )

.pfimMm1InfusionVmKmV = list(
  duringInfusion = list( "Deriv_C1" = "-Vm*C1/(Km+C1) + dose_RespPK/V/Tinf_RespPK" ),
  afterInfusion  = list( "Deriv_C1" = "-Vm*C1/(Km+C1)" )
)

#' Analytic PK equation: Linear2BolusSingleDose_ClQV1V2.
#'
#' Two-compartment IV bolus (Cl, Q, V1, V2 parameterization).
#' @name Linear2BolusSingleDose_ClQV1V2
#' @return Character string: closed-form PK equation.
#' @export

Linear2BolusSingleDose_ClQV1V2 = function() .pfimLinear2Bolus( "ClQ" )

#' Analytic PK equation: Linear2BolusSingleDose_kk12k21V (micro-constant form).
#' @name Linear2BolusSingleDose_kk12k21V
#' @return Character string: closed-form PK equation.
#' @export

Linear2BolusSingleDose_kk12k21V = function() .pfimLinear2Bolus( "micro" )

#' Analytic PK equation: Linear2BolusSteadyState_ClQV1V2tau.
#' @name Linear2BolusSteadyState_ClQV1V2tau
#' @return Character string: closed-form PK equation.
#' @export

Linear2BolusSteadyState_ClQV1V2tau = function() .pfimLinear2Bolus( "ClQ", ss = TRUE )

#' Analytic PK equation: Linear2BolusSteadyState_kk12k21Vtau.
#' @name Linear2BolusSteadyState_kk12k21Vtau
#' @return Character string: closed-form PK equation.
#' @export

Linear2BolusSteadyState_kk12k21Vtau = function() .pfimLinear2Bolus( "micro", ss = TRUE )

#' Analytic PK equation: Linear2FirstOrderSingleDose_kaClQV1V2.
#' @name Linear2FirstOrderSingleDose_kaClQV1V2
#' @return Character string: closed-form PK equation.
#' @export

Linear2FirstOrderSingleDose_kaClQV1V2 = function() .pfimLinear2FirstOrder( "ClQ" )

#' Analytic PK equation: Linear2FirstOrderSingleDose_kakk12k21V.
#' @name Linear2FirstOrderSingleDose_kakk12k21V
#' @return Character string: closed-form PK equation.
#' @export

Linear2FirstOrderSingleDose_kakk12k21V = function() .pfimLinear2FirstOrder( "micro" )

#' Analytic PK equation: Linear2FirstOrderSteadyState_kaClQV1V2tau.
#' @name Linear2FirstOrderSteadyState_kaClQV1V2tau
#' @return Character string: closed-form PK equation.
#' @export

Linear2FirstOrderSteadyState_kaClQV1V2tau = function() .pfimLinear2FirstOrder( "ClQ", ss = TRUE )

#' Analytic PK equation: Linear2FirstOrderSteadyState_kakk12k21Vtau.
#' @name Linear2FirstOrderSteadyState_kakk12k21Vtau
#' @return Character string: closed-form PK equation.
#' @export

Linear2FirstOrderSteadyState_kakk12k21Vtau = function() .pfimLinear2FirstOrder( "micro", ss = TRUE )

#' Analytic PK equation: Linear2InfusionSingleDose_kk12k21V.
#' @name Linear2InfusionSingleDose_kk12k21V
#' @return Named list with \code{duringInfusion} / \code{afterInfusion} equation strings.
#' @export

Linear2InfusionSingleDose_kk12k21V = function() .pfimLinear2Infusion( "micro" )

#' Analytic PK equation: Linear2InfusionSingleDose_ClQV1V2.
#' @name Linear2InfusionSingleDose_ClQV1V2
#' @return Named list with \code{duringInfusion} / \code{afterInfusion} equation strings.
#' @export

Linear2InfusionSingleDose_ClQV1V2 = function() .pfimLinear2Infusion( "ClQ" )

#' Analytic PK equation: Linear2InfusionSteadyState_kk12k21Vtau.
#' @name Linear2InfusionSteadyState_kk12k21Vtau
#' @return Named list with \code{duringInfusion} / \code{afterInfusion} equation strings.
#' @export

Linear2InfusionSteadyState_kk12k21Vtau = function() .pfimLinear2Infusion( "micro", ss = TRUE )

#' Analytic PK equation: Linear2InfusionSteadyState_ClQV1V2tau.
#' @name Linear2InfusionSteadyState_ClQV1V2tau
#' @return Named list with \code{duringInfusion} / \code{afterInfusion} equation strings.
#' @export

Linear2InfusionSteadyState_ClQV1V2tau = function() .pfimLinear2Infusion( "ClQ", ss = TRUE )

# Named catalogue assembled into pkModelLibrary (analytic strings + ODE/MM forms).
models = list(
  # Linear

  # 1. One compartment

  # 1.1 IV bolus

  # 1.1.1 Single dose
  "Linear1BolusSingleDose_kV" = list( "RespPK" = "dose_RespPK/V*(exp(-k*t))" ) ,
  "Linear1BolusSingleDose_ClV" = list( "RespPK" = "dose_RespPK/V * (exp(-Cl/V* t))" ),

  # 1.1.2 Steady state
  "Linear1BolusSteadyState_ClVtau" = list( "RespPK" = "dose_RespPK/V * ( exp(-Cl/V*t)/(1-exp(-Cl/V*tau)))" ),
  "Linear1BolusSteadyState_kVtau" = list("RespPK" = " dose_RespPK/V * ( exp( -k*t )/( 1-exp( -k*tau ) ) )"),

  # 1.2 Infusion

  # 1.2.1 Single dose
  "Linear1InfusionSingleDose_ClV" = list( duringInfusion = list( "RespPK" = "dose_RespPK/Tinf_RespPK/Cl * (1 - exp(-Cl/V * t ) )" ) ,
                                          afterInfusion  = list( "RespPK" = "dose_RespPK/Tinf_RespPK/Cl * (1 - exp(-Cl/V * Tinf_RespPK)) * (exp(-Cl/V * (t - Tinf_RespPK)))")),

  "Linear1InfusionSingleDose_kV" = list( duringInfusion = list( "RespPK" = "dose_RespPK/Tinf_RespPK/(k*V) * (1 - exp(-k * t ) )" ) ,
                                         afterInfusion  = list( "RespPK" = "(dose_RespPK/Tinf_RespPK)/(k*V) * (1 - exp(-k * Tinf_RespPK)) * (exp(-k * (t - Tinf_RespPK)))") ),

  "Linear1InfusionSteadyState_kVtau" = list( duringInfusion = list( "RespPK" = "dose_RespPK/Tinf_RespPK/(k*V) * ( (1 - exp(-k * t)) + exp(-k*tau) * ( (1 - exp(-k*Tinf_RespPK)) * exp(-k*(t-Tinf_RespPK)) / (1-exp(-k*tau) ) ) )" ) ,
                                             afterInfusion  = list( "RespPK" = "dose_RespPK/Tinf_RespPK/(k*V) * ( (1 - exp(-k*Tinf_RespPK ) ) * exp(-k*(t-Tinf_RespPK)) / (1-exp(-k*tau ) ) )") ),

  "Linear1InfusionSteadyState_ClVtau" = list( duringInfusion = list( "RespPK" = "dose_RespPK/Tinf_RespPK/Cl * ( ( 1 - exp(-(Cl/V) * t)) + exp(-(Cl/V)*tau) * ( (1 - exp(-(Cl/V)*Tinf_RespPK)) * exp(-(Cl/V)*(t-Tinf_RespPK)) / (1-exp(-(Cl/V)*tau) ) ) )" ) ,
                                              afterInfusion  = list( "RespPK" = "dose_RespPK/Tinf_RespPK/Cl * ( ( 1 - exp(-(Cl/V)*Tinf_RespPK ) ) * exp(-(Cl/V)*(t-Tinf_RespPK)) / (1-exp(-(Cl/V)*tau ) ) )" ) ),


  # 1.3 First order

  # 1.3.1 Single dose
  "Linear1FirstOrderSingleDose_kaClV" = list( "RespPK" = "dose_RespPK/V * ka/(ka - Cl/V) * (exp(-Cl/V * t) - exp(-ka * t))" ),

  "Linear1FirstOrderSingleDose_kakV" = list( "RespPK" = "dose_RespPK/V * ka/(ka - k) * (exp(-k * t) - exp(-ka * t))" ),

  # 1.3.2 Steady state
  "Linear1FirstOrderSteadyState_kaClVtau" = list("RespPK" = "dose_RespPK/V * ka/(ka - Cl/V) * (exp(-Cl/V * t)/(1-exp(-Cl/V * tau)) - exp(-ka * t)/(1-exp(-ka * tau)))"),
  "Linear1FirstOrderSteadyState_kakVtau" = list("RespPK" = "dose_RespPK/V * ka/(ka - k) * (exp(-k * t)/(1-exp(-k * tau)) - exp(-ka * t)/(1-exp(-ka * tau)))"),

  # 2. Two compartments

  # 2.1 Bolus

  # 2.1.1 Single dose

  # Linear2BolusSingleDose_ClQV1V2
  "Linear2BolusSingleDose_ClQV1V2" = list( "RespPK" = Linear2BolusSingleDose_ClQV1V2() ),

  # Linear2BolusSingleDose_kk12k21V
  "Linear2BolusSingleDose_kk12k21V" = list( "RespPK" = Linear2BolusSingleDose_kk12k21V() ),

  # 2.1.2 Steady state

  # Linear2BolusSteadyState_ClQV1V2tau
  "Linear2BolusSteadyState_ClQV1V2tau" = list( "RespPK" = Linear2BolusSteadyState_ClQV1V2tau() ),

  # Linear2BolusSteadyState_kk12k21Vtau
  "Linear2BolusSteadyState_kk12k21Vtau" = list( "RespPK" = Linear2BolusSteadyState_kk12k21Vtau() ),

  # 2.2 First order

  # 2.2.1 Single dose

  # Linear2First orderSingleDose_kaClQV1V2
  "Linear2FirstOrderSingleDose_kaClQV1V2" = list( "RespPK" = Linear2FirstOrderSingleDose_kaClQV1V2() ),

  # Linear2First orderSingleDose_kakk12k21V
  "Linear2FirstOrderSingleDose_kakk12k21V" = list( "RespPK" = Linear2FirstOrderSingleDose_kakk12k21V() ),

  # 2.2.2 Steady state

  # Linear2FirstOrderSteadyState_kaClQV1V2tau
  "Linear2FirstOrderSteadyState_kaClQV1V2tau" = list( "RespPK" = Linear2FirstOrderSteadyState_kaClQV1V2tau() ),

  # Linear2First orderSteadyState_kakk12k21V
  "Linear2FirstOrderSteadyState_kakk12k21Vtau" = list( "RespPK" = Linear2FirstOrderSteadyState_kakk12k21Vtau() ),

  # 2.3 Infusion

  # 2.3.1 Single dose

  # Linear2InfusionSingleDose_kk12k21V
  "Linear2InfusionSingleDose_kk12k21V"  = Linear2InfusionSingleDose_kk12k21V() ,

  # Linear2InfusionSingleDose_ClQV1V2
  "Linear2InfusionSingleDose_ClQV1V2" = Linear2InfusionSingleDose_ClQV1V2(),

  # 2.3.2 Steady state

  # Linear2InfusionSteadyState_kk12k21Vtau
  "Linear2InfusionSteadyState_kk12k21Vtau" = Linear2InfusionSteadyState_kk12k21Vtau(),

  # Linear2InfusionSteadyState_ClQV1V2tau
  "Linear2InfusionSteadyState_ClQV1V2tau" = Linear2InfusionSteadyState_ClQV1V2tau(),

  # Michaelis-Menten elimination

  # 1. One compartment

  # 1.1 IV bolus (vignette name lists V for the IC D/V; Deriv_C1 does not use V)
  "MichaelisMenten1BolusSingleDose_VmKmV" = .pfimMm1BolusVmKm,
  "MichaelisMenten1BolusSingleDose_VmKm" = .pfimMm1BolusVmKm,


  # 1.2 First order
  "MichaelisMenten1FirstOrderSingleDose_kaVmKmV" = list("Deriv_C1" = "-Vm*C1/(Km+C1) + dose_RespPK/V*ka*exp(-ka*t)"),

  # 1.3 Infusion
  "MichaelisMenten1InfusionSingleDose_VmKmV" = .pfimMm1InfusionVmKmV,
  "MichaelisMenten1InfusionSingleDose_VmKmk12k21V1V2" = .pfimMm1InfusionVmKmV,
  # 2. Two compartments

  # 2.1 IV bolus (central volume is V1, matching the catalogue suffix)
  "MichaelisMenten2BolusSingleDose_VmKmk12k21V1V2" = list(
    "Deriv_C1" = "-Vm*C1/(Km+C1) - k12*C1 + k21*(V2/V1)*C2",
    "Deriv_C2" = "k12*(V1/V2)*C1 - k21*C2"
  ),

  # 2.2 First order
  "MichaelisMenten2FirstOrderSingleDose_kaVmKmk12k21V1V2" = list(
    "Deriv_C1" = "-Vm*C1/(Km+C1) - k12*C1 + k21*(V2/V1)*C2 + dose_RespPK/V1*ka*exp(-ka*t)",
    "Deriv_C2" = "k12*(V1/V2)*C1 - k21*C2"
  ),

  # 2.3 Infusion
  "MichaelisMenten2InfusionSingleDose_VmKmk12k21V1V2" = list(
    duringInfusion = list(
      "Deriv_C1" = "-Vm*C1/(Km+C1) - k12*C1 + k21*(V2/V1)*C2 + dose_RespPK/V1/Tinf_RespPK",
      "Deriv_C2" = "k12*(V1/V2)*C1 - k21*C2"
    ),
    afterInfusion = list(
      "Deriv_C1" = "-Vm*C1/(Km+C1) - k12*C1 + k21*(V2/V1)*C2",
      "Deriv_C2" = "k12*(V1/V2)*C1 - k21*C2"
    )
  )
)  # end named PK equation catalogue (analytic + ODE / MM forms)

#' Default built-in PK model library instance.
#' @format An S7 object of class \code{LibraryOfPKModels}.
#' @return A \code{LibraryOfPKModels} object.
#' @examples
#' \donttest{
#' length(prop(pkModelLibrary, "models"))
#' }
#' @seealso \code{\link{LibraryOfPKModels}} for the class constructor.
#' @export
pkModelLibrary = LibraryOfPKModels( models )








