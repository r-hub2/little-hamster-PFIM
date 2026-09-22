#' @title LibraryOfPDModels
#' @description
#' Built-in pharmacodynamic (PD) model library (Emax, turnover, etc.).
#'
#' Entries are named lists of equation strings (\code{RespPD} for immediate
#' response, \code{Deriv_E} for turnover ODEs). Combined with a PK model via
#' \code{definePKPDModel()}.
#' @param models Named list of built-in PD model equation definitions.
#' @return An S7 object of class \code{LibraryOfPDModels}.
#' @examples
#' \donttest{
#' names(prop(pdModelLibrary, "models"))[1:5]
#' }
#' @include LibraryOfModels.R
#' @export

LibraryOfPDModels = new_class("LibraryOfPDModels",
                              package = "PFIM",
                              parent = LibraryOfModels )

# Named PD equation catalogue (immediate response + turnover Rin/kout families).
models = list(
  # Immediate Response Models

  # Drug action models

  # linear
  "ImmediateDrugLinear_S0Alin" = list( "RespPD" = "S0 + RespPK * Alin" ),
  # quadratic
  "ImmediateDrugImaxQuadratic_S0AlinAquad" = list( "RespPD" = "S0 + RespPK * Alin + Aquad * (RespPK)^2" ),
  # logarithmic
  "ImmediateDrugImaxLogarithmic_S0Alog" = list( "RespPD" = "S0 + Alog * log(RespPK)" ),
  # Emax
  "ImmediateDrugEmax_S0EmaxC50" = list( "RespPD" = "S0 + Emax*RespPK/(RespPK+C50)" ),
  # Sigmoid Emax
  "ImmediateDrugSigmoidEmax_S0EmaxC50gamma" = list( "RespPD" = "S0 + Emax*(RespPK^gamma)/(RespPK^gamma+C50^gamma)" ),
  # Imax
  "ImmediateDrugImax_S0ImaxC50" = list( "RespPD" = "S0 * (1 - Imax * RespPK/( RespPK + C50 ))" ),
  # Sigmoid Imax
  "ImmediateDrugImax_S0ImaxC50_gamma" = list( "RespPD" = "S0 * (1 - Imax * RespPK^gamma/(RespPK^gamma + C50^gamma ) )" ),
  # Full Imax / sigmoid full Imax (Imax = 1)
  "ImmediateDrugFullImax_S0C50" = list( "RespPD" = "S0 * (1 - RespPK/(RespPK + C50))" ),
  "ImmediateDrugSigmoidFullImax_S0C50gamma" = list( "RespPD" = "S0 * (1 - RespPK^gamma/(RespPK^gamma + C50^gamma))" ),

  # Baseline/disease models

  # Constant
  "ImmediateBaselineConstant_S0" = list( "RespPD" = "S0" ),
  # Linear
  "ImmediateBaselineLinear_S0kprog" = list( "RespPD" = "S0 + kprog*t" ),
  # Exponential disease increase (0 -> S0) / decrease (S0 -> 0)
  "ImmediateBaselineExponentialincrease_S0kprog" = list( "RespPD" = "S0*(1-exp(-kprog*t))" ),
  "ImmediateBaselineExponentialdecrease_S0kprog" = list( "RespPD" = "S0*exp(-kprog*t)" ),

  # Turnover Models

  # Models with impact on the input (Rin)

  # Emax
  "TurnoverRinEmax_RinEmaxCC50koutE" = list( "Deriv_E" = "Rin*(1+(Emax*RespPK)/(RespPK+C50))-kout*E" ),
  # Sigmoid Emax
  "TurnoverRinSigmoidEmax_RinEmaxCC50koutEgamma" = list( "Deriv_E" = "Rin*(1+(Emax*RespPK^gamma)/(RespPK^gamma+C50^gamma))-kout*E" ),
  # Imax
  "TurnoverRinImax_RinImaxCC50koutE" = list( "Deriv_E" = "Rin*(1-(Imax*RespPK)/(RespPK+C50))-kout*E" ),
  # Sigmoid Imax
  "TurnoverRinSigmoidImax_RinImaxCC50koutEgamma" = list( "Deriv_E" = "Rin*(1-(Imax*RespPK^gamma)/(RespPK^gamma+C50^gamma))-kout*E" ),
  # Full Imax
  "TurnoverRinFullImax_RinCC50koutE" = list( "Deriv_E" = "(Rin*(1-(RespPK)/(RespPK+C50))-kout*E)" ),
  # Sigmoid Full Imax
  "TurnoverRinSigmoidFullImax_RinCC50koutEgamma" = list( "Deriv_E" = "Rin*(1-RespPK^gamma/(RespPK^gamma+C50^gamma))-kout*E" ),

  # Models with impact on the output (kout)

  # Emax
  "TurnoverkoutEmax_RinEmaxCC50koutE" = list( "Deriv_E" = "Rin-kout*(1+(Emax*RespPK)/(RespPK+C50))*E" ),

  # Sigmoid Emax
  "TurnoverkoutSigmoidEmax_RinEmaxCC50koutEgamma" = list("Deriv_E" = "Rin-kout*(1+(Emax*RespPK^gamma)/(RespPK^gamma+C50^gamma))*E"),
  # Imax
  "TurnoverkoutImax_RinImaxCC50koutE" = list("Deriv_E" = "Rin-kout*(1-Imax*RespPK/(RespPK+C50))*E"),
  # Sigmoid Imax
  "TurnoverkoutSigmoidImax_RinImaxCC50koutEgamma" = list("Deriv_E" = "Rin-kout*(1-Imax*RespPK^gamma/(RespPK^gamma+C50^gamma))*E"),
  # Full Imax
  "TurnoverkoutFullImax_RinCC50koutE" = list("Deriv_E" = "Rin-kout*(1-RespPK/(RespPK+C50))*E"),
  # Sigmoid Full Imax
  "TurnoverkoutSigmoidFullImax_RinCC50koutEgamma" = list("Deriv_E" = "Rin-kout*(1-RespPK^gamma/(RespPK^gamma+C50^gamma))*E")
)

.pfimPdNameAlias = c(
  TurnoverRinSigmoidEmax_RinEmaxCC50koutE       = "TurnoverRinSigmoidEmax_RinEmaxCC50koutEgamma",
  TurnoverRinSigmoidImax_RinImaxCC50koutE       = "TurnoverRinSigmoidImax_RinImaxCC50koutEgamma",
  TurnoverRinSigmoidFullImax_RinCC50koutE       = "TurnoverRinSigmoidFullImax_RinCC50koutEgamma",
  TurnoverkoutSigmoidFullImax_RinCC50koutE      = "TurnoverkoutSigmoidFullImax_RinCC50koutEgamma"
)
models = c( models, stats::setNames( models[ unname( .pfimPdNameAlias ) ], names( .pfimPdNameAlias ) ) )

#' Default built-in PD model library instance.
#' @format An S7 object of class \code{LibraryOfPDModels}.
#' @return A \code{LibraryOfPDModels} object.
#' @examples
#' \donttest{
#' length(prop(pdModelLibrary, "models"))
#' }
#' @seealso \code{\link{LibraryOfPDModels}} for the class constructor.
#' @export
pdModelLibrary = LibraryOfPDModels( models )






