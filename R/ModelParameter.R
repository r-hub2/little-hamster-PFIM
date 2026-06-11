#' @title ModelParameter
#' @description
#' One population model parameter: distribution (mu, omega), optional IOV (gamma), and fix flags.
#' @param name          Character string: name of the parameter.
#' @param gamma         Numeric: the SD for inter-occasion variability (IOV).
#'                      \code{gamma^2} is the IOV variance component; \code{0} means no IOV.
#' @param distribution  An object of class \code{Distribution} for this parameter.
#' @param fixedMu       Logical: TRUE if \code{mu} is fixed (not estimated).
#' @param fixedOmega    Logical: TRUE if \code{omega} is fixed (not estimated).
#' @param value         Numeric: the value of the covariate effects.
#' @include Distribution.R
#' @export

ModelParameter = new_class( "ModelParameter",
                            package = "PFIM",
                            properties = list(
                              name         = new_property(class_character, default = character(0)),
                              gamma        = new_property(class_double,    default = 0),
                              distribution = new_property(Distribution,    default = NULL),
                              fixedMu      = new_property(class_logical,   default = FALSE),
                              fixedOmega   = new_property(class_logical,   default = FALSE),
                              value        = new_property(class_double,    default = numeric(0))
                            ))

getModelParametersData = new_generic( "getModelParametersData", c( "modelParameter" ) )

#' Parameter table for reports
#' @name getModelParametersData
#' @export

method( getModelParametersData, ModelParameter ) = function( modelParameter ) {

  list( modelParameter ) |>
    map( function( parameter ) {
      dist = prop( parameter, "distribution" )
      list(
        Parameters   = prop( parameter, "name"      ),
        mu           = as.character( prop( dist, "mu"    ) ),
        omega2       = as.character( prop( dist, "omega" )^2 ),
        gamma2       = as.character( prop( parameter, "gamma" )^2 ),
        Distribution = str_remove( class( dist )[1], "PFIM::" ),
        mu_fixed     = as.character( prop( parameter, "fixedMu"    ) ),
        omega2_fixed = as.character( prop( parameter, "fixedOmega" ) )
      )
    }) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
}
