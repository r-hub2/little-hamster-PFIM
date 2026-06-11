#' @title CovariateSamplingMethods
#' @description
#' The class \code{CovariateSamplingMethods} specifies how covariate values
#' are drawn for Monte Carlo integration over the covariate distribution.
#'
#' Three approaches are supported:
#' \enumerate{
#'   \item \strong{Empirical}: supply a \code{data.frame} of observed covariate
#'     values and draw from it with or without replacement.
#'   \item \strong{Parametric}: supply a named list of marginal
#'     \code{distributions} and sample from each independently.
#'   \item \strong{Copula}: reserved for future joint-distribution modelling
#'     (slot currently inactive — copula support not yet implemented).
#' }
#' @param data          A \code{data.frame} of observed covariate values.
#' @param distributions Named list of marginal distribution objects.
#' @param replace       Logical — sample from \code{data} with replacement.
#' @param MCSamples     Integer — number of Monte Carlo samples to draw.
#' @name CovariateSamplingMethods
#' @export

CovariateSamplingMethods = new_class( "CovariateSamplingMethods", package = "PFIM",
  properties = list(
    data          = new_property( class_data.frame, default = data.frame() ),
    distributions = new_property( class_list,       default = list()       ),
    replace       = new_property( class_logical,    default = FALSE        ),
    MCSamples     = new_property( class_double,     default = numeric( 0 ) )
  ),
  constructor = function( data = data.frame(),
                          distributions = list(),
                          replace = FALSE,
                          MCSamples = numeric( 0 ) ) {
    new_object( .parent = S7_object(),
                data          = data,
                distributions = distributions,
                replace       = replace,
                MCSamples     = MCSamples )
  }
)
