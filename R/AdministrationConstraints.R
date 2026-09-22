#' @title AdministrationConstraints
#' @description
#' Allowed dose values per outcome when optimizing or reporting designs.
#'
#' Used by \code{generateDosesCombination()} to build the discrete dose grid for
#' Fedorov-Wynn / multiplicative algorithms (and constraint reports).
#' @param outcome Character string: outcome name.
#' @param doses List of admissible dose levels for that outcome.
#' @return An S7 object of class \code{AdministrationConstraints}.
#' @export

AdministrationConstraints = new_class( "AdministrationConstraints",
                                       package = "PFIM",
                                       properties = list(
                                         outcome = new_property(class_character, default = character(0)),
                                         doses = new_property(class_list, default = list())
                                       ),
                                       validator = function( self ) {
                                         o = prop( self, "outcome" )
                                         if ( length( o ) > 1L )
                                           return( "AdministrationConstraints: outcome must be a single string." )
                                         if ( .pfimIsBlankScalar( o ) )
                                           return( "AdministrationConstraints: outcome must be a non-empty string." )
                                         doses = prop( self, "doses" )
                                         if ( !length( doses ) )
                                           return( NULL )
                                         bad = which( !vapply( doses, function( x ) {
                                           is.numeric( x ) && length( x ) && all( is.finite( x ) ) && all( x >= 0 )
                                         }, logical( 1L ) ) )
                                         if ( length( bad ) )
                                           return( paste0(
                                             "AdministrationConstraints: doses[[", bad[ 1L ],
                                             "]] must be finite non-negative numeric."
                                           ) )
                                         NULL
                                       })
