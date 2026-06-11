#' @title AdministrationConstraints
#' @description
#' Allowed dose values per outcome when optimizing or reporting designs.
#' @param outcome Character string: outcome name.
#' @param doses List of admissible dose levels for that outcome.
#' @export

AdministrationConstraints = new_class( "AdministrationConstraints",
                                       package = "PFIM",
                                       properties = list(
                                         outcome = new_property(class_character, default = character(0)),
                                         doses = new_property(class_list, default = list())
                                       ))
