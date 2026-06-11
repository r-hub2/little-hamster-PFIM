#' @title LibraryOfModels
#' @description Registry of built-in PK and PD model definitions.
#' @param models List of PK and PD model metadata objects.
#' @export

LibraryOfModels = new_class("LibraryOfModels", package = "PFIM",

                       properties = list(
                         models = new_property( class_list, default = list())))

