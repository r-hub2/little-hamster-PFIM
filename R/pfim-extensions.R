#' Extension listing helpers (model classes, FIM types, optimizers).
#'
#' Thin wrappers over the registries in \code{pfim-registry.R} and
#' \code{pfim-model-registry.R} for interactive discovery of registered plugins.
#' @name pfim-extensions
#' @include pfim-registry.R
#' @include pfim-model-registry.R
#' @keywords internal
NULL

#' Registered model class names
#' @return A character vector of registered class names.
#' @examples
#' \donttest{
#' pfim_registered_model_classes()
#' pfim_registered_fim_types()
#' pfim_registered_optimizers()
#' }
#' @export
pfim_registered_model_classes = function() {
  .pfimModelClassOrder()
}

#' @rdname pfim_registered_model_classes
#' @export
pfim_registered_fim_types = function() {
  ls( .pfimFimTypeRegistry, all.names = TRUE )
}

#' @rdname pfim_registered_model_classes
#' @export
pfim_registered_optimizers = function() {
  ls( .pfimOptimizerRegistry, all.names = TRUE )
}

#' List registered extensions
#' @return A list of registered extension names.
#' @export
pfim_list_extensions = function() {
  list(
    modelClass = pfim_registered_model_classes(),
    fimType = pfim_registered_fim_types(),
    optimizer = pfim_registered_optimizers()
  )
}
