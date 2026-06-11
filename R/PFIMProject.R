#' @title PFIMProject
#' @description
#' Base S7 class for PFIM design evaluation and optimization projects.
#'
#' The number of study occasions (IOV) is **not** set on the project: it is inferred
#' when the model is built (\code{defineModelType()}) from IOV covariate sequences
#' and/or random IOV (\code{gamma > 0}) on \code{ModelParameter} objects, and stored
#' on the \code{Model} as \code{numberOfOccasions}.
#'
#' All \code{fimType} values support covariates and IOV; \code{population} adds
#' \code{omega} and \code{gamma} blocks to the FIM.
#' @param name Character string: project name.
#' @param modelEquations List of model equations (or empty if using the model library).
#' @param modelFromLibrary List selecting a built-in PK/PD model.
#' @param modelParameters List of \code{ModelParameter} objects.
#' @param modelCovariates List of covariate objects (from \code{Covariate()} factory).
#' @param modelCovariatesEquation Character: \code{"additive"} or \code{"exponential"}.
#' @param modelError List of residual error model objects.
#' @param optimizer Character: optimization algorithm name (for \code{Optimization}).
#' @param optimizerParameters List of algorithm-specific settings.
#' @param outputs Named list mapping internal to user output names.
#' @param designs List of \code{Design} objects.
#' @param fimType Character: \code{"population"}, \code{"individual"}, or \code{"Bayesian"}.
#' @param fim \code{Fim} object filled after \code{run()}.
#' @param odeSolverParameters List with \code{atol} and \code{rtol} for ODE solvers.
#' @include Fim.R
#' @include CovariateModelEquation.R
#' @include pfim-project-access.R
#' @include model-type-dispatch.R
#' @export

PFIMProject = new_class("PFIMProject", package = "PFIM",
                        properties = list(
                          name = new_property(class_character, default = character(0)),
                          modelEquations = new_property(class_list, default = list()),
                          modelCovariatesEquation = new_property(class_character, default = character(0)),
                          modelFromLibrary = new_property(class_list, default = list()),
                          modelParameters = new_property(class_list, default = list()),
                          modelCovariates = new_property(class_list, default = list()),
                          modelError = new_property(class_list, default = list()),
                          optimizer = new_property(class_character, default = character(0)),
                          optimizerParameters = new_property(class_list, default = list()),
                          outputs = new_property(class_list, default = list()),
                          designs = new_property(class_list, default = list()),
                          fimType = new_property(class_character, default = character(0)),
                          fim = new_property(Fim, default = NULL),
                          odeSolverParameters = new_property(class_list, default = list())
                        ))

run = new_generic( "run", "pfimproject" )
defineFim = new_generic( "defineFim", c( "pfimproject" ) )
plotEvaluation = new_generic( "plotEvaluation", c( "pfimproject" ) )
plotSensitivityIndices = new_generic( "plotSensitivityIndices", c( "pfimproject" ) )
plotSE = new_generic( "plotSE", c( "pfimproject" ) )
plotRSE = new_generic( "plotRSE", c( "pfimproject" ) )

# Use S7's external-generic mechanism to wrap methods::show.
# This keeps S7 dispatch while preserving full S4 compatibility (print, autoprint).
# Do NOT create a new_generic("show", ...) - that would shadow methods::show and
# break print() and REPL auto-display.
#' Show methods for PFIM objects
#'
#' Console display for PFIM S7 classes. Dispatches via \code{methods::show}.
#'
#' @param object A PFIM object (\code{Evaluation}, \code{Optimization},
#'   \code{CovariateTest}, etc.).
#' @usage show(object)
#' @docType methods
#' @format NULL
#' @keywords methods
#' @name show-methods
#' @aliases show show-methods show,PFIM::Optimization-method
#'   show,PFIM::Evaluation-method show,PFIM::CovariateTest-method
#' @section Methods:
#' \describe{
#'   \item{\code{Optimization}}{Display initial and optimal design results.}
#'   \item{\code{Evaluation}}{Display evaluation results (FIM, SE, RSE, etc.).}
#'   \item{\code{CovariateTest}}{Display covariate-test results.}
#' }
#' @export
show = S7::new_external_generic( "methods", "show", "object" )

Report = new_generic( "Report", c( "pfimproject" ) )

getFisherMatrix = new_generic( "getFisherMatrix", c( "pfimproject" ) )
getSE = new_generic( "getSE", c( "pfimproject" ) )
getRSE = new_generic( "getRSE", c( "pfimproject" ) )
getShrinkage = new_generic( "getShrinkage", c( "pfimproject" ) )
getDeterminant = new_generic( "getDeterminant", c( "pfimproject" ) )
getDcriterion = new_generic( "getDcriterion", c( "pfimproject" ) )
getCorrelationMatrix = new_generic( "getCorrelationMatrix", c( "pfimproject" ) )

defineModelType = new_generic( "defineModelType", c( "pfimproject" ) )
defineModelEquationsFromLibraryOfModel = new_generic( "defineModelEquationsFromLibraryOfModel", c( "pfimproject" ) )

# tost is a plain alias for covariateTest (defined in CovariateTest.R).
# Do NOT declare it as new_generic here - CovariateTest.R would overwrite the
# generic with a plain function, destroying S7 dispatch.

#' Build the FIM object from project settings
#' @name defineFim
#' @export

method( defineFim, PFIMProject ) = function( pfimproject )
{
  raw = projectProp( pfimproject, "fimType" )
  if ( length( raw ) != 1L || !is.character( raw ) )
    stop( "fimType must be a single character string.", call. = FALSE )
  fimType = tolower( raw )
  switch(
    fimType,
    population = PopulationFim(),
    individual = IndividualFim(),
    bayesian   = BayesianFim(),
    stop(
      sprintf(
        "Invalid fimType '%s'. Use 'population', 'individual', or 'Bayesian'.",
        raw
      ),
      call. = FALSE
    )
  )
}

#' Infer the model class from equations and dosing
#' @name defineModelType
#' @export

method( defineModelType, PFIMProject ) = function( pfimproject )
{
  equations       = projectProp( pfimproject, "modelEquations" )
  parameters      = projectProp( pfimproject, "modelParameters" )
  modelCovariates = projectProp( pfimproject, "modelCovariates" )

  feats = .detectModelFeatures(
    equations,
    .initialConditionsFromProject( pfimproject )
  )
  model = .instantiateModelClass( .selectModelClass( feats ) )

  prop( model, "modelParameters") = parameters
  prop( model, "modelCovariates") = modelCovariates
  prop( model, "odeSolverParameters") = projectProp( pfimproject, "odeSolverParameters" )
  prop( model, "modelError") = projectProp( pfimproject, "modelError" )
  prop( model, "modelEquations") = equations
  prop( model, "numberOfOccasions" ) = inferNumberOfOccasions( modelCovariates, parameters )

  modelCovariatesEquation = projectProp( pfimproject, "modelCovariatesEquation")

  if ( length( modelCovariatesEquation ) != 0L ) {
    modelCovariatesEquation = switch( modelCovariatesEquation,
                                      "additive" = Additive(),
                                      "exponential" = Exponential() )
    prop( model, "modelCovariatesEquation") = modelCovariatesEquation
  }

  model
}

#' Load model equations from the PK or PD library
#' @name defineModelEquationsFromLibraryOfModel
#' @export

method( defineModelEquationsFromLibraryOfModel, PFIMProject ) = function( pfimproject )
{
  equations = prop( pfimproject, "modelFromLibrary" )
  outputs = prop( pfimproject, "outputs" )

  pkModelName = equations[["PKModel"]]
  pdModelName = equations[["PDModel"]]

  pkModels = prop( pkModelLibrary, "models")
  prop( pfimproject, "modelEquations" ) = pkModels[[equations[["PKModel"]]]]

  pkModel = defineModelType( pfimproject )

  # pkpd model
  if ( !is.null( pdModelName ) ) {
    pdModels = prop( pdModelLibrary, "models")
    prop( pfimproject, "modelEquations" ) = pdModels[[equations[["PDModel"]]]]
    pdModel = defineModelType( pfimproject )
    pkpdModelEquations = definePKPDModel( pkModel, pdModel, pfimproject )
  }else{
    pkpdModelEquations = definePKModel( pkModel, pfimproject )
  }

  return( pkpdModelEquations )
}
