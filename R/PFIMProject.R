#' @title PFIMProject
#' @description
#' Base S7 class for PFIM design evaluation and optimization projects.
#'
#' Set \code{numberOfOccasions} on \code{Evaluation} or \code{Optimization} (optional:
#' \code{NA} infers from covariates and IOV). When set explicitly, it must match
#' \code{\link{inferNumberOfOccasions}}. The resolved value is stored on the
#' \code{Model} in \code{defineModelType()}.
#'
#' All \code{fimType} values support covariates and IOV; \code{population} adds
#' \code{omega} and \code{gamma} blocks to the FIM.
#' @param name Character string: project name.
#' @param modelEquations List of model equations (or empty if using the model library).
#' @param modelClass Model S7 class name; filled by \code{defineModelType()} when empty.
#'   Custom classes: \code{\link{pfim_register_model_class}}.
#' @param modelFromLibrary List selecting a built-in PK/PD model.
#' @param modelParameters List of \code{ModelParameter} objects.
#' @param modelCovariates List of covariate objects (from \code{Covariate()} factory).
#' @param modelCovariatesEquation Character: \code{"additive"} or \code{"exponential"}.
#' @param modelError List of residual error model objects.
#' @param optimizer Character: optimization algorithm name (for \code{Optimization}).
#' @param optimizerParameters List of algorithm-specific settings passed to
#'   \code{run(Optimization)}. Names and types are checked at construction for
#'   built-in optimizers; see \code{\link{MultiplicativeAlgorithm}},
#'   \code{\link{FedorovWynnAlgorithm}}, etc. Do not set these on the
#'   \code{*Algorithm} object - that object only stores outputs after \code{run()}.
#' @param outputs Named list mapping internal to user output names.
#' @param designs List of \code{Design} objects.
#' @param fimType Character: \code{"population"}, \code{"individual"}, or \code{"Bayesian"}.
#' @param fim \code{Fim} object filled after \code{run()}.
#' @param odeSolverParameters List with \code{atol} and \code{rtol} for ODE solvers.
#'   Finite-difference steps assume parameter mus are scaled to O(1).
#' @param numberOfOccasions Integer number of study occasions; \code{NA} to infer
#'   (see \code{\link{inferNumberOfOccasions}}). If set, must be consistent with
#'   IOV covariate sequences and \code{gamma} on parameters.
#' @param cacheScope Reserved (internal FIM cache id).
#' @include Fim.R
#' @include CovariateModelEquation.R
#' @include pfim-project-access.R
#' @include model-type-dispatch.R
#' @return An S7 object of class \code{PFIMProject}.
#' @export

PFIMProject = new_class("PFIMProject", package = "PFIM",
                        properties = list(
                          name = new_property(class_character, default = character(0)),
                          modelClass = new_property(class_character, default = character(0)),
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
                          fim = new_property(NULL | Fim, default = NULL),
                          odeSolverParameters = new_property(class_list, default = list()),
                          numberOfOccasions = new_property(class_numeric, default = NA_real_),
                          cacheScope = new_property( class_character, default = "" )
                        ),
                        validator = function( self ) {
                          checks = list(
                            list( prop( self, "designs" ), Design, "PFIMProject:designs", "Design" ),
                            list( prop( self, "modelParameters" ), ModelParameter,
                                  "PFIMProject:modelParameters", "ModelParameter" ),
                            list( prop( self, "modelCovariates" ), Covariate,
                                  "PFIMProject:modelCovariates", "Covariate" ),
                            list( prop( self, "modelError" ), ModelError,
                                  "PFIMProject:modelError", "ModelError" )
                          )
                          msgs = compact( map( checks, function( chk )
                            do.call( .validateS7List, chk ) ) )
                          if ( length( msgs ) )
                            return( msgs[[ 1L ]] )
                          NULL
                        })

#' Run a PFIM evaluation or optimization project
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name run
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' getDeterminant(ev)
#' }
#' @export
run = new_generic( "run", "pfimproject" )

#' Build the FIM object from project settings
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name defineFim
#' @export
defineFim = new_generic( "defineFim", c( "pfimproject" ) )

#' Plot evaluation outputs for a project
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Method-specific arguments. For \code{Evaluation} and
#'   \code{Optimization}, \code{plotOptions} list with \code{unitTime},
#'   \code{unitOutcomes}. Response curves are densified for display
#'   with sampling markers at design times.
#' @usage plotEvaluation(pfimproject, ...)
#' @return A nested list of \code{ggplot2} plot objects (by design / arm / outcome).
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' names(plotEvaluation(ev, list()))
#' }
#' @name plotEvaluation
#' @export
plotEvaluation = new_generic( "plotEvaluation", c( "pfimproject" ) )

#' Plot sensitivity indices for a project
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Method-specific arguments. For \code{Evaluation}, \code{plotOptions} list with
#'   \code{unitTime}, \code{unitOutcomes}. Gradients are plotted at design sampling times.
#' @usage plotSensitivityIndices(pfimproject, ...)
#' @return A \code{ggplot2} plot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' names(plotSensitivityIndices(ev, list()))
#' }
#' @name plotSensitivityIndices
#' @export
plotSensitivityIndices = new_generic( "plotSensitivityIndices", c( "pfimproject" ) )

#' Plot standard errors for project results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @usage plotSE(pfimproject, ...)
#' @return A \code{ggplot2} plot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' inherits(plotSE(ev), "ggplot")
#' }
#' @name plotSE
#' @export
plotSE = new_generic( "plotSE", c( "pfimproject" ) )

#' Plot relative standard errors for project results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @usage plotRSE(pfimproject, ...)
#' @return A \code{ggplot2} plot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' inherits(plotRSE(ev), "ggplot")
#' }
#' @name plotRSE
#' @export
plotRSE = new_generic( "plotRSE", c( "pfimproject" ) )

# methods::show via S7 external generic (NOT exported from PFIM).
# Use show(x) or methods::show(x) - never PFIM::show(x).
#' Show methods for PFIM objects
#'
#' Console display for PFIM S7 classes. Dispatches through the
#' \code{methods::show} generic (S7 external generic). \code{show} is
#' \strong{not} exported by PFIM: prefer \code{show(object)} after
#' \code{library(PFIM)} or \code{methods::show(object)}. Do not call
#' \code{PFIM::show(object)} (fragile / unbound).
#'
#' @param object A PFIM object (\code{Evaluation}, \code{Optimization},
#'   \code{CovariateTest}, etc.).
#' @usage show(object)
#' @docType methods
#' @format NULL
#' @keywords methods
#' @name show-methods
#' @aliases show-methods show,PFIM::Optimization-method show,PFIM::Evaluation-method show,PFIM::CovariateTest-method
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' show(ev)            # OK - methods::show
#' methods::show(ev)   # OK - explicit
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{Optimization}}{Display initial and optimal design results.}
#'   \item{\code{Evaluation}}{Display evaluation results (FIM, SE, RSE, etc.).}
#'   \item{\code{CovariateTest}}{Display covariate-test results.}
#' }
#' @keywords internal
#' @importFrom methods show
show = S7::new_external_generic( "methods", "show", "object" )

#' Generate an HTML evaluation or optimization report
#'
#' Requires \pkg{rmarkdown} and pandoc (\url{https://pandoc.org}).
#' @param pfimproject A \code{PFIMProject}, \code{Evaluation}, or \code{Optimization} object after \code{run()}.
#' @param ... \code{outputPath}, \code{outputFile}, \code{plotOptions} for HTML reports.
#' @name Report
#' @rdname Report
#' @return Invisibly, the path or result from the report generator.
#' @examples
#' \dontrun{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' tmp = tempfile(fileext = ".html")
#' Report(ev, outputPath = tempdir(), outputFile = basename(tmp), plotOptions = list())
#' }
#' @export
Report = new_generic( "Report", c( "pfimproject" ) )

#' Retrieve the Fisher information matrix from results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getFisherMatrix
#' @return A list with \code{fisherMatrix}, \code{fixedEffects},
#'   \code{varianceEffects}, and \code{singularFim}.
#' @export
getFisherMatrix = new_generic( "getFisherMatrix", c( "pfimproject" ) )

#' Retrieve standard errors from results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getSE
#' @return A list or matrix of standard errors for the FIM parameters.
#' @export
getSE = new_generic( "getSE", c( "pfimproject" ) )

#' Retrieve relative standard errors from results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getRSE
#' @return A list or matrix of relative standard errors (percent).
#' @export
getRSE = new_generic( "getRSE", c( "pfimproject" ) )

#' Retrieve shrinkage metrics from results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getShrinkage
#' @return Numeric vector or list of Bayesian shrinkage values.
#' @export
getShrinkage = new_generic( "getShrinkage", c( "pfimproject" ) )

#' Retrieve determinant of the Fisher matrix
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getDeterminant
#' @return Numeric determinant of the Fisher information matrix.
#' @export
getDeterminant = new_generic( "getDeterminant", c( "pfimproject" ) )

#' Retrieve D-criterion from project results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getDcriterion
#' @return Numeric D-optimality criterion value.
#' @export
getDcriterion = new_generic( "getDcriterion", c( "pfimproject" ) )

#' Retrieve correlation matrix from results
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name getCorrelationMatrix
#' @return Correlation matrix of the parameter estimates.
#' @export
getCorrelationMatrix = new_generic( "getCorrelationMatrix", c( "pfimproject" ) )

#' Infer and instantiate the concrete model class
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @name defineModelType
#' @keywords internal
defineModelType = new_generic( "defineModelType", c( "pfimproject" ) )

#' Build model equations from PK/PD model library entries
#' @param pfimproject A \code{PFIMProject} object.
#' @param ... Optional method arguments.
#' @return List of model equations resolved from the library.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' nrow(getFisherMatrix(ev)$fisherMatrix)
#' }
#' @name defineModelEquationsFromLibraryOfModel
#' @keywords internal
defineModelEquationsFromLibraryOfModel = new_generic( "defineModelEquationsFromLibraryOfModel", c( "pfimproject" ) )

#' Build the FIM object from project settings.
#'
#' Instantiates \code{PopulationFim}, \code{IndividualFim}, or \code{BayesianFim}
#' from the project's \code{fimType} string via the type registry.
#' @name defineFim
#' @export

method( defineFim, PFIMProject ) = function( pfimproject )
{
  raw = projectProp( pfimproject, "fimType" )
  if ( length( raw ) != 1L || !is.character( raw ) )
    .pfimStop( "fimType must be a single character string." )
  .pfimInstantiateFimType( raw )
}

#' Infer the model class from equations and dosing.
#'
#' Resolves \code{modelClass} (explicit or legacy token detection), clones
#' parameters/covariates onto a fresh model instance, then attaches residual
#' error, ODE tolerances, occasions, and the additive/exponential covariate
#' structural equation when requested.
#' @name defineModelType
#' @keywords internal

method( defineModelType, PFIMProject ) = function( pfimproject )
{
  equations       = projectProp( pfimproject, "modelEquations" )
  parameters      = projectProp( pfimproject, "modelParameters" )
  modelCovariates = projectProp( pfimproject, "modelCovariates" )

  classInfo = .pfimResolveModelClassInfo( pfimproject )
  className = classInfo$class
  .pfimSyncModelClass( pfimproject, className )
  if ( identical( classInfo$source, "legacy" ) &&
       isTRUE( pfim_get_option( "model.detect.legacy_warn", FALSE ) ) )
    warning(
      "Model class '", className,
      "' inferred from equation tokens; set modelClass explicitly.",
      call. = FALSE
    )
  model = .pfimInstantiateModelClass( className )

  # Clone so arm-level mutations during evaluation cannot mutate the project.
  prop( model, "modelParameters" ) = map( parameters, .pfimCloneS7 )
  prop( model, "modelCovariates" ) = map( modelCovariates, .pfimCloneS7 )
  prop( model, "odeSolverParameters") = projectProp( pfimproject, "odeSolverParameters" )
  prop( model, "modelError") = projectProp( pfimproject, "modelError" )
  prop( model, "modelEquations") = equations
  prop( model, "numberOfOccasions" ) = resolveNumberOfOccasions(
    projectProp( pfimproject, "numberOfOccasions" ),
    modelCovariates,
    parameters
  )

  modelCovariatesEquation = projectProp( pfimproject, "modelCovariatesEquation")

  # Character labels from the user API -> concrete Additive / Exponential objects.
  if ( length( modelCovariatesEquation ) != 0L ) {
    modelCovariatesEquation = switch( modelCovariatesEquation,
                                      "additive" = Additive(),
                                      "exponential" = Exponential() )
    prop( model, "modelCovariatesEquation") = modelCovariatesEquation
  }

  model
}

#' Load model equations from the PK or PD library.
#'
#' Looks up \code{modelFromLibrary$PKModel} (and optional \code{PDModel}),
#' temporarily writes each library equation onto the project to infer model
#' class, then returns the combined PK or PK/PD equation list.
#' @name defineModelEquationsFromLibraryOfModel
#' @keywords internal

method( defineModelEquationsFromLibraryOfModel, PFIMProject ) = function( pfimproject )
{
  equations = prop( pfimproject, "modelFromLibrary" )

  pkModelName = equations[["PKModel"]]
  pdModelName = equations[["PDModel"]]

  pkModels = prop( pkModelLibrary, "models" )
  prop( pfimproject, "modelEquations" ) = .pfimLibraryEquation( pkModels, pkModelName, "PKModel" )

  pkModel = defineModelType( pfimproject )

  # Optional PD library entry -> combined PK/PD equations; else PK-only.
  if ( !is.null( pdModelName ) ) {
    pdModels = prop( pdModelLibrary, "models" )
    prop( pfimproject, "modelEquations" ) = .pfimLibraryEquation( pdModels, pdModelName, "PDModel" )
    pdModel = defineModelType( pfimproject )
    pkpdModelEquations = definePKPDModel( pkModel, pdModel, pfimproject )
  }else{
    pkpdModelEquations = definePKModel( pkModel, pfimproject )
  }

  return( pkpdModelEquations )
}
