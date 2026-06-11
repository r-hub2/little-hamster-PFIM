#' @title Model
#' @description
#' Internal model object built from a \code{PFIMProject}: parameters, covariates,
#' equations, and evaluation hooks (gradients, variance).
#' @param name Model name.
#' @param modelParameters List of \code{ModelParameter} objects.
#' @param modelParametersWithCovariates Per-combination/occasion parameters (internal).
#' @param modelCovariatesEquation \code{Additive} or \code{Exponential} covariate link.
#' @param omegaWithIOV IIV/IOV variance vector for population FIM (internal).
#' @param modelCovariates List of covariate objects.
#' @param covariatesCombination Covariate combination table (internal).
#' @param covariatesEffect Nested covariate effect structure (internal).
#' @param samplings Sampling times for the current arm.
#' @param modelEquations Model equations (analytic or ODE).
#' @param wrapper Compiled or R function wrapper for predictions.
#' @param outputFormula Output expressions per outcome (character strings).
#' @param outputNames Character vector of outcome names.
#' @param variableNames State variable names.
#' @param outcomesWithAdministration Outcomes linked to dosing.
#' @param outcomesWithNoAdministration Outcomes without dosing.
#' @param modelError List of residual error models.
#' @param odeSolverParameters \code{atol} and \code{rtol} for \code{deSolve}.
#' @param parametersForComputingGradient Parameters used in finite-difference Hessians.
#' @param initialConditions ODE initial conditions.
#' @param functionArguments Names passed to the model function.
#' @param functionArgumentsSymbol Symbolic argument list for parsing.
#' @param numberOfOccasions Inferred number of study occasions (stored on the model after
#'   \code{defineModelType()}); see \code{\link{inferNumberOfOccasions}}. Not a user argument
#'   on \code{Evaluation} or \code{Optimization}.
#' @include CovariateModelEquation.R
#' @export

Model = new_class( "Model", package = "PFIM",

                   properties = list(
                     name                           = new_property( class_character,        default = character(0) ),
                     modelParameters                = new_property( class_list,             default = list()       ),
                     modelCovariatesEquation        = new_property( CovariateModelEquation, default = NULL         ),
                     modelCovariates                = new_property( class_list,             default = list()       ),
                     covariatesEffect               = new_property( class_list,             default = list()       ),
                     covariatesCombination          = new_property( class_list,             default = list()       ),
                     modelParametersWithCovariates  = new_property( class_list,             default = list()       ),
                     numberOfOccasions              = new_property( class_numeric,          default = 1            ),
                     omegaWithIOV                   = new_property( class_double,           default = numeric(0)   ),
                     samplings                      = new_property( class_numeric,          default = numeric(0)   ),
                     modelEquations                 = new_property( class_list,             default = list()       ),
                     wrapper                        = new_property( class_function,         default = NULL         ),
                     outputFormula                  = new_property( class_list,             default = list()       ),
                     outputNames                    = new_property( class_character,        default = character(0) ),
                     variableNames                  = new_property( class_character,        default = character(0) ),
                     outcomesWithAdministration     = new_property( class_character,        default = character(0) ),
                     outcomesWithNoAdministration   = new_property( class_character,        default = character(0) ),
                     modelError                     = new_property( class_list,             default = list()       ),
                     odeSolverParameters            = new_property( class_list,             default = list()       ),
                     parametersForComputingGradient = new_property( class_list,             default = list()       ),
                     initialConditions              = new_property( class_double,           default = numeric(0)   ),
                     functionArguments              = new_property( class_character,        default = character(0) ),
                     functionArgumentsSymbol        = new_property( class_list,             default = list()       )
                   ) )

#' Compile model equations and output mapping before evaluation
#' @param model A \code{Model} object.
#' @param ... \code{evaluation} or \code{arm}, depending on the method.
#' @name defineModelWrapper
#' @export
defineModelWrapper                = new_generic( "defineModelWrapper",                c( "model" ) )

#' Bind dosing schedules and sampling grids to the model
#' @param model A \code{Model} object.
#' @param ... \code{arm} for model methods.
#' @name defineModelAdministration
#' @export
defineModelAdministration         = new_generic( "defineModelAdministration",         c( "model" ) )

#' Predict model responses at arm sampling times
#' @param model A \code{Model} object.
#' @param ... \code{arm} for model methods.
#' @name evaluateModel
#' @export
evaluateModel                     = new_generic( "evaluateModel",                     c( "model" ) )

#' Gradients of model responses with respect to parameters
#' @param model A \code{Model} object.
#' @param ... \code{arm} for model methods.
#' @name evaluateModelGradient
#' @export
evaluateModelGradient             = new_generic( "evaluateModelGradient",             c( "model" ) )

#' Residual variance at sampling times for an arm
#' @param model A \code{Model} object.
#' @param ... \code{arm} for model methods.
#' @name evaluateModelVariance
#' @export
evaluateModelVariance             = new_generic( "evaluateModelVariance",             c( "model" ) )

#' Evaluate ODE initial conditions from the arm
#' @param model A \code{Model} object.
#' @param ... \code{arm} and optionally \code{doseEvent} for bolus ODE models.
#' @name evaluateInitialConditions
#' @export
evaluateInitialConditions         = new_generic( "evaluateInitialConditions",         c( "model" ) )

#' Finite-difference perturbations for gradient computation
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name finiteDifferenceHessian
#' @export
finiteDifferenceHessian           = new_generic( "finiteDifferenceHessian",           c( "model" ) )
evaluateCovariatesEffects         = new_generic( "evaluateCovariatesEffects",         c( "model" ) )
evaluateOmegaMatrixFromCovariates = new_generic( "evaluateOmegaMatrixFromCovariates", c( "model" ) )
generateCovariatesCombination     = new_generic( "generateCovariatesCombination",     c( "model" ) )

#' Remap library PK equations onto project compartments
#' @param pkModel Library PK model object.
#' @param pfimproject A \code{PFIMProject}, \code{Evaluation}, or \code{Optimization} object.
#' @param ... Optional method arguments.
#' @name definePKModel
#' @export
definePKModel                     = new_generic( "definePKModel",   c( "pkModel", "pfimproject" ) )

#' Combine PK and PD library equations for the project
#' @param pkModel Library PK model object.
#' @param pdModel Library PD model object.
#' @param pfimproject A \code{PFIMProject}, \code{Evaluation}, or \code{Optimization} object.
#' @param ... Optional method arguments.
#' @name definePKPDModel
#' @export
definePKPDModel                   = new_generic( "definePKPDModel", c( "pkModel", "pdModel", "pfimproject" ) )
modelParametersWithCovariates     = new_generic( "modelParametersWithCovariates",     c( "model" ) )
defineCovariatesData              = new_generic( "defineCovariatesData",              c( "model" ) )
hasCovariates                     = new_generic( "hasCovariates",                     c( "model" ) )
evaluateModelWithCovariates       = new_generic( "evaluateModelWithCovariates",       c( "model" ) )


# Integer codes for the covariate-equation type (avoids magic numbers throughout).
modelCovEqExponential = 1L   # theta = mu * exp(beta * cov)
modelCovEqAdditive    = 2L   # theta = mu * (1 + beta * cov)

# Two-character separator used to build / parse reversible combination names.
# Must not appear in any covariate name or category label.
combinationSep = "::"

#' Extract the number of occasions from IOV covariate sequences.
#'
#' Reads sequence lengths from every \code{CategoricalCovariateWithIOV} and
#' requires them to match.
#'
#' @param modelCovariates List of covariate objects.
#' @return Integer occasion count (1 if no IOV covariate is present).
#' @keywords internal
getOccasionsFromIOVCovariates = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L ) return( 1L )

  covariateWithIov = .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )
  if ( length( covariateWithIov ) == 0L ) return( 1L )

  sequenceLengths = covariateWithIov |>
    map( ~ map_int( prop( .x, "sequences" ), length ) ) |>
    list_c() |>
    unique()

  if ( length( sequenceLengths ) > 1L )
    stop( "All occasion-based covariate sequences must use the same number of occasions." )

  as.integer( sequenceLengths[[1L]] )
}

#' Infer the number of study occasions for a model.
#'
#' Rules (in order):
#' \enumerate{
#'   \item If any \code{CategoricalCovariateWithIOV} is present, use the common
#'         sequence length (e.g. 2, 3, or 4 periods).
#'   \item Else if any parameter has \code{gamma > 0} (random IOV), use 2 occasions.
#'   \item Else use 1 occasion.
#' }
#'
#' @param modelCovariates List of covariate objects.
#' @param modelParameters List of \code{ModelParameter} objects.
#' @return Integer number of occasions.
#' @keywords internal
inferNumberOfOccasions = function( modelCovariates = list(), modelParameters = list() ) {
  has_occasion_cov = length(
    .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )
  ) > 0L

  if ( has_occasion_cov )
    return( getOccasionsFromIOVCovariates( modelCovariates ) )

  gamma_vals = map_dbl( modelParameters, function( p ) {
    g = prop( p, "gamma" )
    if ( length( g ) == 0L || is.na( g ) ) 0 else g
  })
  if ( any( gamma_vals > 0, na.rm = TRUE ) ) return( 2L )

  1L
}

#' Occasion count for a built \code{Model} object.
#'
#' Uses the value stored on the model when set (e.g. temporary models in the
#' gradient core force \code{numberOfOccasions = 1}); otherwise calls
#' \code{inferNumberOfOccasions()}.
#'
#' @param model A \code{Model} object.
#' @return Integer number of occasions.
#' @keywords internal
getNumberOfOccasionsForModel = function( model ) {
  n = prop( model, "numberOfOccasions" )
  if ( length( n ) > 0L && !is.na( n[[1L]] ) && n[[1L]] >= 1L )
    return( as.integer( n[[1L]] ) )

  inferNumberOfOccasions(
    prop( model, "modelCovariates" ),
    prop( model, "modelParameters" )
  )
}

#' Whether the combination-by-occasion evaluation path is required.
#'
#' Returns \code{TRUE} when the model has covariates or more than one occasion
#' (including random IOV with \code{gamma > 0} and no occasion covariate).
#'
#' @param model A \code{Model} object.
#' @return Logical scalar.
#' @keywords internal
usesCovariateOccasionStructure = function( model ) {
  if ( length( prop( model, "modelCovariates" ) ) > 0L ) return( TRUE )
  getNumberOfOccasionsForModel( model ) > 1L
}

# Split a flat covariate list by short class name (strips the "PFIM::" prefix).
.splitCovariatesByClass = function( covariates ) {
  if ( length( covariates ) == 0L ) return( list() )
  split( covariates, map_chr( covariates, ~ str_remove( pluck( class( .x ), 1L ), "^PFIM::" ) ) )
}

# Return only covariates matching a given short class name.
.filterCovariatesByClass = function( covariates, shortClass ) {
  pluck( .splitCovariatesByClass( covariates ), shortClass, .default = list() )
}

# Builds a fully reversible combination name from a named list of
# (covariateName -> categoryOrSequenceLabel) pairs.
# Uses combinationSep ("::") and "=" so that underscores inside names are safe.
.encodeCombinationName = function( parts ) {
  if ( length( parts ) == 0L ) return( "Reference" )
  paste( names( parts ), unlist( parts ), sep = "=", collapse = combinationSep )
}

# Inverse of .encodeCombinationName.
.decodeCombinationName = function( combinationName ) {
  if ( combinationName == "Reference" ) return( list() )

  tokens = str_split( combinationName, fixed( combinationSep ) )[[1L]]
  pairs  = str_split( tokens, fixed( "=" ) )

  if ( any( map_int( pairs, length ) != 2L ) )
    stop( sprintf( "Malformed combination name: '%s'", combinationName ) )

  set_names(
    map( pairs, ~ .x[[2L]] ),
    map_chr( pairs, ~ .x[[1L]] )
  )
}

# Recovers the full covariate-value mapping for a combination name, injecting
# reference levels (first category / first sequence) for absent covariates.
parseCombinationName = function( combinationName, modelCovariates ) {
  decoded     = .decodeCombinationName( combinationName )
  covNames    = map_chr( modelCovariates, ~ prop( .x, "name" ) )
  missingMask = map_lgl( covNames, ~ is.null( decoded[[ .x ]] ) )

  refDefaults = map( modelCovariates[ missingMask ], function( cov ) {
    if ( S7::S7_inherits( cov, CategoricalCovariate ) ) {
      prop( cov, "categories" )[[1L]]
    } else {
      seqs = prop( cov, "sequences" )
      if ( is.null( names( seqs ) ) ) "sequence_1" else names( seqs )[[1L]]
    }
  }) |> set_names( covNames[ missingMask ] )

  c( decoded, refDefaults )
}
