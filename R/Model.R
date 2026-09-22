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
#'   Finite-difference gradients assume parameter mus are O(1); values with
#'   \eqn{|\mu| < 10^{-4}} use an absolute step floor.
#' @param parametersForComputingGradient Internal FD stencil (\code{eps^{1/3}}
#'   relative steps).
#' @param initialConditions ODE initial conditions.
#' @param functionArguments Names passed to the model function.
#' @param functionArgumentsSymbol Symbolic argument list for parsing.
#' @param numberOfOccasions Number of study occasions (set on the \code{Model} in
#'   \code{defineModelType()} from \code{Evaluation}/\code{Optimization}; validated against
#'   \code{\link{inferNumberOfOccasions}}).
#' @include CovariateModelEquation.R
#' @return An S7 object of class \code{Model}.
#' @export

Model = new_class( "Model", package = "PFIM",

                   properties = list(
                     name                           = new_property( class_character,        default = character(0) ),
                     modelParameters                = new_property( class_list,             default = list()       ),
                     modelCovariatesEquation        = new_property( NULL | CovariateModelEquation, default = NULL ),
                     modelCovariates                = new_property( class_list,             default = list()       ),
                     covariatesEffect               = new_property( class_list,             default = list()       ),
                     covariatesCombination          = new_property( class_list,             default = list()       ),
                     modelParametersWithCovariates  = new_property( class_list,             default = list()       ),
                     numberOfOccasions              = new_property( class_numeric,          default = 1            ),
                     omegaWithIOV                   = new_property( class_double,           default = numeric(0)   ),
                     samplings                      = new_property( class_numeric,          default = numeric(0)   ),
                     modelEquations                 = new_property( class_list,             default = list()       ),
                     wrapper                        = new_property( class_function | NULL,  default = NULL         ),
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
                   ),
                   validator = function( self ) {
                     checks = list(
                       list( prop( self, "modelParameters" ), ModelParameter, "Model:modelParameters", "ModelParameter" ),
                       list( prop( self, "modelCovariates" ), Covariate, "Model:modelCovariates", "Covariate" ),
                       list( prop( self, "modelError" ), ModelError, "Model:modelError", "ModelError" )
                     )
                     msgs = compact( map( checks, function( chk )
                       do.call( .validateS7List, chk ) ) )
                     if ( length( msgs ) )
                       return( msgs[[ 1L ]] )
                     NULL
                   } )

#' Compile model equations and output mapping before evaluation
#' @param model A \code{Model} object.
#' @param ... PFIMProject passed to methods.
#' @usage defineModelWrapper(model, ...)
#' @name defineModelWrapper
#' @keywords internal
defineModelWrapper                = new_generic( "defineModelWrapper",                c( "model" ) )

#' Bind dosing schedules and sampling grids to the model
#' @param model A \code{Model} object.
#' @param ... Arm passed to methods.
#' @usage defineModelAdministration(model, ...)
#' @name defineModelAdministration
#' @keywords internal
defineModelAdministration         = new_generic( "defineModelAdministration",         c( "model" ) )

#' Predict model responses at arm sampling times
#' @param model A \code{Model} object.
#' @param ... Arm passed to methods.
#' @usage evaluateModel(model, ...)
#' @name evaluateModel
#' @keywords internal
evaluateModel                     = new_generic( "evaluateModel",                     c( "model" ) )

#' Gradients of model responses with respect to parameters
#' @param model A \code{Model} object.
#' @param ... Arm passed to methods.
#' @usage evaluateModelGradient(model, ...)
#' @name evaluateModelGradient
#' @keywords internal
evaluateModelGradient             = new_generic( "evaluateModelGradient",             c( "model" ) )

#' Residual variance at sampling times for an arm
#' @param model A \code{Model} object.
#' @param ... Arm passed to methods.
#' @usage evaluateModelVariance(model, ...)
#' @name evaluateModelVariance
#' @keywords internal
evaluateModelVariance             = new_generic( "evaluateModelVariance",             c( "model" ) )

#' Evaluate ODE initial conditions from the arm
#' @param model A \code{Model} object.
#' @param ... Arm and optional bolus dose-event data (see methods).
#' @usage evaluateInitialConditions(model, ...)
#' @name evaluateInitialConditions
#' @keywords internal
evaluateInitialConditions         = new_generic( "evaluateInitialConditions",         c( "model" ) )

#' Finite-difference perturbations for gradient computation
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name finiteDifferenceHessian
#' @keywords internal
finiteDifferenceHessian           = new_generic( "finiteDifferenceHessian",           c( "model" ) )

#' Build covariate effect vectors for model parameters
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name evaluateCovariatesEffects
#' @keywords internal
evaluateCovariatesEffects         = new_generic( "evaluateCovariatesEffects",         c( "model" ) )

#' Construct omega matrix with IOV structure from covariates
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name evaluateOmegaMatrixFromCovariates
#' @keywords internal
evaluateOmegaMatrixFromCovariates = new_generic( "evaluateOmegaMatrixFromCovariates", c( "model" ) )

#' Generate covariate combination grid and proportions
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name generateCovariatesCombination
#' @keywords internal
generateCovariatesCombination     = new_generic( "generateCovariatesCombination",     c( "model" ) )

#' Remap library PK equations onto project compartments
#' @param pkModel Library PK model object.
#' @param pfimproject A \code{PFIMProject}, \code{Evaluation}, or \code{Optimization} object.
#' @param ... Optional method arguments.
#' @return List of PK model equations for the project compartments.
#' @examples
#' \dontrun{
#' vignette("LibraryOfModels")
#' }
#' @name definePKModel
#' @export
definePKModel                     = new_generic( "definePKModel",   c( "pkModel", "pfimproject" ) )

#' Combine PK and PD library equations for the project
#' @param pkModel Library PK model object.
#' @param pdModel Library PD model object.
#' @param pfimproject A \code{PFIMProject}, \code{Evaluation}, or \code{Optimization} object.
#' @param ... Optional method arguments.
#' @return Combined PK/PD equation list for the project.
#' @examples
#' \dontrun{
#' vignette("LibraryOfModels")
#' }
#' @name definePKPDModel
#' @export
definePKPDModel                   = new_generic( "definePKPDModel", c( "pkModel", "pdModel", "pfimproject" ) )

#' Build occasion-specific parameters for each covariate combination
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name modelParametersWithCovariates
#' @keywords internal
modelParametersWithCovariates     = new_generic( "modelParametersWithCovariates",     c( "model" ) )

#' Prepare all model covariate-derived data structures
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name defineCovariatesData
#' @keywords internal
defineCovariatesData              = new_generic( "defineCovariatesData",              c( "model" ) )

#' Report whether a model includes covariates
#' @param model A \code{Model} object.
#' @param ... Optional method arguments.
#' @name hasCovariates
#' @keywords internal
hasCovariates                     = new_generic( "hasCovariates",                     c( "model" ) )

#' Evaluate model outputs over covariate combinations and occasions
#' @param model A \code{Model} object.
#' @param ... Arm and evaluation core function (see methods).
#' @usage evaluateModelWithCovariates(model, ...)
#' @name evaluateModelWithCovariates
#' @keywords internal
evaluateModelWithCovariates       = new_generic( "evaluateModelWithCovariates",       c( "model" ) )


# Covariate-equation type codes used by gradient chain-rule branches.
modelCovEqExponential = 1L   # theta = mu * exp(beta * cov)
modelCovEqAdditive    = 2L   # theta = mu * (1 + beta * cov)

# Combination name encoding (cov=cat pairs joined by "::").
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
#'   \item Else if any parameter has \code{gamma > 0} (random IOV), default to 2
#'         occasions; set \code{numberOfOccasions >= 2} on the project for more periods.
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

#' Resolve \code{numberOfOccasions} for a project: infer when unset, else validate.
#' @param userN Value from \code{Evaluation}/\code{Optimization} (\code{NA} = infer).
#' @param modelCovariates List of covariate objects.
#' @param modelParameters List of \code{ModelParameter} objects.
#' @return Integer number of occasions.
#' @keywords internal
resolveNumberOfOccasions = function( userN, modelCovariates, modelParameters ) {
  expected = inferNumberOfOccasions( modelCovariates, modelParameters )
  if ( length( userN ) == 0L || ( length( userN ) == 1L && is.na( userN ) ) )
    return( expected )
  if ( length( userN ) != 1L )
    stop( "numberOfOccasions must be a single integer.", call. = FALSE )
  userN = as.integer( userN )
  if ( userN < 1L )
    stop( "numberOfOccasions must be >= 1.", call. = FALSE )

  has_iov_cov = length(
    .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )
  ) > 0L
  gamma_vals = map_dbl( modelParameters, function( p ) {
    g = prop( p, "gamma" )
    if ( length( g ) == 0L || is.na( g ) ) 0 else g
  })
  gamma_only_iov = !has_iov_cov && any( gamma_vals > 0, na.rm = TRUE )

  if ( gamma_only_iov ) {
    if ( userN < 2L )
      stop( "numberOfOccasions must be >= 2 when gamma (IOV) is set.", call. = FALSE )
    return( userN )
  }

  if ( userN != expected )
    stop(
      sprintf(
        "numberOfOccasions (%d) is inconsistent with covariates/IOV (expected %d).",
        userN, expected
      ),
      call. = FALSE
    )
  userN
}

#' Occasion count for a built \code{Model} object.
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

# Split covariates by short class name (strip "PFIM::" prefix).
#' Split covariates by unqualified class name.
#' @param covariates List of covariate objects.
#' @return Named list grouping covariates by short class.
#' @noRd
#' @keywords internal
.splitCovariatesByClass = function( covariates ) {
  if ( length( covariates ) == 0L ) return( list() )
  split( covariates, map_chr( covariates, ~ str_remove( pluck( class( .x ), 1L ), "^PFIM::" ) ) )
}

#' Filter covariates by a short class label.
#' @param covariates List of covariate objects.
#' @param shortClass Character scalar class key without package prefix.
#' @return List of covariates matching \code{shortClass}.
#' @noRd
#' @keywords internal
.filterCovariatesByClass = function( covariates, shortClass ) {
  pluck( .splitCovariatesByClass( covariates ), shortClass, .default = list() )
}

# covName=categoryOrSequence tokens, joined by combinationSep.
#' Encode covariate combination parts into one name.
#' @param parts Named list of \code{covariate=value} entries.
#' @return Character scalar combination name.
#' @noRd
#' @keywords internal
.encodeCombinationName = function( parts ) {
  if ( length( parts ) == 0L ) return( "Reference" )
  paste( names( parts ), unlist( parts ), sep = "=", collapse = combinationSep )
}

# Inverse of .encodeCombinationName.
#' Decode a combination name into covariate-value pairs.
#' @param combinationName Character scalar generated by \code{.encodeCombinationName()}.
#' @return Named list of decoded covariate values.
#' @noRd
#' @keywords internal
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

# fill reference categories for covariates omitted from the combination name.
#' Parse combination name and fill missing reference values.
#' @param combinationName Character scalar combination identifier.
#' @param modelCovariates List of covariate objects used to infer defaults.
#' @return Named list of covariate values including references.
#' @keywords internal
parseCombinationName = function( combinationName, modelCovariates ) {
  decoded     = .decodeCombinationName( combinationName )
  covNames    = map_chr( modelCovariates, ~ prop( .x, "name" ) )
  missingMask = map_lgl( covNames, ~ is.null( decoded[[ .x ]] ) )

  refDefaults = map( modelCovariates[ missingMask ], function( cov ) {
    if ( S7::S7_inherits( cov, CategoricalCovariate ) ) {
      prop( cov, "categories" )[[1L]]
    } else {
      seqs     = prop( cov, "sequences" )
      seqNames = names( seqs ) %||% paste0( "sequence_", seq_along( seqs ) )
      seqNames[[1L]]
    }
  }) |> set_names( covNames[ missingMask ] )

  c( decoded, refDefaults )
}
