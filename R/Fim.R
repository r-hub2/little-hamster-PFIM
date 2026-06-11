#' @title Fim
#' @description
#' Base class for Fisher information matrices (\code{PopulationFim},
#' \code{IndividualFim}, \code{BayesianFim}). All three support covariates and IOV;
#' \code{PopulationFim} additionally includes IIV/IOV variance blocks (\code{omega},
#' \code{gamma}).
#'
#' @section Choosing \code{fimType} and reading SEs:
#' \describe{
#'   \item{\code{population}}{Design FIM for the nonlinear mixed-effects model: estimable
#'     \eqn{\mu}, covariate \eqn{\beta}, IIV \eqn{\omega^2}, IOV \eqn{\gamma^2}, and \eqn{\sigma}.
#'     Use for classical PFIM design criteria and full \code{covariateTest} (all slots).}
#'   \item{\code{individual}}{Expected subject-level information with covariate/IOV averaging.
#'     Matrix \eqn{(\mu,\beta)} plus \eqn{\sigma}; no \eqn{\omega^2}/\eqn{\gamma^2} block. SEs on \eqn{\mu}
#'     are not shrunk; interpret \code{covariateTest} slots 2--3 for \eqn{\beta} only.}
#'   \item{\code{Bayesian}}{Shrunk design FIM on typical PK parameters with prior on \eqn{\mu};
#'     optional \eqn{\beta} block with covariates. Smaller SEs on \eqn{\mu} than individual; \eqn{\sigma}
#'     not in the displayed matrix. See \code{?covariateTest} for power on covariate effects.}
#' }
#' @param fisherMatrix Labelled n x n FIM matrix.
#' @param fixedEffects Fixed-effects sub-block.
#' @param varianceEffects Variance-effects sub-block.
#' @param SEAndRSE List of data frames: SE, RSE, and combined SE/RSE.
#' @param condNumberFixedEffects Condition number of the fixed-effects block.
#' @param condNumberVarianceEffects Condition number of the variance-effects block.
#' @param shrinkage Named numeric vector of Bayesian shrinkage (percent) per parameter.
#' @export

Fim = new_class( "Fim", package = "PFIM",
                 properties = list(
                   fisherMatrix              = new_property( class_double, default = numeric( 0 ) ),
                   fixedEffects              = new_property( class_double, default = numeric( 0 ) ),
                   varianceEffects           = new_property( class_double, default = numeric( 0 ) ),
                   SEAndRSE                  = new_property( class_list,   default = list()        ),
                   condNumberFixedEffects    = new_property( class_double, default = 0.0           ),
                   condNumberVarianceEffects = new_property( class_double, default = 0.0           ),
                   shrinkage                 = new_property( class_double, default = numeric( 0 )  )
                 )
)

# --- S7 generics ---
evaluateFim                = new_generic( "evaluateFim",                c( "fim", "model", "arm" ) )
evaluateVarianceFIM        = new_generic( "evaluateVarianceFIM",        c( "fim", "model", "arm" ) )
setEvaluationFim           = new_generic( "setEvaluationFim",           "fim"                      )
setOptimalArms             = new_generic( "setOptimalArms",             c( "fim", "optimizationAlgorithm" ) )
Dcriterion                 = new_generic( "Dcriterion",                 "fim"                      )
showFIM                    = new_generic( "showFIM",                    "fim"                      )
#' Barplot of standard errors from a FIM
#' @name plotSEFIM
#' @param fim        A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} object.
#' @param \dots      Additional arguments passed to methods.
#' @export
plotSEFIM                  = new_generic( "plotSEFIM",                  c( "fim", "evaluation" )   )

#' Barplot of relative standard errors from a FIM
#' @name plotRSEFIM
#' @param fim        A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} object.
#' @param \dots      Additional arguments passed to methods.
#' @export
plotRSEFIM                 = new_generic( "plotRSEFIM",                 c( "fim", "evaluation" )   )

#' plotShrinkage: plot Bayesian shrinkage bars.
#' @name plotShrinkage
#' @param fim        A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} object.
#' @param \dots      Additional arguments passed to methods.
#' @export
plotShrinkage              = new_generic( "plotShrinkage",              c( "fim", "evaluation" )   )

#' FIM tables for HTML reports
#' @name tablesForReport
#' @param fim        A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} object.
#' @param \dots      Additional arguments passed to methods.
#' @export
tablesForReport            = new_generic( "tablesForReport",            c( "fim", "evaluation" )   )

#' Render the evaluation HTML report
#' @name generateReportEvaluation
#' @param fim        A \code{Fim} object.
#' @param \dots      Additional arguments passed to methods.
#' @export
generateReportEvaluation   = new_generic( "generateReportEvaluation",   "fim"                      )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @param fim                   A \code{Fim} object.
#' @param optimizationAlgorithm A \code{Optimization} object.
#' @param \dots                 Additional arguments passed to methods.
#' @export
generateReportOptimization = new_generic( "generateReportOptimization", c( "fim", "optimizationAlgorithm" ) )

# Greek symbol tables
# Single source of truth used by all three FIM subclasses.
# Defined here (common ancestor) so @include Fim.R is sufficient for every file
# that needs them (no load-order dependency on PopulationFim.R).

#' Console / plain-text display (Unicode, separator already included)
.greekConsole = c(
  mu    = "\u03bc_",
  beta  = "\u03b2_",
  omega = "\u03c9\u00B2_",
  gamma = "\u03b3\u00B2_",
  sigma = "\u03c3"
)

#' ggplot axis labels (bare Unicode, no separator)
.greekPlot = c(
  mu    = "\u03bc",
  omega = "\u03c9\u00B2",
  gamma = "\u03b3\u00B2",
  sigma = "\u03c3"
)

#' LaTeX for kableExtra HTML tables
.greekLatex = c(
  mu    = "$\\mu_{",
  beta  = "$\\beta_{",
  omega = "$\\omega^2_{",
  gamma = "$\\gamma^2_{",
  sigma = "$\\sigma^2_{"
)

# Shared computational helpers
#' Build a (nTimes x nParams) gradient matrix from an arm's gradients.
#' Aggregates combination-by-occasion gradients when \code{model} is supplied and
#' the design uses covariates or IOV.
#'
#' @param arm   An \code{Arm} object.
#' @param cols  Character vector of column names to select.
#' @param model Optional \code{Model} object (required for nested gradient layout).
#' @return Numeric matrix (nTimes x length(cols)).
#' @keywords internal
.gradientMatrix = function( arm, cols, model = NULL ) {
  mat = if ( !is.null( model ) && usesCovariateOccasionStructure( model ) ) {
    getArmEvaluationGradientsMatrix( model, arm )
  } else {
    raw = prop( arm, "evaluationGradients" )
    df  = if ( is.data.frame( raw ) ) raw else do.call( rbind, raw )
    as.matrix( df )
  }
  mat[ , cols, drop = FALSE ]
}

#' Compute the symmetric (nLambda x nLambda) variance-FIM block.
#'
#' Entry (i, j) is the Cramer-Rao contribution of the (i, j) pair of variance
#' parameters:
#' \deqn{M_{\lambda}[i,j] = \frac{1}{2}\,\mathrm{Tr}\!\left(
#'   V^{-1}\frac{\partial V}{\partial\lambda_i}
#'   V^{-1}\frac{\partial V}{\partial\lambda_j}\right)
#'   = \frac{1}{2}\langle T_i,\,\frac{\partial V}{\partial\lambda_j}\rangle_F}
#' where \eqn{T_k = V^{-1}(\partial V/\partial\lambda_k)V^{-1}} and
#' \eqn{\langle A, B\rangle_F = \sum_{ij} A_{ij}B_{ij}} for symmetric matrices.
#'
#' @param V_inv    Dense numeric matrix \eqn{V^{-1}} computed via Cholesky.
#' @param dV_list  Named list of \eqn{\partial V/\partial\lambda_k} matrices.
#' @return Symmetric (n x n) numeric matrix.
#' @keywords internal
.computeMFVar = function( V_inv, dV_list ) {
  T_mats = map( dV_list, ~ V_inv %*% .x %*% V_inv )
  n      = length( T_mats )
  outer(
    seq_len( n ), seq_len( n ),
    Vectorize( function( i, j ) 0.5 * sum( T_mats[[ i ]] * dV_list[[ j ]] ) )
  )
}

#' @keywords internal
.setOptimalArmsMultiplicative = function( optimizationAlgorithm ) {
  out          = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
  weights      = out$multiplicativeAlgorithmOutput[[ "weights" ]]
  weightsIndex = which( weights > out$weightThreshold )
  armList = map( weightsIndex, function( idx ) {
    arm = pluck( out$armFims[[ idx ]], 1L )
    prop( arm, "size" ) = 1.0
    prop( arm, "name" ) = paste0( "Arm", idx )
    arm
  })
  armList[ rev( order( map_dbl( armList, ~ prop( .x, "size" ) ) ) ) ]
}

#' @keywords internal
.setOptimalArmsFedorovWynn = function( optimizationAlgorithm ) {
  out = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
  imap( out$listArms, ~ {
    prop( .x$arm, "name" ) = paste0( "Arm", .y )
    prop( .x$arm, "size" ) = 1.0
    .x
  })
}

#' Absolute path to an Rmd skeleton template in the PFIM package.
#' Called lazily at render time, never at package load time.
#' @param filename Basename of the Rmd template file.
#' @keywords internal
.reportTemplatePath = function( filename ) {
  dir = file.path( system.file( package = "PFIM" ),
                    "rmarkdown", "templates", "skeleton" )
  path = file.path( dir, filename )
  if ( file.exists( path ) ) return( path )
  stem = sub( "\\.[rR]md$", "", filename )
  hits = list.files( dir,
                      pattern = paste0( "^", stem, "\\.[rR]md$" ),
                      full.names = TRUE, ignore.case = TRUE )
  if ( length( hits ) >= 1L ) return( hits[[ 1L ]] )
  path
}

#' create a \code{generateReportOptimization} method body.
#'
#' Returns a closure that calls \code{rmarkdown::render} with the given
#' template.  The template string is captured at factory-call time via
#' \code{force()}, ensuring each method gets its own distinct value.
#'
#' Use:
#' \code{method(generateReportOptimization, list(MyFim, MyAlgo)) =
#'   .renderReport("MyTemplate.rmd")}
#'
#' @param template Basename of the optimisation Rmd template.
#' @return A function with the standard \code{generateReportOptimization} signature.
#' @keywords internal
.renderReport = function( template ) {
  force( template )
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath )
    rmarkdown::render(
      input       = .reportTemplatePath( template ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
}

#' create a \code{generateReportEvaluation} method body.
#' @param template Basename of the evaluation Rmd template.
#' @return A function with the standard \code{generateReportEvaluation} signature.
#' @keywords internal
.renderEvalReport = function( template ) {
  force( template )
  function( fim, tablesForReport, outputFile, outputPath )
    rmarkdown::render(
      input       = .reportTemplatePath( template ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
}

#' D-criterion from a Fisher information matrix
#' @name Dcriterion
#' @export

method( Dcriterion, Fim ) = function( fim ) {
  M      = prop( fim, "fisherMatrix" )
  p      = nrow( M )
  logDet = as.numeric( determinant( M, logarithm = TRUE )$modulus )
  if ( is.na( logDet ) || is.infinite( logDet ) )
    return( as.numeric( det( M )^( 1 / p ) ) )
  as.numeric( exp( logDet / p ) )
}

.conditionNumber = function( M ) {
  if ( nrow( M ) == 0L ) return( Inf )
  ev = eigen( M, symmetric = TRUE, only.values = TRUE )$values
  max( abs( ev ) ) / min( abs( ev ) )
}

#' Console table of parameter values, SE and RSE (cat-based for Windows / RStudio).
#' @param fim A \code{Fim} object.
#' @param evaluation Optional \code{PFIMProject}; when supplied, \code{setEvaluationFim}
#'   is called before printing.
#' @keywords internal

.printFimSeAndRse = function( fim, evaluation = NULL ) {
  if ( !is.null( evaluation ) )
    fim = setEvaluationFim( fim, evaluation )

  se_list = as.list( prop( fim, "SEAndRSE" ) )
  se_tbl  = se_list[[ "table" ]] %||% se_list[[ "SEAndRSE" ]]

  if ( !is.data.frame( se_tbl ) || nrow( se_tbl ) == 0L ) {
    cat( "  (SE/RSE table unavailable)\n\n" )
    flush.console()
    return( invisible( NULL ) )
  }

  rn = rownames( se_tbl )
  if ( is.null( rn ) ) rn = paste0( "p", seq_len( nrow( se_tbl ) ) )

  fmt = function( x, w ) format( x, width = w, justify = "right", trim = TRUE )
  wn  = max( 12L, max( nchar( rn, type = "width" ), na.rm = TRUE ) )

  cat(
    sprintf( paste0( "%-", wn, "s %14s %12s %10s\n" ),
             "Parameter", "Value", "SE", "RSE(%)" )
  )
  for ( i in seq_len( nrow( se_tbl ) ) ) {
    cat(
      sprintf(
        paste0( "%-", wn, "s %14s %12s %10s\n" ),
        rn[i],
        fmt( se_tbl$parametersValues[i], 14 ),
        fmt( se_tbl$SE[i], 12 ),
        fmt( se_tbl$RSE[i], 10 )
      )
    )
  }
  cat( "\n" )
  flush.console()
  invisible( se_tbl )
}

#' @keywords internal
.duplicateFim = function( fim ) {
  cls = S7::S7_class( fim )@name
  if ( !cls %in% c( "PopulationFim", "IndividualFim", "BayesianFim" ) )
    stop(
      sprintf( "Cannot duplicate FIM object of class '%s'.", cls ),
      call. = FALSE
    )
  clone = unserialize( serialize( fim, NULL ) )
  S7::set_props(
    clone,
    fisherMatrix              = numeric( 0 ),
    fixedEffects              = numeric( 0 ),
    varianceEffects           = numeric( 0 ),
    SEAndRSE                  = list(),
    condNumberFixedEffects    = 0.0,
    condNumberVarianceEffects = 0.0,
    shrinkage                 = numeric( 0 )
  )
}
