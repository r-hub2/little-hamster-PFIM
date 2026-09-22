# FIM class, generics, numeric helpers (D-criterion, condition number).
# Labels/plots: pfim-fim-labels.R; discrete arms: pfim-fim-optimal-arms.R;
# HTML reports: pfim-fim-report-render.R.

#' Square-matrix guard for FIM properties (empty default is allowed).
#' @noRd
#' @keywords internal
.fimValidateSquareMatrix = function( value ) {
  if ( length( value ) == 0L ) return( NULL )
  if ( !is.matrix( value ) ) return( "must be a matrix" )
  if ( nrow( value ) != ncol( value ) ) return( "must be square" )
  if ( !isSymmetric( unname( value ), tol = 1e-8 ) ) return( "must be symmetric" )
  NULL
}

#' @title Fim
#' @description
#' Base class for Fisher information matrices (\code{PopulationFim},
#' \code{IndividualFim}, \code{BayesianFim}).
#' @param fisherMatrix Labelled FIM matrix.
#' @param fixedEffects Fixed-effects sub-block.
#' @param varianceEffects Variance-effects sub-block.
#' @param SEAndRSE List with \code{SE}, \code{RSE}, and combined tables.
#' @param condNumberFixedEffects Condition number of the fixed-effects block.
#' @param condNumberVarianceEffects Condition number of the variance-effects block.
#' @param shrinkage Named shrinkage (percent) per parameter (Bayesian FIM).
#' @param singularFim \code{TRUE} when SE/RSE used a pseudo-inverse (singular FIM).
#'   Set in \code{.fimStoreEvaluationResult} from \code{.fimBuildSeAndRse}; D-criterion
#'   / log-det still use the raw matrix (may be 0 / \code{-Inf} when singular).
#' @return An S7 object holding FIM matrices, SE/RSE tables, and condition numbers.
#' @export
Fim = new_class( "Fim", package = "PFIM",
                 properties = list(
                   fisherMatrix              = new_property(
                     class_double,
                     default   = numeric( 0 ),
                     validator = .fimValidateSquareMatrix
                   ),
                   fixedEffects              = new_property( class_double,  default = numeric( 0 ) ),
                   varianceEffects           = new_property( class_double,  default = numeric( 0 ) ),
                   SEAndRSE                  = new_property( class_list,    default = list()        ),
                   condNumberFixedEffects    = new_property( class_double,  default = 0.0           ),
                   condNumberVarianceEffects = new_property( class_double,  default = 0.0           ),
                   shrinkage                 = new_property( class_double,  default = numeric( 0 )  ),
                   singularFim               = new_property( class_logical, default = FALSE         )
                 )
)

# S7 generics
#' @name evaluateFim
#' @return The updated \code{Fim} object after evaluation.
#' @keywords internal
evaluateFim                = new_generic( "evaluateFim",                c( "fim", "model", "arm" ) )

#' @name evaluateVarianceFIM
#' @return A variance contribution matrix for the FIM.
#' @keywords internal
evaluateVarianceFIM        = new_generic( "evaluateVarianceFIM",        c( "fim", "model", "arm" ) )

#' Attach evaluated FIM results to a project
#' @param fim A \code{Fim} object.
#' @param ... Method arguments (see methods).
#' @usage setEvaluationFim(fim, ...)
#' @name setEvaluationFim
#' @return The modified \code{Fim} object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' PFIM:::setEvaluationFim(prop(ev, "fim"), ev)
#' }
#' @keywords internal
setEvaluationFim           = new_generic( "setEvaluationFim",           "fim"                      )

#' Optimal arms from design optimization
#' @param fim A \code{Fim} object.
#' @param optimizationAlgorithm Optimizer object from \code{run(Optimization)}.
#' @param ... Not used by current methods.
#' @usage setOptimalArms(fim, optimizationAlgorithm, ...)
#' @name setOptimalArms
#' @return The modified \code{Fim} object.
#' @keywords internal
setOptimalArms             = new_generic( "setOptimalArms",             c( "fim", "optimizationAlgorithm" ) )
Dcriterion                 = new_generic( "Dcriterion",                 "fim"                      )

#' Print FIM summaries to the console
#' @param fim A \code{Fim} object.
#' @param ... Not used by current methods.
#' @usage showFIM(fim, ...)
#' @name showFIM
#' @return Invisibly, the input \code{Fim} object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' showFIM(prop(ev, "fim"))
#' }
#' @export
showFIM                    = new_generic( "showFIM",                    "fim"                      )

#' Plot standard errors from a FIM
#' @param fim A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} evaluation object.
#' @param ... Graphics parameters passed to \code{ggplot2}.
#' @usage plotSEFIM(fim, evaluation, ...)
#' @name plotSEFIM
#' @return A \code{ggplot2} plot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' PFIM:::plotSEFIM(prop(ev, "fim"), ev)
#' }
#' @keywords internal
plotSEFIM                  = new_generic( "plotSEFIM",                  c( "fim", "evaluation" )   )

#' Plot relative standard errors from a FIM
#' @param fim A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} evaluation object.
#' @param ... Graphics parameters passed to \code{ggplot2}.
#' @usage plotRSEFIM(fim, evaluation, ...)
#' @name plotRSEFIM
#' @return A \code{ggplot2} plot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' PFIM:::plotRSEFIM(prop(ev, "fim"), ev)
#' }
#' @keywords internal
plotRSEFIM                 = new_generic( "plotRSEFIM",                 c( "fim", "evaluation" )   )

#' Plot Bayesian shrinkage from a FIM
#' @param fim A \code{BayesianFim} object.
#' @param evaluation A \code{PFIMProject} evaluation object.
#' @param ... Graphics parameters passed to \code{ggplot2}.
#' @usage plotShrinkage(fim, evaluation, ...)
#' @name plotShrinkage
#' @export
plotShrinkage              = new_generic( "plotShrinkage",              c( "fim", "evaluation" )   )

#' FIM tables for HTML reports
#' @param fim A \code{Fim} object.
#' @param evaluation A \code{PFIMProject} evaluation object.
#' @param ... Not used by current methods.
#' @usage tablesForReport(fim, evaluation, ...)
#' @name tablesForReport
#' @return Report tables as a list or data frame.
#' @keywords internal
tablesForReport            = new_generic( "tablesForReport",            c( "fim", "evaluation" )   )

#' Generate an evaluation report from FIM results
#' @param fim A \code{Fim} object.
#' @param ... Report paths and options (see methods).
#' @usage generateReportEvaluation(fim, ...)
#' @name generateReportEvaluation
#' @return Invisibly, the generated report path or render result.
#' @keywords internal
generateReportEvaluation   = new_generic( "generateReportEvaluation",   "fim"                      )

#' Generate an optimization report from FIM results
#' @param fim A \code{Fim} object.
#' @param optimizationAlgorithm Optimizer object with optimal design outputs.
#' @param ... Report paths and options (see methods).
#' @usage generateReportOptimization(fim, optimizationAlgorithm, ...)
#' @name generateReportOptimization
#' @return Invisibly, the generated report path or render result.
#' @keywords internal
generateReportOptimization = new_generic( "generateReportOptimization", c( "fim", "optimizationAlgorithm" ) )

#' Deep copy an S7 object (serialize-based; isolates mutable model state per arm).
#' @noRd
#' @keywords internal
.pfimCloneS7 = function( object, reset = list() ) {
  # Deep-copy via native serialize (xdr=FALSE is faster on the same machine).
  clone = unserialize( serialize( object, NULL, xdr = FALSE ) )
  nms = names( reset )
  if ( !length( nms ) ) return( clone )
  reduce( nms, function( obj, nm ) {
    prop( obj, nm ) = reset[[ nm ]]
    obj
  }, .init = clone )
}

#' Coerce a FIM payload to a square matrix; empty if not reconstructible.
#' @noRd
#' @keywords internal
.fimCoerceSquareMatrix = function( M ) {
  if ( is.null( M ) || length( M ) == 0L )
    return( matrix( numeric( 0 ), 0L, 0L ) )
  if ( is.matrix( M ) ) {
    if ( nrow( M ) != ncol( M ) )
      return( matrix( numeric( 0 ), 0L, 0L ) )
    return( M )
  }
  n = length( M )
  p = as.integer( round( sqrt( n ) ) )
  if ( p * p != n )
    return( matrix( numeric( 0 ), 0L, 0L ) )
  matrix( M, nrow = p, ncol = p )
}

#' Log-determinant of a symmetric FIM (0 or negative sign -> -Inf).
#' @noRd
#' @keywords internal
.fimLogDeterminant = function( M ) {
  M = .fimCoerceSquareMatrix( M )
  p = nrow( M )
  if ( is.null( p ) || p == 0L ) return( -Inf )
  ld = determinant( M, logarithm = TRUE )
  # Negative sign means det(M) <= 0 (singular or indefinite for our purposes).
  if ( ld$sign < 0 ) return( -Inf )
  as.numeric( ld$modulus )
}

#' Determinant of a Fisher matrix (0 when indefinite/singular; may be Inf if overflow).
#'
#' Prefer \code{.fimLogDeterminant()} for display when \eqn{|\log\det|} is large.
#' @noRd
#' @keywords internal
.fimDeterminant = function( M ) {
  ld = .fimLogDeterminant( M )
  # Map non-positive / singular FIMs to 0 rather than NaN from exp(-Inf).
  if ( !is.finite( ld ) ) return( 0 )
  if ( ld > log( .Machine$double.xmax ) ) return( Inf )
  if ( ld < log( .Machine$double.xmin ) ) return( 0 )
  exp( ld )
}

#' D-criterion from a Fisher matrix (\eqn{\det(M)^{1/p}}).
#'
#' Geometric mean of eigenvalues; returns 0 when \eqn{\log\det} is non-finite
#' (singular / indefinite). Prefer this over raw \code{det()} for optimizer fitness.
#' @noRd
#' @keywords internal
.fimDcriterionFromMatrix = function( M ) {
  if ( !is.matrix( M ) || nrow( M ) == 0L ) return( 0 )
  p = nrow( M )
  ld = .fimLogDeterminant( M )
  if ( !is.finite( ld ) ) return( 0 )
  # Geometric mean of eigenvalues: det(M)^{1/p} = exp(logdet / p).
  as.numeric( exp( ld / p ) )
}

#' D-criterion
#' @name Dcriterion
#' @return Numeric D-optimality criterion value.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' PFIM:::Dcriterion(prop(ev, "fim"))
#' }
#' @keywords internal

method( Dcriterion, Fim ) = function( fim ) {
  # D-opt = det(F)^{1/p} from the labelled fisherMatrix property.
  .fimDcriterionFromMatrix( prop( fim, "fisherMatrix" ) )
}

#' @noRd
#' @keywords internal
.conditionNumber = function( M ) {
  M = .fimCoerceSquareMatrix( M )
  if ( nrow( M ) == 0L ) return( Inf )
  ev = eigen( M, symmetric = TRUE, only.values = TRUE )$values
  scale = max( abs( ev ), 0 )
  tol   = max( scale * .Machine$double.eps * max( length( ev ), 1L ), 0 )
  # Tiny negative eigenvalues from roundoff are ignored; clear negatives => Inf.
  if ( any( ev < -tol ) ) return( Inf )
  # Near-zero eigenvalues => Inf kappa (do not drop them - that would hide singularity).
  if ( any( ev <= tol ) ) return( Inf )
  max( ev ) / min( ev )
}

#' Print a console note when SE/RSE used a pseudo-inverse (singular FIM).
#' @noRd
#' @keywords internal
.printFimSingularNote = function( fim ) {
  if ( !isTRUE( prop( fim, "singularFim" ) ) ) return( invisible( NULL ) )
  cat(
    "\n*** singularFim = TRUE ***\n",
    "Fisher information matrix was singular or not positive definite;\n",
    "SE/RSE used a Moore-Penrose pseudo-inverse (see warnings).\n",
    sep = ""
  )
  invisible( NULL )
}

#' Always print the \code{singularFim} flag (long note stays with SE/RSE).
#' @noRd
#' @keywords internal
.fimPrintSingularStatus = function( fim ) {
  cat( "singularFim:", isTRUE( prop( fim, "singularFim" ) ), "\n" )
  invisible( NULL )
}

#' Store labelled FIM blocks and SE/RSE on any \code{Fim} subclass.
#'
#' Shared by Population / Individual / Bayesian \code{setEvaluationFim} so the
#' on-object shape of \code{SEAndRSE} and condition numbers stays identical.
#' @param fim \code{Fim} object.
#' @param M Labelled Fisher matrix.
#' @param fixedEffects Fixed-effect block.
#' @param se List from \code{.fimBuildSeAndRse}.
#' @param varianceEffects Optional variance block (\code{NULL} for Bayesian).
#' @param shrinkage Optional Bayesian shrinkage matrix.
#' @return Updated \code{fim}.
#' @noRd
#' @keywords internal
.fimStoreEvaluationResult = function( fim, M, fixedEffects, se,
                                      varianceEffects = NULL,
                                      shrinkage = NULL ) {
  prop( fim, "fisherMatrix"           ) = M
  prop( fim, "fixedEffects"           ) = fixedEffects
  prop( fim, "condNumberFixedEffects" ) = .conditionNumber( fixedEffects )
  if ( !is.null( varianceEffects ) ) {
    prop( fim, "varianceEffects"           ) = varianceEffects
    prop( fim, "condNumberVarianceEffects" ) = .conditionNumber( varianceEffects )
  }
  if ( !is.null( shrinkage ) )
    prop( fim, "shrinkage" ) = shrinkage
  # Always the full list from .fimBuildSeAndRse (SE, RSE, table, SEAndRSE, singular).
  # singularFim tracks SE path only; D-criterion / kappa are computed separately above.
  prop( fim, "SEAndRSE"    ) = se
  prop( fim, "singularFim" ) = isTRUE( se$singular )
  fim
}

#' HTML note for singular FIM in evaluation / optimization reports.
#' @return Character HTML (empty when not singular).
#' @noRd
#' @keywords internal
.fimSingularHtmlNote = function( fim ) {
  if ( !isTRUE( prop( fim, "singularFim" ) ) ) return( "" )
  paste0(
    '<p><strong>singularFim = TRUE</strong>: Fisher information matrix was ',
    'singular or not positive definite; SE/RSE used a Moore-Penrose ',
    'pseudo-inverse. Interpret SE/RSE and design criteria with caution.</p>\n'
  )
}

#' Standard Pop/Ind report tables (criteria + FE/VE + SE/RSE kables).
#' @noRd
#' @keywords internal
.fimTablesForReportStandard = function( fim, evaluation ) {
  fim      = setEvaluationFim( fim, evaluation )
  SEAndRSE = prop( fim, "SEAndRSE" )$SEAndRSE %||% prop( fim, "SEAndRSE" )$table
  M        = prop( fim, "fisherMatrix" )
  fe       = as.matrix( prop( fim, "fixedEffects" ) )
  ve       = as.matrix( prop( fim, "varianceEffects" ) )
  cn1      = prop( fim, "condNumberFixedEffects" )
  cn2      = prop( fim, "condNumberVarianceEffects" )

  feLabels = .fimFixedEffectLatexLabels( evaluation, fixedEffects = fe )
  veLabels = .fimConsoleToLatexLabels( rownames( ve ) )
  dimnames( fe ) = list( feLabels, feLabels )
  dimnames( ve ) = list( veLabels, veLabels )

  paramLabels = .fimSeRseReportLabels(
    evaluation, SEAndRSE, fixedEffects = fe, varianceEffects = ve
  )
  list(
    fixedEffectsTable    = .kblReportStyled( fe ),
    varianceEffectsTable = .kblReportStyled( ve ),
    FIMCriteriaTable     = .fimCriteriaKable(
      .fimLogDeterminant( M ), Dcriterion( fim ), cn1, cn2,
      singularFim = isTRUE( prop( fim, "singularFim" ) )
    ),
    SEAndRSETable        = .fimSeRseKable( paramLabels, SEAndRSE ),
    singularFimNote      = .fimSingularHtmlNote( fim ),
    singularFim          = isTRUE( prop( fim, "singularFim" ) )
  )
}

#' Console legend for beta / gamma symbols when present in FIM rownames.
#' @noRd
#' @keywords internal
.fimShowSymbolLegend = function( rn, include_omega = FALSE, include_sigma = FALSE,
                                 mu_label = "\u03bc  = fixed effects (population means)",
                                 sigma_label = "\u03c3  = residual error SD" ) {
  has_beta  = any( grepl( "\u03b2_", rn, fixed = TRUE ) )
  has_gamma = any( grepl( "\u03b3\u00B2", rn, fixed = TRUE ) )
  # Only when covariate / IOV symbols appear (same rule as former showFIM methods).
  if ( !has_beta && !has_gamma )
    return( invisible( NULL ) )
  cat( "\n*************************************** \n Legend: \n" )
  cat( " ", mu_label, "\n", sep = "" )
  if ( has_beta )
    cat( " \u03b2  = covariate effects\n" )
  if ( include_omega )
    cat( " \u03c9\u00B2 = inter-individual variability (IIV)\n" )
  if ( has_gamma )
    cat( " \u03b3\u00B2 = inter-occasion variability (IOV)\n" )
  if ( include_sigma )
    cat( " ", sigma_label, "\n", sep = "" )
  cat( "*************************************** \n\n" )
  invisible( NULL )
}

#' @noRd
#' @keywords internal
.printFimSeAndRse = function( fim, evaluation = NULL ) {
  # Optionally refresh SE/RSE from a project evaluation before printing.
  if ( !is.null( evaluation ) )
    fim = setEvaluationFim( fim, evaluation )

  .printFimSingularNote( fim )

  se_list = as.list( prop( fim, "SEAndRSE" ) )
  # Older objects used key "SEAndRSE"; current code stores "table".
  se_tbl  = se_list[[ "table" ]] %||% se_list[[ "SEAndRSE" ]]

  if ( !is.data.frame( se_tbl ) || nrow( se_tbl ) == 0L ) {
    cat( "  (SE/RSE table unavailable)\n\n" )
    flush.console()
    return( invisible( NULL ) )
  }

  rn = rownames( se_tbl )
  if ( is.null( rn ) ) rn = paste0( "p", seq_len( nrow( se_tbl ) ) )

  fmt = function( x, w ) format( x, width = w, justify = "right", trim = TRUE )
  # Column width for names: at least 12, or widest parameter label.
  wn  = max( 12L, max( nchar( rn, type = "width" ), na.rm = TRUE ) )

  cat(
    sprintf( paste0( "%-", wn, "s %14s %12s %10s\n" ),
             "Parameter", "Value", "SE", "RSE(%)" )
  )
  row_txt = sprintf(
    paste0( "%-", wn, "s %14s %12s %10s\n" ),
    rn,
    fmt( se_tbl$parametersValues, 14 ),
    fmt( se_tbl$SE, 12 ),
    fmt( se_tbl$RSE, 10 )
  )
  cat( paste( row_txt, collapse = "" ), "\n", sep = "" )
  flush.console()
  invisible( se_tbl )
}

#' @noRd
#' @keywords internal
.duplicateFim = function( fim ) {
  cls = S7::S7_class( fim )@name
  # Subclass-specific copy factories are registered via pfim_register_fim_type.
  if ( !exists( cls, envir = .pfimFimDuplicatorRegistry, inherits = FALSE ) )
    stop(
      sprintf(
        "Cannot duplicate FIM object of class '%s' (register with pfim_register_fim_type).",
        cls
      ),
      call. = FALSE
    )
  get( cls, envir = .pfimFimDuplicatorRegistry )( fim )
}
