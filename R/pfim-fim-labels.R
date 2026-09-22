# FIM display labels (console / plotmath) and ggplot helpers.

#' @include Fim.R
#' @noRd
#' @keywords internal
NULL

# Greek labels (mu, beta, omega, gamma, sigma) for all FIM subclasses.

.greekConsole = c(
  mu    = "\u03bc_",
  beta  = "\u03b2_",
  omega = "\u03c9\u00B2_",
  gamma = "\u03b3\u00B2_",
  sigma = "\u03c3_"
)

.greekPlotmath = c(
  mu    = "mu",
  beta  = "beta",
  omega = "omega^2",
  gamma = "gamma^2",
  sigma = "sigma"
)

#' Plotmath facet labels for SE/RSE bar charts (device-independent Greek).
#' @noRd
#' @keywords internal
.pfimSeRseFacetLabel = function( metric, key ) {
  g = .greekPlotmath[[ key ]]
  if ( !.pfimIsNonEmptyScalar( g ) )
    stop( "Unknown Greek key for facet label: ", key, call. = FALSE )
  # plotmath: literal metric text, then Greek symbol (parsed by label_parsed).
  paste0( "'", metric, "'~'  '~", g )
}

#' Plotmath subscript for a parameter suffix (handles underscores).
#' @noRd
#' @keywords internal
.pfimPlotmathSubscript = function( suffix ) {
  # Underscores break bare plotmath tokens; quote the whole subscript.
  if ( grepl( "_", suffix, fixed = TRUE ) )
    paste0( "['", suffix, "']" )
  else
    paste0( "[", suffix, "]" )
}

#' Normalize internal/console parameter names to ASCII tokens for plotmath.
#' @noRd
#' @keywords internal
.pfimParamNameAscii = function( parameterName ) {
  # Map console Greek prefixes to ASCII tokens that plotmath can parse.
  name = as.character( parameterName )
  name = gsub( "\u03bc_",    "mu_",      name, fixed = TRUE )
  name = gsub( "\u03b2_",    "beta_",    name, fixed = TRUE )
  name = gsub( "\u03c9\u00b2_", "omega^2_", name, fixed = TRUE )
  name = gsub( "\u03b3\u00b2_", "gamma^2_", name, fixed = TRUE )
  name = gsub( "\u03c3\u00b2_", "sigma_", name, fixed = TRUE )
  name = gsub( "\u03c3_",       "sigma_", name, fixed = TRUE )
  name
}

#' Plotmath label for a FIM parameter name (mu, beta, omega^2, ...).
#' @noRd
#' @keywords internal
.pfimParamPlotmathLabel = function( parameterName ) {
  name = .pfimParamNameAscii( parameterName )
  # Longer prefixes first so "omega^2_" is not matched as a bare "omega".
  prefixes = c( "omega^2_", "gamma^2_", "sigma^2_", "sigma_", "mu_", "beta_" )
  greek    = c( "omega^2",  "gamma^2",  "sigma^2",  "sigma",  "mu",  "beta" )
  i = detect_index( prefixes, ~ startsWith( name, .x ) )
  if ( i == 0L ) return( name )
  suffix = substring( name, nchar( prefixes[[ i ]] ) + 1L )
  paste0( greek[[ i ]], .pfimPlotmathSubscript( suffix ) )
}

#' Plotmath y-axis label for sensitivity plots: df/d parameter.
#' @noRd
#' @keywords internal
.pfimSensitivityYLab = function( parameterName ) {
  # plotmath fraction df / d(parameter) for sensitivity figure y-axis.
  paste0( "frac(df, d*", .pfimParamPlotmathLabel( parameterName ), ")" )
}

#' Parsed plotmath expression from a label string.
#' @noRd
#' @keywords internal
.pfimParsePlotmath = function( label ) {
  # First (only) expression from parse(); keep.source=FALSE for speed.
  parse( text = label, keep.source = FALSE )[[ 1L ]]
}

#' Caption / x-axis annotation for sensitivity plots with Greek parameter names.
#' @noRd
#' @keywords internal
.pfimPlotmathEscape = function( x )
  gsub( "'", "\\\\'", x, fixed = TRUE )

#' @noRd
#' @keywords internal
.pfimSensitivityCaption = function( designName, armName, outputName, parameterName ) {
  paste0(
    "paste('Design: ", .pfimPlotmathEscape( gsub( "_", " ", designName, fixed = TRUE ) ),
    "   Arm: ", .pfimPlotmathEscape( armName ),
    "   Output: ", .pfimPlotmathEscape( outputName ),
    "   Parameter: ', ", .pfimParamPlotmathLabel( parameterName ), ")"
  )
}

#' Full x-axis plotmath label: time unit atop sensitivity caption.
#' @noRd
#' @keywords internal
.pfimSensitivityXLab = function( unitXAxis, designName, armName, outputName, parameterName ) {
  paste0(
    "atop('Time (", .pfimPlotmathEscape( unitXAxis ), ")', ",
    .pfimSensitivityCaption( designName, armName, outputName, parameterName ),
    ")"
  )
}

#' Require complete SE/RSE rownames (shared by pop / Bayesian plot + report paths).
#' @noRd
#' @keywords internal
.pfimSeRownames = function( seDF, context ) {
  rn = rownames( seDF )
  if ( !length( rn ) || length( rn ) != nrow( seDF ) )
    stop( sprintf( "%s: SEAndRSE rownames missing or incomplete.", context ), call. = FALSE )
  rn
}

.greekLatex = c(
  mu    = "$\\mu_{",
  beta  = "$\\beta_{",
  omega = "$\\omega^2_{",
  gamma = "$\\gamma^2_{",
  sigma = "$\\sigma_{"
)

#' Strip a Greek name prefix from FIM column labels.
#' @noRd
#' @keywords internal
.stripGreekPrefix = function( names, prefix )
  sub( paste0( "^", prefix ), "", names )

#' Residual-error parameter labels (\code{inter_}/\code{slope_} blocks).
#' @noRd
#' @keywords internal
.sigmaNames = function( modelError, greekSigma ) {
  # Filter must match .pfimSigmaIsEstimable / residual_error_derivatives.
  modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( .pfimSigmaIsEstimable( prop( err, "sigmaInter" ), prop( err, "sigmaInterFixed" ) ) )
          paste0( greekSigma, "inter_", out ),
        if ( .pfimSigmaIsEstimable( prop( err, "sigmaSlope" ), prop( err, "sigmaSlopeFixed" ) ) )
          paste0( greekSigma, "slope_", out )
      )
    } ) |> unlist( use.names = FALSE )
}

#' Estimable residual-error values in \code{inter}/\code{slope} column order.
#' @noRd
#' @keywords internal
.sigmaValues = function( modelError ) {
  # Same order and filter as .sigmaNames so values align with column labels.
  modelError |>
    map( function( err ) {
      c(
        if ( .pfimSigmaIsEstimable( prop( err, "sigmaInter" ), prop( err, "sigmaInterFixed" ) ) )
          prop( err, "sigmaInter" ),
        if ( .pfimSigmaIsEstimable( prop( err, "sigmaSlope" ), prop( err, "sigmaSlopeFixed" ) ) )
          prop( err, "sigmaSlope" )
      )
    } ) |> unlist( use.names = FALSE )
}


#' Shared ggplot2 theme for PFIM figures (16 pt).
#' @noRd
#' @keywords internal
.pfimBaseTheme = function( base_size = 16 ) {
  # Grey panel + 16 pt titles/ticks for report figures.
  theme_grey( base_size = base_size ) +
    theme(
      legend.position = "none",
      plot.title   = element_text( size = base_size, hjust = 0.5 ),
      axis.title.x = element_text( size = base_size ),
      axis.title.y = element_text( size = base_size ),
      axis.text.x  = element_text( size = base_size, angle = 90, vjust = 0.5 ),
      axis.text.y  = element_text( size = base_size, angle = 0, vjust = 0.5, hjust = 0.5 ),
      strip.text.x = element_text( size = base_size ),
      plot.margin  = margin( 10, 12, 10, 10 )
    )
}

#' Red secondary axis for design sampling times (response / SI plots).
#' @noRd
#' @keywords internal
.pfimSamplingAxisTheme = function( base_size = 16 ) {
  .pfimBaseTheme( base_size ) +
    theme(
      axis.title.x.top = element_text( color = "red", size = base_size, vjust = 2.0 ),
      axis.text.x.top  = element_text(
        angle = 90, hjust = 0, color = "red", size = base_size
      )
    )
}

#' SE/RSE bar-chart theme: 16 pt (args kept for callers).
#' @noRd
#' @keywords internal
.pfimSeRseTheme = function( n_facets = 5L, max_label_len = 12L ) {
  .pfimBaseTheme( 16 ) +
    theme( panel.spacing.x = unit( 0.6, "lines" ) )
}

#' \code{facet_wrap(space=)} exists only in ggplot2 >= 4.0.0; drop it on 3.5.x.
#' @noRd
#' @keywords internal
.facetWrapCompat = function( ... ) {
  args = list( ... )
  if ( !.pfimGgplot2SupportsFacetSpace() )
    args$space = NULL
  do.call( ggplot2::facet_wrap, args )
}

#' Faceted SE/RSE bar chart shared by FIM plot methods.
#' @noRd
#' @keywords internal
.fimSeRseBarPlot = function( df, metric, facet_levels ) {
  ggplot( df, aes( x = .data[["Parameter"]], y = .data[[ metric ]] ) ) +
    geom_col( position = "dodge", show.legend = FALSE ) +
    .facetWrapCompat(
      # free_x: each Greek category keeps its own scale; space= when ggplot2 >= 4.
      ~ factor( cat, levels = facet_levels ),
      scales   = "free_x",
      space    = "free_x",
      labeller = ggplot2::label_parsed
    ) +
    labs( x = "Parameter", y = metric ) +
    .pfimSeRseTheme( n_facets = length( facet_levels ) )
}

#' Print nested ggplot lists so knitr records one figure per plot.
#'
#' HTML reports store response / SI plots as design -> arm -> outcome (-> parameter)
#' lists. A bare \code{print(list)} does not draw grobs; this walks the tree.
#' @param x A ggplot, or a (possibly nested) list of ggplots.
#' @noRd
#' @keywords internal
.pfimPrintReportPlots = function( x ) {
  if ( inherits( x, "ggplot" ) ) {
    print( x )
  } else if ( is.list( x ) && !is.data.frame( x ) && length( x ) ) {
    walk( x, .pfimPrintReportPlots )
  }
  invisible( NULL )
}

#' @noRd
#' @keywords internal
.gradientMatrix = function( arm, cols, model = NULL ) {
  # Covariate/occasion models store gradients on the model, not the arm.
  mat = if ( !is.null( model ) && usesCovariateOccasionStructure( model ) ) {
    getArmEvaluationGradientsMatrix( model, arm )
  } else {
    raw = prop( arm, "evaluationGradients" )
    # May be one data.frame or a list of per-output frames to stack.
    df  = if ( is.data.frame( raw ) ) raw else do.call( rbind, raw )
    as.matrix( df )
  }
  mat[ , cols, drop = FALSE ]
}

#' Variance-parameter FIM block: \eqn{\frac12 \mathrm{Tr}(V^{-1} \mathrm{d}V_i V^{-1} \mathrm{d}V_j)}.
#' @noRd
#' @keywords internal
.computeMFVar_R = function( V_inv, dV_list ) {
  n = length( dV_list )
  if ( !n ) return( matrix( numeric( 0 ), 0L, 0L ) )
  # Precompute V^{-1} dV_i V^{-1}; Frobenius product with dV_j yields the (i,j) entry.
  T_mats = map( dV_list, ~ V_inv %*% .x %*% V_inv )
  outer(
    seq_len( n ), seq_len( n ),
    Vectorize( function( i, j ) 0.5 * sum( T_mats[[ i ]] * dV_list[[ j ]] ) )
  )
}

#' @rdname dot-computeMFVar_R
#' @noRd
#' @keywords internal
.computeMFVar = function( V_inv, dV_list ) {
  n = length( dV_list )
  if ( !n ) return( matrix( numeric( 0 ), 0L, 0L ) )
  # C++ kernel needs dense matrices; assembly may still hand sparse Diagonal blocks.
  computeMFVar_Rcpp( as.matrix( V_inv ), map( dV_list, as.matrix ) )
}
