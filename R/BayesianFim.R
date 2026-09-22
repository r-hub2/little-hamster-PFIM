#' @title BayesianFim
#' @description
#' Bayesian Fisher information matrix with shrinkage on random effects.
#'
#' First-order linearization on the random-effect scale. Covariate \code{beta} is not a
#' subject parameter (MAP conditions on \eqn{\beta}) and is omitted. Across
#' strata / arms, subject FIMs keep their prior and are aggregated by averaging
#' covariances: \eqn{\bar C = \sum w_s M_s^{-1}}, \eqn{M_{\mathrm{eff}}=\bar C^{-1}}.
#' A positive weight on a non-identifiable protocol yields \code{Inf} SE for the
#' null-space parameters. This covariance mixture is not comparable to summing
#' Fisher matrices (PFIM 6 / PopED on this case).
#'
#' With IOV (\eqn{\gamma>0}), each occasion uses the marginal FO residual
#' \eqn{V_k = R_k + F_k\,\mathrm{diag}(\gamma^2)\,F_k^\top} on the \eqn{\eta}-scale
#' (equivalent to the \eqn{\eta} block of the augmented \eqn{(\eta,\kappa)} system).
#'
#' @inheritParams Fim
#' @return A \code{BayesianFim} object (filled by \code{run()}).
#' @examples
#' \dontrun{
#' vignette("Example01")
#' }
#' @include Fim.R
#' @include pfim-fim-report-render.R
#' @include PFIMProject.R
#' @export

BayesianFim = new_class( "BayesianFim",
                         package    = "PFIM",
                         parent     = Fim
)
S4_register( BayesianFim )

# Estimable mu names with Greek prefix (setEvaluationFim, plots, ...).
#' Build estimable parameter labels for Bayesian outputs.
#' @param parameters List of \code{ModelParameter} objects.
#' @param greekPrefix Character prefix added to each parameter name.
#' @return Character vector of estimable parameter labels.
#' @noRd
#' @keywords internal
.bayesianEstimableParamNames = function( parameters, greekPrefix ) {
  .pfimBayesianEtaParameters( parameters ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    map_chr( ~ paste0( greekPrefix, .x ) )
}

#' PK mu/omega vectors for Bayesian shrinkage blocks.
#'
#' \code{adjustGradient(d, 1, mu)} yields the \eqn{M}-scale diagonal entry
#' (LogNormal -> \eqn{\mu}, Normal -> 1) without touching the data FIM yet.
#' @noRd
#' @keywords internal
.bayesianPkMuOmegaBlock = function( idx, feCols, paramByName ) {
  rows = map( idx, function( j ) {
    col   = feCols[[ j ]]
    pname = sub( "^mu_", "", col )
    p     = paramByName[[ pname ]]
    d     = prop( p, "distribution" )
    list(
      mu    = adjustGradient( d, 1, prop( d, "mu" ) ),
      omega = prop( d, "omega" ),
      fixed = .paramOmegaFixed( p )
    )
  } )
  list(
    mu    = map_dbl( rows, "mu" ),
    omega = map_dbl( rows, "omega" ),
    fixed = map_lgl( rows, "fixed" )
  )
}

#' Prior covariance \eqn{\Omega} on the random-effect (\eqn{\eta}) scale.
#'
#' For the Bayesian FIM, \eqn{\Omega = \mathrm{diag}(\omega^2)} even when
#' the data block is transformed by \eqn{M = \mathrm{diag}(\mu)} (LogNormal).
#' Always returns a matrix (including the 1x1 case) for \code{.safeSolve()}.
#' @noRd
#' @keywords internal
.bayesianOmega = function( omega_pk ) {
  n = length( omega_pk )
  if ( !n )
    return( matrix( numeric( 0L ), 0L, 0L ) )
  diag( omega_pk^2, nrow = n )
}

#' Prior precision \eqn{\Omega^{-1} = \mathrm{diag}(1/\omega^2)}.
#' Prefer this over \code{.safeSolve(.bayesianOmega(...))} - Omega is diagonal.
#' @noRd
#' @keywords internal
.bayesianOmegaInv = function( omega_pk ) {
  n = length( omega_pk )
  if ( !n )
    return( matrix( numeric( 0L ), 0L, 0L ) )
  diag( 1 / ( omega_pk^2 ), nrow = n )
}

#' Shrinkage (%) from PK FIM block and prior covariance on \eqn{\eta}.
#'
#' \eqn{W = M_{\mathrm{BF}}^{-1} \Omega^{-1}} (ratio of posterior
#' to prior variance on the random-effect scale), reported as percent.
#' @noRd
.bayesianShrinkageFromPrior = function( MF_pk, priorVariance ) {
  # priorVariance is Omega (diagonal); invert elementwise rather than via chol.
  omega_inv = if ( is.matrix( priorVariance ) && nrow( priorVariance ) == ncol( priorVariance ) )
    diag( 1 / diag( priorVariance ), nrow = nrow( priorVariance ) )
  else
    .safeCholInv( priorVariance )
  as.vector( diag( .safeCholInv( MF_pk ) %*% omega_inv ) * 100 )
}

#' @noRd
.bayesianShrinkageValues = function( shrinkage ) {
  if ( !is.matrix( shrinkage ) )
    return( as.numeric( shrinkage ) )
  if ( nrow( shrinkage ) == 1L )
    return( as.numeric( shrinkage[ 1L, , drop = TRUE ] ) )
  if ( ncol( shrinkage ) == 1L )
    return( as.numeric( shrinkage[ , 1L, drop = TRUE ] ) )
  as.numeric( shrinkage )
}

#' Console mu-labels aligned with stored Bayesian shrinkage.
#' @noRd
#' @keywords internal
.bayesianShrinkageMuLabels = function( shrinkage, evaluation ) {
  muPrefix = .greekConsole[ "mu" ]
  cols     = character( 0L )
  if ( is.matrix( shrinkage ) ) {
    if ( nrow( shrinkage ) == 1L && !is.null( colnames( shrinkage ) ) )
      cols = colnames( shrinkage )
    else if ( ncol( shrinkage ) == 1L && !is.null( rownames( shrinkage ) ) )
      cols = rownames( shrinkage )
  }
  if ( !length( cols ) || identical( cols, "Shrinkage" ) ) {
    n   = length( .bayesianShrinkageValues( shrinkage ) )
    est = .bayesianEstimableParamNames( prop( evaluation, "modelParameters" ), muPrefix )
    fe  = .fimFixedEffectLabels( evaluation )$columnNamesMu
    cols = if ( length( est ) == n ) est else if ( length( fe ) == n ) fe else est
  }
  .stripGreekPrefix( cols, muPrefix )
}

#' Map one shrinkage label onto an SE/RSE row (exact name, else mu-prefixed bare name).
#' @noRd
#' @keywords internal
.bayesianMatchShrinkageRow = function( name, rn, mu_prefix ) {
  j = match( name, rn )
  if ( !is.na( j ) )
    return( j )
  match( paste0( mu_prefix, .stripGreekPrefix( name, mu_prefix ) ), rn )
}

#' Align Bayesian shrinkage values to SE/RSE table rows (mu only; beta stays NA).
#' @noRd
#' @keywords internal
.bayesianShrinkageReportColumn = function( seDF, shrinkage ) {
  rn  = .pfimSeRownames( seDF, "BayesianFim tablesForReport" )
  out = rep( NA_real_, nrow( seDF ) )
  sh_vals = .bayesianShrinkageValues( shrinkage )
  if ( !length( sh_vals ) )
    return( out )

  greek  = .greekConsole
  mu_idx = which( startsWith( rn, greek[ "mu" ] ) )
  sh_mat = as.matrix( shrinkage )
  sh_names = colnames( sh_mat )

  if ( length( sh_names ) == length( sh_vals ) && length( sh_names ) ) {
    j = vapply(
      sh_names,
      function( nm ) .bayesianMatchShrinkageRow( nm, rn, greek[ "mu" ] ),
      integer( 1L )
    )
    ok = !is.na( j )
    out[ j[ ok ] ] = sh_vals[ ok ]
    unmatched = sum( !ok )
    if ( unmatched > 0L )
      warning(
        sprintf(
          "BayesianFim tablesForReport: %d shrinkage label(s) did not match FIM rows.",
          unmatched
        ),
        call. = FALSE
      )
  } else if ( length( sh_vals ) == length( mu_idx ) ) {
    out[ mu_idx ] = sh_vals
  } else {
    warning(
      sprintf(
        "BayesianFim tablesForReport: shrinkage length (%d) != mu rows (%d); leaving NA.",
        length( sh_vals ), length( mu_idx )
      ),
      call. = FALSE
    )
  }
  out
}

#' @noRd
.bayesianShrinkageMatrix = function( shrinkage, shrinkageCols ) {
  matrix(
    .bayesianShrinkageValues( shrinkage ),
    nrow = 1L,
    dimnames = list( "Shrinkage", shrinkageCols )
  )
}

#' PK prior covariance and shrinkage from a Bayesian FIM block.
#' @noRd
#' @keywords internal
.bayesianShrinkage = function( fisherMatrix, model, arm ) {

  parameters    = .pfimBayesianEtaParameters( prop( model, "modelParameters" ) )
  omega_pk      = vapply(
    parameters,
    function( p ) prop( prop( p, "distribution" ), "omega" ),
    numeric( 1L )
  )
  priorVariance = .bayesianOmega( omega_pk )
  n_eta         = length( omega_pk )
  MF_pk         = as.matrix( fisherMatrix )
  if ( nrow( MF_pk ) != n_eta || ncol( MF_pk ) != n_eta ) {
    # Legacy layouts with trailing beta / fixed columns: take the leading eta block.
    if ( nrow( MF_pk ) >= n_eta )
      MF_pk = MF_pk[ seq_len( n_eta ), seq_len( n_eta ), drop = FALSE ]
  }
  .bayesianShrinkageFromPrior( MF_pk, priorVariance )
}

#' Bayesian shrinkage on a design-level Fisher matrix.
#'
#' Shrinkage is always derived from \code{fisherMatrix} (the assembled design FIM
#' with prior once), so SE/RSE and shrinkage share one convention.
#' @param fisherMatrix Aggregated design Fisher matrix.
#' @param model \code{Model} used for PK layout.
#' @param evaluationArms Evaluated arms (first arm supplies column layout).
#' @return Numeric shrinkage vector (percent).
#' @noRd
#' @keywords internal
.bayesianShrinkageForDesign = function( fisherMatrix, model, evaluationArms ) {
  if ( !length( evaluationArms ) )
    return( numeric( 0 ) )
  .bayesianShrinkage( fisherMatrix, model, evaluationArms[[ 1L ]] )
}

#' Prior precision layout for the PK block of a Bayesian FIM (top-left eta block).
#' @return List with \code{n_pk} and \code{OmegaInv}, or \code{NULL} if no PK.
#' @noRd
#' @keywords internal
.bayesianPriorLayout = function( model, arm ) {
  parameters  = prop( model, "modelParameters" )
  feCols      = .fimFixedEffectColumnNames( model, arm )
  betaIdx     = which( startsWith( feCols, "beta_" ) )
  pkIdx       = if ( any( startsWith( feCols, "mu_" ) ) ) {
    which( startsWith( feCols, "mu_" ) )
  } else {
    seq_along( feCols )[ !seq_along( feCols ) %in% betaIdx ]
  }
  paramByName = set_names( parameters, map_chr( parameters, ~ prop( .x, "name" ) ) )
  pk          = .bayesianPkMuOmegaBlock( pkIdx, feCols, paramByName )
  keepPk      = !pk$fixed
  omega_pk    = pk$omega[ keepPk ]
  n_pk        = length( omega_pk )
  if ( !n_pk ) return( NULL )
  list( n_pk = n_pk, OmegaInv = .bayesianOmegaInv( omega_pk ) )
}

#' Strip one copy of \eqn{\Omega^{-1}} from the PK block (design multi-arm fix).
#' @noRd
#' @keywords internal
.bayesianStripPriorFromMatrix = function( M, OmegaInv, n_pk ) {
  if ( is.null( OmegaInv ) || !n_pk ) return( M )
  idx = seq_len( n_pk )
  M[ idx, idx ] = M[ idx, idx, drop = FALSE ] - OmegaInv
  M
}

#' Add one copy of \eqn{\Omega^{-1}} to the PK block.
#' @noRd
#' @keywords internal
.bayesianAddPriorToMatrix = function( M, OmegaInv, n_pk ) {
  if ( is.null( OmegaInv ) || !n_pk ) return( M )
  idx = seq_len( n_pk )
  M[ idx, idx ] = M[ idx, idx, drop = FALSE ] + OmegaInv
  M
}

#' Variance / data block of the Bayesian FIM (mu only; no beta, no sigma).
#'
#' Flat: \eqn{G^\top V^{-1} G}. Nested cov/IOV: harmonic mean of subject
#' FIMs (prior included per stratum).
#' @name evaluateVarianceFIM
#' @keywords internal

method( evaluateVarianceFIM, list( BayesianFim, Model, Arm ) ) = function( fim, model, arm ) {

  feCols = .pfimBayesianMuCols( model, arm )

  if ( .isNestedArmEvaluation( prop( arm, "evaluationGradients" ) ) ) {
    harm = .evaluateIndBayesMixtureFim( model, arm, feCols, bayesian = TRUE )
    return( list( MFbeta = harm$fisherMatrix, V = NULL, complete = TRUE ) )
  }

  gradient = .gradientMatrix( arm, feCols, model )
  V        = as.matrix( getArmEvaluationVarianceFlat( arm )$errorVariance )
  MFbeta   = crossprod( gradient, .safeCholInv( V ) ) %*% gradient

  list( MFbeta = MFbeta, V = V, complete = FALSE )
}

#' Compute the Bayesian FIM for one arm.
#'
#' FO Bayesian FIM: \eqn{M^\top M_{\mathrm{data}} M + \Omega^{-1}} for every
#' parameter with \eqn{\omega > 0} (\eqn{\beta} omitted; fix flags do not drop
#' eta). Shrinkage (\%) uses \eqn{M_{\mathrm{BF}}^{-1}\Omega^{-1}}.
#' @name evaluateFim
#' @keywords internal

method( evaluateFim, list( BayesianFim, Model, Arm ) ) = function( fim, model, arm ) {

  varBlock = evaluateVarianceFIM( fim, model, arm )
  if ( isTRUE( varBlock$complete ) ) {
    M = as.matrix( varBlock$MFbeta )
  } else {
    M = .pfimBayesianSubjectFim( as.matrix( varBlock$MFbeta ), model )
  }

  parameters    = .pfimBayesianEtaParameters( prop( model, "modelParameters" ) )
  omega_pk      = vapply( parameters, function( p ) prop( prop( p, "distribution" ), "omega" ), numeric( 1L ) )
  priorVariance = .bayesianOmega( omega_pk )

  prop( fim, "fisherMatrix" ) = M
  prop( fim, "shrinkage" )    = .bayesianShrinkageFromPrior( M, priorVariance )
  fim
}

#' Attach evaluated Bayesian FIM results (labels, SE/RSE, shrinkage matrix).
#'
#' SE/RSE use \code{.pfimBayesianSeRse}: for LogNormal, SE on \eqn{\theta} is
#' \eqn{\mu\cdot\mathrm{SE}_\eta} and RSE is \eqn{100\cdot\mathrm{SE}_\eta}.
#' @name setEvaluationFim
#' @keywords internal

method( setEvaluationFim, BayesianFim ) = function( fim, evaluation ) {

  parameters = prop( evaluation, "modelParameters" )
  greek      = .greekConsole
  allNames   = .bayesianEstimableParamNames( parameters, greek[ "mu" ] )

  fisherMatrix = prop( fim, "fisherMatrix" )
  if ( ncol( fisherMatrix ) != length( allNames ) )
    stop( sprintf(
      "BayesianFim setEvaluationFim: FIM dim %d != %d column names.",
      ncol( fisherMatrix ), length( allNames )
    ), call. = FALSE )
  dimnames( fisherMatrix ) = list( allNames, allNames )

  shrinkage = .bayesianShrinkageValues( prop( fim, "shrinkage" ) )
  se        = .pfimBayesianSeRse( fisherMatrix, parameters, allNames )
  .fimStoreEvaluationResult(
    fim, fisherMatrix, fisherMatrix, se,
    shrinkage = .bayesianShrinkageMatrix( shrinkage, allNames )
  )
}

#' Print FIM summaries to the console
#' @name showFIM
#' @export

method( showFIM, BayesianFim ) = function( fim ) {

  SEAndRSE               = prop( fim, "SEAndRSE" )
  fisherMatrix           = prop( fim, "fisherMatrix" )
  fixedEffects           = prop( fim, "fixedEffects" )
  shrinkage              = prop( fim, "shrinkage" )
  condNumberFixedEffects = prop( fim, "condNumberFixedEffects" )
  dcrit                  = Dcriterion( fim )

  cat( "\n*************************************** \n Bayesian Fisher Matrix \n*************************************** \n\n" )
  print( fisherMatrix )
  cat( "\n*************************************** \n Fixed effects \n*************************************** \n\n" )
  print( fixedEffects )
  cat( "\n*********************************************** \n Determinant, condition numbers and D-criterion \n*********************************************** \n\n" )
  cat( c( "log-Determinant:",  as.numeric( .fimLogDeterminant( fisherMatrix ) ) ), "\n" )
  cat( c( "D-criterion:",  as.numeric( dcrit              ) ), "\n" )
  cat( c( "Condition number of the fixed effects:", as.numeric( condNumberFixedEffects ), "\n" ) )
  .fimPrintSingularStatus( fim )
  cat( "\n*************************************** \n Shrinkage \n*************************************** \n\n" )
  print( shrinkage )
  cat( "\n*************************************** \n Parameters estimation \n*************************************** \n\n" )
  .printFimSeAndRse( fim )

  .fimShowSymbolLegend(
    rownames( prop( fim, "fisherMatrix" ) ),
    mu_label = "\u03bc  = fixed effects"
  )
  invisible( fim )
}

#' ggplot data for Bayesian SE/RSE barplots (mu and beta facets).
#' @param evaluation A \code{PFIMProject} object.
#' @param metric Character scalar, \code{"SE"} or \code{"RSE"}.
#' @return List with \code{data} (plot-ready \code{data.frame}) and \code{facet_levels}.
#' @noRd
#' @keywords internal
.bayesianSeRsePlotData = function( evaluation, metric ) {
  fim    = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  seDF   = prop( fim, "SEAndRSE" )$SEAndRSE
  greekC = .greekConsole
  rn     = .pfimSeRownames( seDF, "Bayesian SE/RSE plot" )

  # Labels follow FIM rows (omega-estimable mus), not all mu-estimable names.
  is_mu   = startsWith( rn, greekC[ "mu" ] )
  is_beta = startsWith( rn, greekC[ "beta" ] )
  y_col   = if ( metric == "SE" ) seDF$SE else seDF$RSE

  df = data.frame(
    Parameter = c(
      .stripGreekPrefix( rn[ is_mu ],   greekC[ "mu" ] ),
      .stripGreekPrefix( rn[ is_beta ], greekC[ "beta" ] )
    ),
    y = c( y_col[ is_mu ], y_col[ is_beta ] ),
    cat = c(
      rep( .pfimSeRseFacetLabel( metric, "mu"   ), sum( is_mu ) ),
      rep( .pfimSeRseFacetLabel( metric, "beta" ), sum( is_beta ) )
    ),
    stringsAsFactors = FALSE
  )
  names( df )[ 2L ] = metric
  list( data = df, facet_levels = unique( df$cat ) )
}

.bayesianSeRseBarPlot = function( evaluation, metric ) {
  px = .bayesianSeRsePlotData( evaluation, metric )
  .fimSeRseBarPlot( px$data, metric, px$facet_levels )
}

method( plotSEFIM, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation )
  .bayesianSeRseBarPlot( evaluation, "SE" )

method( plotRSEFIM, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation )
  .bayesianSeRseBarPlot( evaluation, "RSE" )

#' Default method for \code{BayesianFim}.
#' @param fim First argument of generic.
#' @param evaluation \code{PFIMProject} providing model parameter labels.
#' @return \code{ggplot} object showing Bayesian shrinkage by parameter.
#' @name plotShrinkage
#' @export
method( plotShrinkage, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  fim         = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  shrinkage   = prop( fim, "shrinkage" )
  paramLabels = .bayesianShrinkageMuLabels( shrinkage, evaluation )
  shrinkVals  = .bayesianShrinkageValues( shrinkage )

  data = data.frame( Parameter = paramLabels, Shrinkage = shrinkVals )
  ggplot( data, aes( x = Parameter, y = Shrinkage ) ) +
    geom_col( show.legend = FALSE ) +
    labs( x = "Parameter", y = "Shrinkage (%)" ) +
    .pfimBaseTheme()
}

#' FIM tables for HTML reports
#' @name tablesForReport
#' @keywords internal

method( tablesForReport, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {

  fim = setEvaluationFim( fim, evaluation )

  SEAndRSE               = prop( fim, "SEAndRSE" )$SEAndRSE
  fisherMatrix           = prop( fim, "fisherMatrix" )
  fixedEffects           = as.matrix( prop( fim, "fixedEffects" ) )
  shrinkage              = prop( fim, "shrinkage" )
  condNumberFixedEffects = prop( fim, "condNumberFixedEffects" )

  columnNamesFe = .fimFixedEffectLatexLabels( evaluation, fixedEffects = fixedEffects )
  colnames( fixedEffects ) = columnNamesFe
  rownames( fixedEffects ) = columnNamesFe
  shrink_col = .bayesianShrinkageReportColumn( SEAndRSE, shrinkage )

  list(
    fixedEffectsTable = .kblReportStyled( fixedEffects ),
    FIMCriteriaTable  = .fimCriteriaKable(
      .fimLogDeterminant( fisherMatrix ), Dcriterion( fim ), condNumberFixedEffects,
      singularFim = isTRUE( prop( fim, "singularFim" ) )
    ),
    SEAndRSETable     = .fimSeRseKable( columnNamesFe, SEAndRSE, shrink_col ),
    singularFimNote   = .fimSingularHtmlNote( fim ),
    singularFim       = isTRUE( prop( fim, "singularFim" ) )
  )
}

# Report rendering methods

#' Render the evaluation HTML report
#' @name generateReportEvaluation
#' @keywords internal
method( generateReportEvaluation, BayesianFim ) =
  .renderEvalReport( "EvaluationBayesianFIM.Rmd" )
