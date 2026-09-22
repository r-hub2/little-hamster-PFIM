#' @title PopulationFim
#' @description
#' Population Fisher information matrix for nonlinear mixed-effects models.
#'
#' @inheritParams Fim
#' @details
#' The population FIM has fixed-effects (\eqn{\mu}, \eqn{\beta}) and variance-effects
#' (\eqn{\omega^2}, \eqn{\gamma^2}, \eqn{\sigma}) blocks:
#' \deqn{M_P = \mathrm{bdiag}(N \cdot M_\mu,\; N \cdot M_\lambda)}
#' where \eqn{N} is the arm size and \eqn{M_\mu = G_\mu^\top V^{-1} G_\mu}.
#' @return A \code{PopulationFim} object (filled by \code{run()}).
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' S7::S7_inherits(prop(ev, "fim"), PopulationFim)
#' }
#' @include Fim.R
#' @include pfim-fim-report-render.R
#' @include PFIMProject.R
#' @export

PopulationFim = new_class( "PopulationFim", package = "PFIM", parent = Fim )
S4_register( PopulationFim )

#' Compute the population FIM for one arm.
#'
#' Builds fixed-effect and variance-effect blocks, scales both by arm size
#' \eqn{N}, and returns \code{bdiag(N M_beta, N M_lambda)}. Fixed-\eqn{\mu}
#' columns are dropped from the FE block; \code{fixedOmega} (not
#' \code{fixedMu}) drops \eqn{\omega^2} from the variance block.
#'
#' Scaling: both FE and VE blocks are \eqn{\times N} (iid subjects in the arm).
#' Individual FIM does not scale; Bayesian shrinkage is per-subject then aggregated
#' at the design level.
#' @name evaluateFim
#' @keywords internal
method( evaluateFim, list( PopulationFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameters          = prop( model, "modelParameters" )
  armSize             = prop( arm,   "size"            )
  hasComplexStructure = usesCovariateOccasionStructure( model )

  isFixedMu    = map_lgl( parameters, .paramMuFixed )
  isFixedOmega = map_lgl( parameters, .paramOmegaFixed )
  # Variance-block drop is omega-only: fixedMu must not remove omega².
  isFixedLambda = isFixedOmega

  result = .evaluateVarianceFIMPop( fim, model, arm,
                                    isFixedMu     = isFixedMu,
                                    isFixedOmega  = isFixedOmega,
                                    isFixedLambda = isFixedLambda )

  if ( !hasComplexStructure ) {
    # Trailing columns after omega are residual sigma; always keep them.
    nSigma     = ncol( result$MFVar ) - length( isFixedOmega )
    keepLambda = c( !isFixedLambda, rep( TRUE, nSigma ) )
    MFVar      = result$MFVar[ keepLambda, keepLambda, drop = FALSE ]
  } else {
    MFVar = result$MFVar
  }

  # M_P = bdiag(N M_beta, N M_lambda); N multiplies both blocks (same subject count).
  prop( fim, "fisherMatrix" ) = as.matrix(
    bdiag( result$MFbeta * armSize, MFVar * armSize )
  )
  fim
}

#' Greek-prefixed labels for estimable \code{omega} parameters.
#' @noRd
#' @keywords internal
.omegaNames = function( parameters, greek ) {
  .paramLabelNames( parameters, .paramOmegaEstimable, greek )
}

#' Greek-prefixed labels for estimable \code{gamma} parameters.
#' @noRd
#' @keywords internal
.gammaNames = function( parameters, greek ) {
  .paramLabelNames( parameters, .paramGammaEstimable, greek )
}

#' Build SE/RSE bar plots for population FIM blocks.
#' @param fim \code{PopulationFim} object.
#' @param evaluation \code{PFIMProject} with fitted model metadata.
#' @param metric Character scalar, \code{"SE"} or \code{"RSE"}.
#' @return \code{ggplot} bar chart object.
#' @noRd
#' @keywords internal
.plotFimBars = function( fim, evaluation, metric ) {
  fim    = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  seDF   = prop( fim, "SEAndRSE" )$SEAndRSE
  greekC = .greekConsole
  facet  = function( key ) .pfimSeRseFacetLabel( metric, key )
  rn     = .pfimSeRownames( seDF, "plotFimBars" )

  # SE/RSE barplots: mu / beta / omega / gamma / sigma (sensitivity plots omit beta).
  keep = startsWith( rn, greekC[ "mu" ] ) |
    startsWith( rn, greekC[ "beta" ] ) |
    startsWith( rn, greekC[ "omega" ] ) |
    startsWith( rn, greekC[ "gamma" ] ) |
    startsWith( rn, greekC[ "sigma" ] )
  rn_k   = rn[ keep ]
  y_vals = if ( metric == "SE" ) seDF$SE[ keep ] else seDF$RSE[ keep ]
  params = c(
    .stripGreekPrefix( rn_k[ startsWith( rn_k, greekC[ "mu" ]    ) ], greekC[ "mu" ]    ),
    .stripGreekPrefix( rn_k[ startsWith( rn_k, greekC[ "beta" ]  ) ], greekC[ "beta" ]  ),
    .stripGreekPrefix( rn_k[ startsWith( rn_k, greekC[ "omega" ] ) ], greekC[ "omega" ] ),
    .stripGreekPrefix( rn_k[ startsWith( rn_k, greekC[ "gamma" ] ) ], greekC[ "gamma" ] ),
    .stripGreekPrefix( rn_k[ startsWith( rn_k, greekC[ "sigma" ] ) ], greekC[ "sigma" ] )
  )
  cats   = c(
    rep( facet( "mu"    ), sum( startsWith( rn_k, greekC[ "mu" ]    ) ) ),
    rep( facet( "beta"  ), sum( startsWith( rn_k, greekC[ "beta" ]  ) ) ),
    rep( facet( "omega" ), sum( startsWith( rn_k, greekC[ "omega" ] ) ) ),
    rep( facet( "gamma" ), sum( startsWith( rn_k, greekC[ "gamma" ] ) ) ),
    rep( facet( "sigma" ), sum( startsWith( rn_k, greekC[ "sigma" ] ) ) )
  )
  df = data.frame(
    Parameter = params,
    y         = y_vals,
    cat       = cats,
    stringsAsFactors = FALSE
  )
  names( df )[ 2L ] = metric

  facet_levels = c(
    facet( "mu" ),
    if ( any( startsWith( rn_k, greekC[ "beta" ] ) ) ) facet( "beta" ),
    facet( "omega" ),
    if ( any( startsWith( rn_k, greekC[ "gamma" ] ) ) ) facet( "gamma" ),
    facet( "sigma" )
  )
  .fimSeRseBarPlot( df, metric, facet_levels )
}

#' Attach evaluated population FIM results to a project.
#'
#' Labels fixed/variance blocks with Greek console prefixes, splits SE/RSE, and
#' stores condition numbers for report / \code{showFIM} consumers.
#'
#' Column contract: rownames/colnames must match
#' \code{c(mu, beta, omega, gamma?, sigma)} in that order, same length as
#' \code{ncol(fisherMatrix)}. RSE uses \code{abs(value)} denominators
#' (\code{abs_denominator = TRUE}) so negative betas stay finite.
#' @name setEvaluationFim
#' @keywords internal
method( setEvaluationFim, PopulationFim ) = function( fim, evaluation ) {

  parameters = prop( evaluation, "modelParameters" )
  modelError = prop( evaluation, "modelError"       )
  greek      = .greekConsole

  has_IOV = .hasIovParameters( parameters )
  fe      = .fimFixedEffectLabels( evaluation, greek )

  columnNamesMu    = .estimableMuNames( parameters, greek[ "mu" ] )
  columnNamesBeta  = fe$columnNamesBeta
  muValues         = fe$muValues
  betaValues       = fe$betaValues
  columnNamesOmega = .omegaNames( parameters, greek[ "omega" ] )
  columnNamesGamma = if ( has_IOV ) .gammaNames( parameters, greek[ "gamma" ] ) else character( 0L )
  columnNamesSigma = .sigmaNames( modelError, greek[ "sigma" ] )
  # VE values are omega² / gamma² (variance scale), matching FIM parameterisation.
  omegaValues      = .paramOmegaSqValues( parameters )
  gammaValues      = if ( has_IOV ) .paramGammaSqValues( parameters ) else numeric( 0L )
  sigmaValues      = .sigmaValues( modelError )

  M        = prop( fim, "fisherMatrix" )
  allNames = c( columnNamesMu, columnNamesBeta,
                columnNamesOmega, columnNamesGamma, columnNamesSigma )

  if ( ncol( M ) != length( allNames ) )
    stop( sprintf(
      "setEvaluationFim: FIM dim %d \u2260 %d column names.\nNames: %s",
      ncol( M ), length( allNames ), paste( allNames, collapse = ", " )
    ))

  dimnames( M ) = list( allNames, allNames )
  feNames       = c( columnNamesMu, columnNamesBeta )
  veNames       = c( columnNamesOmega, columnNamesGamma, columnNamesSigma )
  fixedEffects  = M[ feNames, feNames, drop = FALSE ]
  varEffects    = M[ veNames, veNames, drop = FALSE ]

  pVals = c( muValues, betaValues, omegaValues, gammaValues, sigmaValues )
  se    = .fimBuildSeAndRse( M, allNames, pVals, abs_denominator = TRUE )
  .fimStoreEvaluationResult(
    fim, M, fixedEffects, se, varianceEffects = varEffects
  )
}

#' Print population FIM sections to the console.
#' @param fim \code{PopulationFim} object.
#' @param evaluation Optional \code{PFIMProject} used to refresh summaries.
#' @return Invisibly returns \code{fim}.
#' @noRd
#' @keywords internal
.showPopulationFimConsole = function( fim, evaluation = NULL ) {

  .hdr = function( t ) { cat( "\n*************************************** \n" )
    cat( " ", t, "\n" )
    cat( "*************************************** \n\n" ) }

  .hdr( "Population Fisher Matrix" );                     print( prop( fim, "fisherMatrix"   ) )
  .hdr( "Fixed effects (\u03bc)" );                       print( prop( fim, "fixedEffects"   ) )
  .hdr( "Variance components (\u03c9\u00B2, \u03b3\u00B2, \u03c3)" )
  print( prop( fim, "varianceEffects" ) )

  cat( "\n********************************************* \n",
       " log-Determinant, condition numbers and D-criterion \n",
       "*********************************************** \n\n" )
  cat( "log-Determinant:",  .fimLogDeterminant( prop( fim, "fisherMatrix" ) ), "\n" )
  cat( "D-criterion:",  Dcriterion( fim ),                  "\n" )
  cat( "Condition number (fixed effects):",      prop( fim, "condNumberFixedEffects"    ), "\n" )
  cat( "Condition number (variance components):", prop( fim, "condNumberVarianceEffects"), "\n" )
  .fimPrintSingularStatus( fim )

  .hdr( "Parameters estimation" )
  .printFimSeAndRse( fim, evaluation )

  .fimShowSymbolLegend(
    rownames( prop( fim, "fisherMatrix" ) ),
    include_omega = TRUE,
    include_sigma = TRUE
  )
  invisible( fim )
}

#' Print FIM summaries to the console
#' @name showFIM
#' @export
method( showFIM, PopulationFim ) = function( fim ) {
  .showPopulationFimConsole( fim, evaluation = NULL )
}

method( plotSEFIM,  list( PopulationFim, PFIMProject ) ) =
  function( fim, evaluation ) .plotFimBars( fim, evaluation, "SE"  )

method( plotRSEFIM, list( PopulationFim, PFIMProject ) ) =
  function( fim, evaluation ) .plotFimBars( fim, evaluation, "RSE" )

#' FIM tables for HTML reports
#' @name tablesForReport
#' @keywords internal
method( tablesForReport, list( PopulationFim, PFIMProject ) ) = function( fim, evaluation ) {
  .fimTablesForReportStandard( fim, evaluation )
}

#' Render the evaluation HTML report
#' @name generateReportEvaluation
#' @keywords internal
method( generateReportEvaluation, PopulationFim ) =
  .renderEvalReport( "EvaluationPopulationFIM.Rmd" )
