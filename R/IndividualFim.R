#' @title IndividualFim
#' @description
#' Individual (subject-level) Fisher information matrix.
#'
#' Supports covariates and IOV: gradients and residual variance are averaged over
#' covariate combinations and occasions (same expectation as for population designs).
#'
#' @inheritParams Fim
#' @details
#' Block-diagonal structure:
#' \deqn{M_I = \mathrm{bdiag}(M_\mu, M_\sigma)}
#' where \eqn{M_\mu = G^\top V^{-1} G} and \eqn{M_\sigma} is the variance-effects block.
#' @include Fim.R
#' @export

IndividualFim = new_class( "IndividualFim", package = "PFIM", parent = Fim )
S4_register( IndividualFim )

#' Variance block of the population FIM
#' @name evaluateVarianceFIM
#' @export

method( evaluateVarianceFIM, list( IndividualFim, Model, Arm ) ) = function( fim, model, arm ) {

  ev    = getArmEvaluationVarianceFlat( arm )
  V     = as.matrix( ev$errorVariance )
  V_inv = .safeCholInv( V )

  list( MFVar = .computeMFVar( V_inv, ev$sigmaDerivatives ), V = V )
}

#' Compute the FIM for one arm
#' @name evaluateFim
#' @export

method( evaluateFim, list( IndividualFim, Model, Arm ) ) = function( fim, model, arm ) {

  ev    = getArmEvaluationVarianceFlat( arm )
  V     = as.matrix( ev$errorVariance )
  V_inv = .safeCholInv( V )
  MFVar = .computeMFVar( V_inv, ev$sigmaDerivatives )

  feCols = .fimFixedEffectColumnNames( model, arm )
  G      = .gradientMatrix( arm, feCols, model )
  MFbeta = crossprod( G, V_inv ) %*% G

  prop( fim, "fisherMatrix" ) = as.matrix( bdiag( MFbeta, MFVar ) )
  fim
}

# setOptimalArms
#' Store optimal arm designs on the FIM object
#' @name setOptimalArms
#' @export
method( setOptimalArms, list( IndividualFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm )
    .setOptimalArmsMultiplicative( optimizationAlgorithm )

#' Store optimal arm designs on the FIM object
#' @name setOptimalArms
#' @export
method( setOptimalArms, list( IndividualFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm )
    .setOptimalArmsFedorovWynn( optimizationAlgorithm )

#' Attach evaluated FIM results to a project
#' @name setEvaluationFim
#' @export

method( setEvaluationFim, IndividualFim ) = function( fim, evaluation ) {

  greek  = .greekConsole
  fe     = .fimFixedEffectLabels( evaluation, greek )
  sigma  = .fimSigmaBlockLabels( evaluation, greek )

  allNames = c( fe$columnNamesMu, fe$columnNamesBeta, sigma$columnNamesSigma )
  pVals    = c( fe$muValues, fe$betaValues, sigma$sigmaValues )

  M = prop( fim, "fisherMatrix" )

  # if ( ncol( M ) != length( allNames ) )
  #   stop( sprintf(
  #     "IndividualFim setEvaluationFim: FIM dim %d != %d column names.",
  #     ncol( M ), length( allNames )
  #   ), call. = FALSE )

  dimnames( M ) = list( allNames, allNames )

  feNames         = c( fe$columnNamesMu, fe$columnNamesBeta )
  fixedEffects    = M[ feNames, feNames, drop = FALSE ]
  varianceEffects = M[ sigma$columnNamesSigma, sigma$columnNamesSigma, drop = FALSE ]

  se = .fimBuildSeAndRse( M, allNames, pVals, abs_denominator = FALSE )

  prop( fim, "fisherMatrix"              ) = M
  prop( fim, "fixedEffects"              ) = fixedEffects
  prop( fim, "varianceEffects"           ) = varianceEffects
  prop( fim, "condNumberFixedEffects"    ) = .conditionNumber( fixedEffects )
  prop( fim, "condNumberVarianceEffects" ) = .conditionNumber( varianceEffects )
  prop( fim, "SEAndRSE"                  ) = se

  fim
}

#' Print FIM summaries to the console
#' @name showFIM
#' @export

method( showFIM, IndividualFim ) = function( fim ) {

  .hdr = function( t ) cat( sprintf(
    "\n*************************************** \n %s \n*************************************** \n\n", t
  ))

  .hdr( "Individual Fisher Matrix" );        print( prop( fim, "fisherMatrix"   ) )
  .hdr( "Fixed effects (\u03bc)" );          print( prop( fim, "fixedEffects"   ) )
  .hdr( "Variance components (\u03c3)" );    print( prop( fim, "varianceEffects") )

  M   = prop( fim, "fisherMatrix" )
  cn1 = prop( fim, "condNumberFixedEffects"    )
  cn2 = prop( fim, "condNumberVarianceEffects" )

  cat( "\n*********************************************** \n",
       " Determinant, condition numbers and D-criterion \n",
       "*********************************************** \n\n" )
  cat( "Determinant:",  as.numeric( det( M ) ),       "\n" )
  cat( "D-criterion:",  as.numeric( Dcriterion(fim) ), "\n" )
  cat( "Conditional number (fixed effects):",    as.numeric( cn1 ), "\n" )
  cat( "Conditional number (variance effects):", as.numeric( cn2 ), "\n" )
  .hdr( "Parameters estimation" )
  .printFimSeAndRse( fim )

  invisible( fim )
}

# SE / RSE bar charts
.individualSEPlot = function( fim, evaluation, metric ) {

  parameters = prop( evaluation, "modelParameters" )
  modelError = prop( evaluation, "modelError"       )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se         = prop( fim, "SEAndRSE" )
  greek      = .greekPlot   # defined in Fim.R

  paramsMu = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) )

  paramsSigma = modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "_inter_", out ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "_slope_", out )
      )
    }) |> unlist( use.names = FALSE )

  y_vals = if ( metric == "SE" ) se$SE$SE else se$RSE$RSE
  cats   = paste0( metric, " ",
                   c( rep( greek[ "mu"    ], length( paramsMu    ) ),
                      rep( greek[ "sigma" ], length( paramsSigma ) ) ) )

  df = data.frame( Parameter = c( paramsMu, paramsSigma ), y = y_vals, cat = cats )
  names( df )[ 2L ] = metric

  ggplot( df, aes( x = .data[["Parameter"]], y = .data[[ metric ]] ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~ factor( cat, levels = unique( cats ) ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

#' Barplot of standard errors from a FIM
#' @name plotSEFIM
#' @export
method( plotSEFIM,  list( IndividualFim, PFIMProject ) ) =
  function( fim, evaluation ) .individualSEPlot( fim, evaluation, "SE"  )

#' Barplot of relative standard errors from a FIM
#' @name plotRSEFIM
#' @export
method( plotRSEFIM, list( IndividualFim, PFIMProject ) ) =
  function( fim, evaluation ) .individualSEPlot( fim, evaluation, "RSE" )

#' FIM tables for HTML reports
#' @name tablesForReport
#' @export

method( tablesForReport, list( IndividualFim, PFIMProject ) ) = function( fim, evaluation ) {

  parameters                = prop( evaluation, "modelParameters" )
  modelError                = prop( evaluation, "modelError"       )
  SEAndRSE                  = prop( fim, "SEAndRSE" )$SEAndRSE
  M                         = prop( fim, "fisherMatrix"            )
  fe                        = as.matrix( prop( fim, "fixedEffects"    ) )
  ve                        = as.matrix( prop( fim, "varianceEffects" ) )
  condNumberFixedEffects    = prop( fim, "condNumberFixedEffects"    )
  condNumberVarianceEffects = prop( fim, "condNumberVarianceEffects" )
  greek                     = .greekLatex   # defined in Fim.R

  columnNamesMu = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~  prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ paste0( greek[ "mu" ], prop( .x, "name" ), "}$" ) )

  columnNamesSigma = modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "{inter}_{", out, "}$" ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "{slope}_{", out, "}$" )
      )
    }) |> unlist( use.names = FALSE )

  dimnames( fe ) = list( columnNamesMu,    columnNamesMu    )
  dimnames( ve ) = list( columnNamesSigma, columnNamesSigma )

  .kbl_styled = function( df )
    kbl( df ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  FIMCriteriaTable = data.frame(
    determinant               = det( M ),
    dcriterion                = Dcriterion( fim ),
    condNumberFixedEffects    = condNumberFixedEffects,
    condNumberVarianceEffects = condNumberVarianceEffects
  ) |>
    kbl( col.names = c( "", "", "Fixed effects", "Variance effects" ),
         align = "c", format = "html" ) |>
    add_header_above( c( "Determinant" = 1, "D-criterion" = 1,
                         "Condition number" = 2 ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  SEAndRSETable = data.frame(
    c( columnNamesMu, columnNamesSigma ), round( SEAndRSE, 3 )
  ) |>
    (\( df ) { row.names( df ) = NULL; df })() |>
    kbl( col.names = c( "Parameters", "Parameter values", "SE", "RSE (%)"),
         align = "c" ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  list( fixedEffectsTable    = .kbl_styled( fe ),
        varianceEffectsTable = .kbl_styled( ve ),
        FIMCriteriaTable     = FIMCriteriaTable,
        SEAndRSETable        = SEAndRSETable )
}

# Report rendering
# .renderReport and .renderEvalReport are defined in Fim.R.
# .reportIndividualPath is REMOVED â€” .reportTemplatePath (Fim.R) replaces it.

#' Render the evaluation HTML report
#' @name generateReportEvaluation
#' @export
method( generateReportEvaluation, IndividualFim ) =
  .renderEvalReport( "EvaluationIndividualFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( IndividualFim, MultiplicativeAlgorithm ) ) =
  .renderReport( "OptimizationMultiplicativeAlgorithmIndividualFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( IndividualFim, FedorovWynnAlgorithm ) ) =
  .renderReport( "OptimizationFedorovWynnAlgorithmIndividualFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( IndividualFim, SimplexAlgorithm ) ) =
  .renderReport( "OptimizationSimplexAlgorithmIndividualFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( IndividualFim, PSOAlgorithm ) ) =
  .renderReport( "OptimizationPSOAlgorithmIndividualFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( IndividualFim, PGBOAlgorithm ) ) =
  .renderReport( "OptimizationPGBOAlgorithmIndividualFIM.Rmd" )
