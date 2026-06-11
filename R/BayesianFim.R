#' @title BayesianFim
#' @description
#' Bayesian Fisher information matrix with shrinkage on random effects.
#'
#' Supports covariates and IOV via expectation over combinations and occasions.
#' PK parameters use the Bayesian shrinkage formula; covariate \code{beta} effects
#' enter as extra fixed-effect columns without IIV shrinkage.
#'
#' @inheritParams Fim
#' @include Fim.R
#' @include MultiplicativeAlgorithm.R
#' @include FedorovWynnAlgorithm.R
#' @export

BayesianFim = new_class( "BayesianFim",
                         package    = "PFIM",
                         parent     = Fim
)
S4_register( BayesianFim )

# â”€â”€ Private helper â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# Filter estimable (non-fixed, non-zero) parameters and return their names
# with the console Greek-mu prefix â€” shared by setEvaluationFim, plotSEFIM, etc.
.bayesianEstimableParamNames = function( parameters, greekPrefix ) {
  parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    map_chr( ~ paste0( greekPrefix, .x ) )
}

#' Variance block of the population FIM
#' @name evaluateVarianceFIM
#' @export

method( evaluateVarianceFIM, list( BayesianFim, Model, Arm ) ) = function( fim, model, arm ) {

  feCols   = .fimFixedEffectColumnNames( model, arm )
  gradient = .gradientMatrix( arm, feCols, model )
  V        = as.matrix( getArmEvaluationVarianceFlat( arm )$errorVariance )
  MFbeta   = crossprod( gradient, .safeCholInv( V ) ) %*% gradient

  list( MFbeta = MFbeta, V = V )
}

#' Compute the FIM for one arm
#' @name evaluateFim
#' @export

method( evaluateFim, list( BayesianFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameters = prop( model, "modelParameters" )
  feCols     = .fimFixedEffectColumnNames( model, arm )
  gradient   = .gradientMatrix( arm, feCols, model )
  V          = as.matrix( getArmEvaluationVarianceFlat( arm )$errorVariance )
  V_inv      = .safeCholInv( V )

  MFbeta_full = crossprod( gradient, V_inv ) %*% gradient

  cols    = feCols
  betaIdx = which( startsWith( cols, "beta_" ) )
  pkIdx   = if ( any( startsWith( cols, "mu_" ) ) ) {
    which( startsWith( cols, "mu_" ) )
  } else {
    seq_along( cols )[ !seq_along( cols ) %in% betaIdx ]
  }

  paramByName = set_names( parameters, map_chr( parameters, ~ prop( .x, "name" ) ) )

  buildPkMuOmega = function( idx ) {

    mu_v = omega_v = numeric( length( idx ) )
    fixed = logical( length( idx ) )

    for ( i in seq_along( idx ) ) {
      j     = idx[[ i ]]
      col   = cols[[ j ]]
      pname = sub( "^mu_", "", col )
      p     = paramByName[[ pname ]]
      d     = prop( p, "distribution" )

      mu_v[[ i ]]    = if ( S7::S7_inherits( d, Normal ) ) 1 else prop( d, "mu" )
      omega_v[[ i ]] = prop( d, "omega" )
      fixed[[ i ]]   = isTRUE( prop( p, "fixedMu" ) ) || isTRUE( prop( p, "fixedOmega" ) ) ||
        mu_v[[ i ]] == 0 || omega_v[[ i ]] == 0
    }
    list( mu = mu_v, omega = omega_v, fixed = fixed )
  }

  if ( length( betaIdx ) > 0L ) {
    pk  = buildPkMuOmega( pkIdx )
    keepPk = !pk$fixed

    pkIdxKeep = pkIdx[ keepPk ]
    mu_pk     = pk$mu[ keepPk ]
    omega_pk  = pk$omega[ keepPk ]

    MF_pk = MFbeta_full[ pkIdxKeep, pkIdxKeep, drop = FALSE ]
    mu_d  = if ( length( mu_pk ) == 1L ) mu_pk[[ 1L ]] else diag( mu_pk )
    om_d  = if ( length( omega_pk ) == 1L ) omega_pk[[ 1L ]]^2 else diag( omega_pk^2 )

    priorVariance = mu_d %*% om_d %*% mu_d
    MF_bayes_pk   = t( mu_d ) %*% MF_pk %*% mu_d + .safeSolve( priorVariance )

    MF_beta = MFbeta_full[ betaIdx, betaIdx, drop = FALSE ]
    MFbeta  = as.matrix( bdiag( MF_bayes_pk, MF_beta ) )

    prop( fim, "shrinkage" ) = as.vector( diag( .safeCholInv( MF_bayes_pk ) %*% .safeCholInv( priorVariance ) ) * 100 )
  } else {
    pk = buildPkMuOmega( pkIdx )
    indexFixed = which( pk$fixed )
    mu_vec    = pk$mu
    omega_vec = pk$omega

    if ( length( indexFixed ) > 0L ) {
      mu_vec    = mu_vec[ -indexFixed ]
      omega_vec = omega_vec[ -indexFixed ]
      MFbeta_full = MFbeta_full[ -indexFixed, -indexFixed, drop = FALSE ]
    }

    mu    = if ( length( mu_vec ) == 1L ) mu_vec[[ 1L ]] else diag( mu_vec )
    omega = if ( length( omega_vec ) == 1L ) omega_vec[[ 1L ]]^2 else diag( omega_vec^2 )

    priorVariance = mu %*% omega %*% mu
    MFbeta        = t( mu ) %*% MFbeta_full %*% mu + .safeSolve( priorVariance )
    prop( fim, "shrinkage" ) = as.vector( diag( .safeCholInv( MFbeta ) %*% .safeCholInv( priorVariance ) ) * 100 )
  }
  prop( fim, "fisherMatrix" ) = MFbeta
  fim
}

#' Store optimal arm designs on the FIM object
#' @name setOptimalArms
#' @export

method( setOptimalArms, list( BayesianFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {

    out             = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
    armFims         = out$armFims
    weights         = out$multiplicativeAlgorithmOutput[[ "weights" ]]
    weightsIndex    = which( weights > out$weightThreshold )

    armList = map( weightsIndex, function( idx ) {
      arm = pluck( armFims[[ idx ]], 1L )
      prop( arm, "size" ) = 1.0
      prop( arm, "name" ) = paste0( "Arm", idx )
      arm
    })

    armList[ rev( order( map_dbl( armList, ~ prop( .x, "size" ) ) ) ) ]
  }

#' Store optimal arm designs on the FIM object
#' @name setOptimalArms
#' @export

method( setOptimalArms, list( BayesianFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {
    out = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
    imap( out$listArms, ~ {
      prop( .x$arm, "name" ) = paste0( "Arm", .y )
      prop( .x$arm, "size" ) = 1.0
      .x
    })
  }

#' Attach evaluated FIM results to a project
#' @name setEvaluationFim
#' @export

method( setEvaluationFim, BayesianFim ) = function( fim, evaluation ) {

  parameters = prop( evaluation, "modelParameters" )
  greek      = .greekConsole
  fe         = .fimFixedEffectLabels( evaluation, greek )

  useCovBeta = fe$has_cov && length( fe$feInternal ) > 0L &&
    any( startsWith( fe$feInternal, "beta_" ) )

  if ( useCovBeta ) {
    allNames      = c( fe$columnNamesMu, fe$columnNamesBeta )
    pVals         = c( fe$muValues, fe$betaValues )
    shrinkageCols = fe$columnNamesMu
  } else {
    allNames      = .bayesianEstimableParamNames( parameters, greek[ "mu" ] )
    pVals         = parameters |>
      keep( ~ !isTRUE( prop( .x, "fixedMu"    ) ) ) |>
      keep( ~ prop( prop( .x, "distribution" ), "mu"    ) != 0 ) |>
      keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
      keep( ~ prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
      map_dbl( ~ prop( prop( .x, "distribution" ), "mu" ) )
    shrinkageCols = allNames[ str_starts( allNames, greek[ "mu" ] ) ]
  }

  fisherMatrix = prop( fim, "fisherMatrix" )
  colnames( fisherMatrix ) = allNames
  rownames( fisherMatrix ) = allNames
  fixedEffects = fisherMatrix

  shrinkage = prop( fim, "shrinkage" )
  se        = .fimBuildSeAndRse( fisherMatrix, allNames, pVals, abs_denominator = TRUE )

  prop( fim, "fisherMatrix"           ) = fisherMatrix
  prop( fim, "fixedEffects"           ) = fixedEffects
  prop( fim, "shrinkage"              ) = t( matrix(
    shrinkage, nrow = 1L,
    dimnames = list( "Shrinkage", shrinkageCols )
  ) )
  prop( fim, "condNumberFixedEffects" ) = .conditionNumber( fixedEffects )
  prop( fim, "SEAndRSE"               ) = se

  fim
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
  cat( c( "Determinant:",  as.numeric( det( fisherMatrix ) ) ), "\n" )
  cat( c( "D-criterion:",  as.numeric( dcrit              ) ), "\n" )
  cat( c( "Conditional number of the fixed effects:", as.numeric( condNumberFixedEffects ), "\n" ) )
  cat( "\n*************************************** \n Shrinkage \n*************************************** \n\n" )
  print( shrinkage )
  cat( "\n*************************************** \n Parameters estimation \n*************************************** \n\n" )
  .printFimSeAndRse( fim )

  invisible( fim )
}

#' Barplot of standard errors from a FIM
#' @name plotSEFIM
#' @export

method( plotSEFIM, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  parameters = prop( evaluation, "modelParameters" )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se         = prop( fim, "SEAndRSE" )
  greek      = .greekPlot

  params  = .bayesianEstimableParamNames( parameters, "" )  # raw names
  data    = data.frame(
    Parameter = params,
    SE        = se$SE$SE,
    cat       = paste0( "SE ", greek[ "mu" ] )
  )
  ggplot( data, aes( x = Parameter, y = SE ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~factor( cat, levels = paste0( "SE ", greek[ "mu" ] ) ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title   = element_text( size = 16, hjust = 0.5 ),
           axis.text.x  = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

method( plotRSEFIM, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  parameters = prop( evaluation, "modelParameters" )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se         = prop( fim, "SEAndRSE" )
  greek      = .greekPlot

  params  = .bayesianEstimableParamNames( parameters, "" )
  data    = data.frame(
    Parameter = params,
    RSE       = se$RSE$RSE,
    cat       = paste0( "RSE ", greek[ "mu" ] )
  )
  ggplot( data, aes( x = Parameter, y = RSE ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~factor( cat, levels = paste0( "RSE ", greek[ "mu" ] ) ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

method( plotShrinkage, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  parameters = prop( evaluation, "modelParameters" )
  shrinkage  = prop( fim, "shrinkage" )
  params     = .bayesianEstimableParamNames( parameters, "" )

  data = data.frame( Parameter = params, Shrinkage = as.vector( shrinkage ) )
  ggplot( data, aes( x = Parameter, y = Shrinkage ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

#' FIM tables for HTML reports
#' @name tablesForReport
#' @export

method( tablesForReport, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {

  parameters             = prop( evaluation, "modelParameters" )
  SEAndRSE               = prop( fim, "SEAndRSE" )$SEAndRSE
  fisherMatrix           = prop( fim, "fisherMatrix"   )
  fixedEffects           = as.matrix( prop( fim, "fixedEffects" ) )
  shrinkage              = prop( fim, "shrinkage" )
  condNumberFixedEffects = prop( fim, "condNumberFixedEffects" )
  greek                  = .greekLatex

  columnNamesMu = .bayesianEstimableParamNames( parameters, greek[ "mu" ] ) |>
    map_chr( ~ paste0( .x, "}$" ) )

  colnames( fixedEffects ) = columnNamesMu
  rownames( fixedEffects ) = columnNamesMu

  fixedEffectsTable = fixedEffects |>
    kbl() |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  FIMCriteriaTable = data.frame(
    determinant  = det( fisherMatrix ),
    dcriterion   = Dcriterion( fim ),
    FixedEffects = condNumberFixedEffects
  ) |>
    kbl( col.names = c( "", "", "Fixed effects" ), align = "c", format = "html" ) |>
    add_header_above( c( "Determinant" = 1, "D-criterion" = 1, "Condition number" = 1 ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  SEAndRSETable = data.frame(
    c( columnNamesMu ), round( SEAndRSE, 3 ), as.vector( shrinkage )
  ) |>
    (\( df ) { row.names( df ) = NULL; df })() |>
    kbl( col.names = c( "Parameters", "Parameter values", "SE", "RSE (%)", "Shrinkage" ),
         align = "c" ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  list( fixedEffectsTable = fixedEffectsTable,
        FIMCriteriaTable  = FIMCriteriaTable,
        SEAndRSETable     = SEAndRSETable )
}

# Report rendering methods

#' Render the evaluation HTML report
#' @name generateReportEvaluation
#' @export
method( generateReportEvaluation, BayesianFim ) =
  .renderEvalReport( "EvaluationBayesianFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( BayesianFim, MultiplicativeAlgorithm ) ) =
  .renderReport( "OptimizationMultiplicativeAlgorithmBayesianFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( BayesianFim, FedorovWynnAlgorithm ) ) =
  .renderReport( "OptimizationFedorovWynnAlgorithmBayesianFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( BayesianFim, SimplexAlgorithm ) ) =
  .renderReport( "OptimizationSimplexAlgorithmBayesianFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( BayesianFim, PSOAlgorithm ) ) =
  .renderReport( "OptimizationPSOAlgorithmBayesianFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( BayesianFim, PGBOAlgorithm ) ) =
  .renderReport( "OptimizationPGBOAlgorithmBayesianFIM.Rmd" )
