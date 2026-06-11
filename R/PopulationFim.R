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
#' @include Fim.R
#' @export

PopulationFim = new_class( "PopulationFim", package = "PFIM", parent = Fim )
S4_register( PopulationFim )

#' Compute the FIM for one arm
#' @name evaluateFim
#' @export
method( evaluateFim, list( PopulationFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameters          = prop( model, "modelParameters" )
  armSize             = prop( arm,   "size"            )
  hasComplexStructure = usesCovariateOccasionStructure( model )

  isFixedMu    = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedMu"    ) ) ||
                            prop( prop( .x, "distribution" ), "mu"    ) == 0 )
  isFixedOmega = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedOmega" ) ) ||
                            prop( prop( .x, "distribution" ), "omega" ) == 0 )

  result = .evaluateVarianceFIMPop( fim, model, arm,
                                    isFixedMu    = isFixedMu,
                                    isFixedOmega = isFixedOmega )

  if ( !hasComplexStructure ) {
    nSigma     = ncol( result$MFVar ) - length( isFixedOmega )
    keepLambda = c( !isFixedOmega, rep( TRUE, nSigma ) )
    MFVar      = result$MFVar[ keepLambda, keepLambda, drop = FALSE ]
  } else {
    MFVar = result$MFVar
  }

  prop( fim, "fisherMatrix" ) = as.matrix(
    bdiag( result$MFbeta * armSize, MFVar * armSize )
  )
  fim
}

.muNames = function( parameters, greek )
  parameters |>
  keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) &&
          prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
  map_chr( ~ paste0( greek, prop( .x, "name" ) ) )

.omegaNames = function( parameters, greek )
  parameters |>
  keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) &&
          prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
  map_chr( ~ paste0( greek, prop( .x, "name" ) ) )

.gammaNames = function( parameters, greek )
  parameters |>
  keep( ~ pluck( .x, "gamma", .default = 0 ) > 0 &&
          !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
  map_chr( ~ paste0( greek, prop( .x, "name" ) ) )

.sigmaNames = function( modelError, greekSigma, sep = "_" )
  modelError |>
  map( function( err ) {
    out = prop( err, "output" )
    c(
      if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
        paste0( greekSigma, sep, "inter_", out ),
      if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
        paste0( greekSigma, sep, "slope_", out )
    )
  }) |> unlist( use.names = FALSE )

#' Store optimal arm designs on the FIM object
#' @name setOptimalArms
#' @export
method( setOptimalArms, list( PopulationFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {

    outputs      = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
    weightsIndex = outputs$weightsIndex
    weights      = outputs$optimalWeights
    N_total      = prop( outputs$armFims[[ weightsIndex[ 1L ] ]][[ 1L ]], "size" )
    nPerGroup    = round( N_total * weights / sum( weights ) )

    armList = map2( weightsIndex, nPerGroup, function( idx, size ) {
      arm = outputs$armFims[[ idx ]][[ 1L ]]
      prop( arm, "size" ) = size
      prop( arm, "name" ) = paste0( "Arm", idx )
      arm
    })
    armList[ order( map_dbl( armList, ~ prop( .x, "size" ) ), decreasing = TRUE ) ]
  }

#' Store optimal arm designs on the FIM object
#' @name setOptimalArms
#' @export
method( setOptimalArms, list( PopulationFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {
    outputs             = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
    numberOfIndividuals = outputs$numberOfIndividuals
    map2( outputs$listArms, seq_along( outputs$listArms ), function( la, i ) {
      prop( la$arm, "size" ) = as.double( numberOfIndividuals[[ i ]] )
      prop( la$arm, "name" ) = paste0( "Arm", i )
      la
    })
  }

.plotFimBars = function( fim, evaluation, metric ) {
  parameters = prop( evaluation, "modelParameters" )
  modelError = prop( evaluation, "modelError"       )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se         = prop( fim, "SEAndRSE" )
  greek      = .greekPlot
  has_IOV    = any( map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) ) > 0 )

  paramsMu    = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) &&
            prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) )
  paramsOmega = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) &&
            prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) )
  paramsGamma = if ( has_IOV )
    parameters |>
    keep( ~ pluck( .x, "gamma", .default = 0 ) > 0 &&
            !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
    map_chr( ~ prop( .x, "name" ) )
  else character( 0L )
  paramsSigma = .sigmaNames( modelError, greek[ "sigma" ] )

  y_vals = if ( metric == "SE" ) se$SE$SE else se$RSE$RSE
  cats   = paste0( metric, " ", c(
    rep( greek[ "mu"    ], length( paramsMu    ) ),
    rep( greek[ "omega" ], length( paramsOmega ) ),
    rep( greek[ "gamma" ], length( paramsGamma ) ),
    rep( greek[ "sigma" ], length( paramsSigma ) )
  ))

  df = data.frame( Parameter        = c( paramsMu, paramsOmega, paramsGamma, paramsSigma ),
                   parametersValues = se$SEAndRSE$parametersValues,
                   y                = y_vals,
                   cat              = cats )
  names( df )[ 3L ] = metric

  facet_levels = paste0( metric, " ", c( greek[ "mu" ], greek[ "omega" ],
                                         if ( has_IOV ) greek[ "gamma" ],
                                         greek[ "sigma" ] ) )

  ggplot( df, aes( x = .data[["Parameter"]], y = .data[[ metric ]] ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~ factor( cat, levels = facet_levels ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.title  = element_text( size = 16 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ),
           axis.text.y = element_text( size = 16 ),
           strip.text.x = element_text( size = 16 ) )
}

#' Attach evaluated FIM results to a project
#' @name setEvaluationFim
#' @export
method( setEvaluationFim, PopulationFim ) = function( fim, evaluation ) {

  parameters = prop( evaluation, "modelParameters" )
  modelError = prop( evaluation, "modelError"       )
  greek      = .greekConsole

  has_IOV = any( map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) ) > 0 )
  fe      = .fimFixedEffectLabels( evaluation, greek )

  columnNamesMu   = .muNames( parameters, greek[ "mu" ] )
  columnNamesBeta = fe$columnNamesBeta
  muValues        = fe$muValues
  betaValues      = fe$betaValues

  columnNamesOmega = .omegaNames( parameters, greek[ "omega" ] )
  columnNamesGamma = if ( has_IOV ) .gammaNames( parameters, greek[ "gamma" ] ) else character( 0L )
  columnNamesSigma = .sigmaNames( modelError, greek[ "sigma" ] )

  omegaValues = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) &&
            prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_dbl( ~ prop( prop( .x, "distribution" ), "omega" )^2 )

  gammaValues = if ( has_IOV )
    parameters |>
    keep( ~ pluck( .x, "gamma", .default = 0 ) > 0 &&
            !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
    map_dbl( ~ pluck( .x, "gamma" )^2 )
  else numeric( 0L )

  sigmaValues = modelError |>
    map( function( err ) {
      v = list()
      if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
        v$sigmaInter = prop( err, "sigmaInter" )
      if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
        v$sigmaSlope = prop( err, "sigmaSlope" )
      v
    }) |> unlist( use.names = FALSE )

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
  SE    = sqrt( diag( .safeCholInv( M ) ) )
  RSE   = SE / abs( pVals ) * 100
  seDF  = data.frame( parametersValues = pVals, SE = SE, RSE = RSE,
                      row.names = allNames )

  prop( fim, "fisherMatrix"              ) = M
  prop( fim, "fixedEffects"              ) = fixedEffects
  prop( fim, "varianceEffects"           ) = varEffects
  prop( fim, "condNumberFixedEffects"    ) = .conditionNumber( fixedEffects )
  prop( fim, "condNumberVarianceEffects" ) = .conditionNumber( varEffects   )
  prop( fim, "SEAndRSE" ) = list(
    SE       = seDF[ , c( "parametersValues", "SE"  ) ],
    RSE      = seDF[ , c( "parametersValues", "RSE" ) ],
    table    = seDF,
    SEAndRSE = seDF
  )
  fim
}

.showPopulationFimConsole = function( fim, evaluation = NULL ) {

  .hdr = function( t ) { cat( "\n*************************************** \n" )
    cat( " ", t, "\n" )
    cat( "*************************************** \n\n" ) }

  .hdr( "Parameters estimation" )
  .printFimSeAndRse( fim, evaluation )

  cat( "\n********************************************* \n",
       " determinant, condition numbers and d-criterion \n",
       "*********************************************** \n\n" )
  cat( "determinant:",  det( prop( fim, "fisherMatrix" ) ), "\n" )
  cat( "D-criterion:",  Dcriterion( fim ),                  "\n" )
  cat( "Conditional number (fixed effects):",      prop( fim, "condNumberFixedEffects"    ), "\n" )
  cat( "Conditional number (variance components):", prop( fim, "condNumberVarianceEffects"), "\n" )

  .hdr( "Population Fisher Matrix" );                     print( prop( fim, "fisherMatrix"   ) )
  .hdr( "Fixed effects (\u03bc)" );                       print( prop( fim, "fixedEffects"   ) )
  .hdr( "Variance components (\u03c9\u00B2, \u03b3\u00B2, \u03c3\u00B2)" )
  print( prop( fim, "varianceEffects" ) )

  if ( any( grepl( "\u03b3\u00B2", rownames( prop( fim, "fisherMatrix" ) ) ) ) ) {
    cat( "\n*************************************** \n Legend: \n" )
    cat( " \u03bc  = fixed effects (population means)\n" )
    cat( " \u03c9\u00B2 = inter-individual variability (IIV)\n" )
    cat( " \u03b3\u00B2 = inter-occasion variability (IOV)\n" )
    cat( " \u03c3\u00B2 = residual error variance\n" )
    cat( "*************************************** \n\n" )
  }

  invisible( fim )
}

#' Print FIM summaries to the console
#' @name showFIM
#' @export
method( showFIM, PopulationFim ) = function( fim ) {
  .showPopulationFimConsole( fim, evaluation = NULL )
}

#' Barplot of standard errors from a FIM
#' @name plotSEFIM
#' @export
method( plotSEFIM,  list( PopulationFim, PFIMProject ) ) =
  function( fim, evaluation ) .plotFimBars( fim, evaluation, "SE"  )

#' Barplot of relative standard errors from a FIM
#' @name plotRSEFIM
#' @export
method( plotRSEFIM, list( PopulationFim, PFIMProject ) ) =
  function( fim, evaluation ) .plotFimBars( fim, evaluation, "RSE" )

#' FIM tables for HTML reports
#' @name tablesForReport
#' @export
method( tablesForReport, list( PopulationFim, PFIMProject ) ) = function( fim, evaluation ) {

  SEAndRSE        = prop( fim, "SEAndRSE"     )$SEAndRSE
  M               = prop( fim, "fisherMatrix" )
  fe              = as.matrix( prop( fim, "fixedEffects"   ) )
  ve              = prop( fim,              "varianceEffects" )
  cn1             = prop( fim, "condNumberFixedEffects"    )
  cn2             = prop( fim, "condNumberVarianceEffects" )
  parameters      = prop( evaluation, "modelParameters" )
  modelError      = prop( evaluation, "modelError"       )
  modelCovariates = prop( evaluation, "modelCovariates"  )
  greek           = .greekLatex

  has_IOV = any( map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) ) > 0 )
  has_cov = length( modelCovariates ) > 0L

  muL    = paste0( .muNames( parameters, greek[ "mu" ] ),    "}$" )
  betaL  = if ( has_cov )
    paste0( greek[ "beta" ],
            str_remove( .betaInternalNamesFromCovariates( modelCovariates ), "^beta_" ), "}$" )
  else character( 0L )
  omegaL = paste0( .omegaNames( parameters, greek[ "omega" ] ), "}$" )
  gammaL = if ( has_IOV ) paste0( .gammaNames( parameters, greek[ "gamma" ] ), "}$" )
  else character( 0L )
  sigmaL = modelError |>
    map( function( err ) {
      out  = prop( err, "output" )
      sig2 = "$\\sigma^2_{"
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( sig2, "inter_", out, "}$" ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( sig2, "slope_", out, "}$" )
      )
    }) |> unlist( use.names = FALSE )

  feNames = c( muL, betaL );  dimnames( fe ) = list( feNames, feNames )
  veNames = c( omegaL, gammaL, sigmaL );
  colnames( ve ) = veNames;  rownames( ve ) = veNames

  .kbl_s = function( df )
    kbl( df ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  FIMCriteriaTable = data.frame( det   = det( M ),
                                 dcrit = Dcriterion( fim ),
                                 fe    = cn1, ve = cn2 ) |>
    kbl( col.names = c( "", "", "Fixed effects", "Variance effects" ),
         align = "c" ) |>
    add_header_above( c( "determinant" = 1, "d-criterion" = 1,
                         "Condition number" = 2 ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  SEAndRSETable = data.frame(
    c( muL, betaL, omegaL, gammaL, sigmaL ), round( SEAndRSE, 3 )
  ) |>
    (\( df ) { row.names( df ) = NULL; df })() |>
    kbl( col.names = c( "Parameters", "Parameter values", "SE", "RSE (%)"),
         align = "c", row.names = FALSE ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  list( fixedEffectsTable    = .kbl_s( fe ),
        varianceEffectsTable = .kbl_s( ve ),
        FIMCriteriaTable     = FIMCriteriaTable,
        SEAndRSETable        = SEAndRSETable )
}

#' Render the evaluation HTML report
#' @name generateReportEvaluation
#' @export
method( generateReportEvaluation, PopulationFim ) =
  .renderEvalReport( "EvaluationPopulationFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, MultiplicativeAlgorithm ) ) =
  .renderReport( "OptimizationMultiplicativeAlgorithmPopulationFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, FedorovWynnAlgorithm ) ) =
  .renderReport( "OptimizationFedorovWynnAlgorithmPopulationFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, SimplexAlgorithm ) ) =
  .renderReport( "OptimizationSimplexAlgorithmPopulationFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, PSOAlgorithm ) ) =
  .renderReport( "OptimizationPSOAlgorithmPopulationFIM.Rmd" )

#' Render the optimization HTML report
#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, PGBOAlgorithm ) ) =
  .renderReport( "OptimizationPGBOAlgorithmPopulationFIM.Rmd" )
