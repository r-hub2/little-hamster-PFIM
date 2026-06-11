# Population FIM variance blocks (Mentré et al., 1997; covariate/IOV extensions).

.outerColPopFim <- function( v ) tcrossprod( v )

#' Standard population path: fixed mu, diagonal IIV, no occasion structure.
#' @keywords internal
.evaluateVarianceFIMPopSimple <- function( model, arm, isFixedMu, isFixedOmega ) {
  parameters     = prop( model, "modelParameters" )
  parameterNames = map_chr( parameters, ~ prop( .x, "name" ) )

  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  outputNames      = prop( model, "outputNames" )

  distributions  = map( parameters, ~ prop( .x, "distribution" ) )
  omega_IIV      = vapply( distributions, function(x) prop( x, "omega" )^2, numeric( 1L ), USE.NAMES = FALSE )
  muValues       = set_names(
    vapply( distributions, function(x) prop( x, "mu" ), numeric( 1L ), USE.NAMES = FALSE ),
    parameterNames
  )
  nOmega        = length( parameterNames )
  errorVariance = as.matrix( varianceResults$errorVariance )

  gradients         = do.call( cbind, map( outputNames, ~ t( allGradientsData[[ .x ]] ) ) )
  OMEGA             = diag( omega_IIV, nrow = nOmega )
  gradientsAdjusted = gradients * muValues

  V     = t( OMEGA %*% gradientsAdjusted ) %*% gradientsAdjusted + errorVariance
  V_inv = .safeCholInv( V )

  MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )
  gradT       = t( gradientsAdjusted )
  dV_omega    = map( seq_len( nOmega ), ~ .outerColPopFim( gradT[ , .x, drop = FALSE ] ) )
  dV_dlambda  = c( dV_omega, varianceResults$sigmaDerivatives )
  MFVar       = .computeMFVar( V_inv, dV_dlambda )

  mu_idx_keep = seq_len( nOmega )[ !isFixedMu ]
  list(
    MFbeta = MFbeta_full[ mu_idx_keep, mu_idx_keep, drop = FALSE ],
    MFVar  = MFVar
  )
}

#' Covariate / occasion path: sum Fisher blocks over covariate combinations.
#' @keywords internal
.evaluateVarianceFIMPopCovariateOccasion <- function( model, arm,
                                                       isFixedMu,
                                                       isFixedOmega ) {
  parameters     = prop( model, "modelParameters" )
  parameterNames = map_chr( parameters, ~ prop( .x, "name" ) )

  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  outputNames      = prop( model, "outputNames" )

  distributions = map( parameters, ~ prop( .x, "distribution" ) )
  omega_IIV     = map_dbl( distributions, ~ prop( .x, "omega" )^2 )
  gamma_values  = vapply( parameters, function(x) pluck( x, "gamma", .default = 0 ), numeric( 1L ), USE.NAMES = FALSE )
  has_IOV       = any( gamma_values > 0 )
  muValues      = set_names( map_dbl( distributions, ~ prop( .x, "mu" ) ), parameterNames )

  nOmega = length( parameterNames )
  nGamma = sum( gamma_values > 0 )
  nSigma = length(
    varianceResults[[ 1L ]]$variances[[ 1L ]]$variance$sigmaDerivatives
  )
  numberOfOccasions = varianceResults[[ 1L ]]$variances |>
    map_chr( "occasion" ) |> unique() |> length()
  nGamma_eff = if ( has_IOV ) nGamma else 0L
  nLambda    = nOmega + nGamma_eff + nSigma

  compute_fisher_one_combination = function( iter ) {
    variance_by_occasion = map(
      seq_len( numberOfOccasions ),
      ~ varianceResults[[ iter ]]$variances[[ .x ]]$variance$errorVariance
    )
    errorVariance = bdiag( variance_by_occasion )

    gradients_by_occasion = map( seq_len( numberOfOccasions ), function( occ ) {
      gpo = map( outputNames, ~ t( allGradientsData[[ iter ]]$gradients[[ occ ]]$gradient[[ .x ]] ) )
      if ( length( gpo ) == 1L ) gpo[[ 1L ]] else do.call( cbind, gpo )
    })
    gradients  = if ( length( gradients_by_occasion ) == 1L ) gradients_by_occasion[[ 1L ]]
    else reduce( gradients_by_occasion, cbind )
    proportion = allGradientsData[[ iter ]]$proportion
    nTotalRows = nrow( gradients )

    gradients_mu   = gradients[ seq_len( nOmega ), , drop = FALSE ]
    gradients_beta = if ( nTotalRows > nOmega )
      gradients[ ( nOmega + 1L ):nTotalRows, , drop = FALSE ]
    else matrix( 0, nrow = 0L, ncol = ncol( gradients ) )

    OMEGA = if ( has_IOV && numberOfOccasions > 1L )
      diag( c( omega_IIV, rep( gamma_values^2, numberOfOccasions ) ) )
    else if ( has_IOV )
      diag( omega_IIV + gamma_values^2 )
    else
      diag( omega_IIV, nrow = length( omega_IIV ) )

    gradientsAdjusted_mu = gradients_mu * muValues
    gradientsAdjusted    = if ( nrow( gradients_beta ) > 0L )
      rbind( gradientsAdjusted_mu, gradients_beta )
    else gradientsAdjusted_mu

    if ( has_IOV && numberOfOccasions > 1L ) {
      n_occ1 = ncol( gradients_by_occasion[[ 1L ]] )
      A_mu   = gradientsAdjusted_mu[ ,              seq_len( n_occ1 ),   drop = FALSE ]
      B_mu   = gradientsAdjusted_mu[ , ( n_occ1 + 1L ):ncol( gradientsAdjusted_mu ), drop = FALSE ]

      scaled_mu = gradientsAdjusted_mu * sqrt( omega_IIV )
      V_iiv     = as.matrix( crossprod( scaled_mu ) )

      A_mu_T = t( A_mu )
      B_mu_T = t( B_mu )

      gamma_indices = which( gamma_values > 0 )
      V_iov = if ( length( gamma_indices ) > 0L )
        reduce(
          map( gamma_indices, ~ gamma_values[ .x ]^2 * as.matrix( bdiag(
            .outerColPopFim( A_mu_T[ , .x, drop = FALSE ] ),
            .outerColPopFim( B_mu_T[ , .x, drop = FALSE ] )
          ))),
          `+`
        )
      else matrix( 0, nrow( errorVariance ), ncol( errorVariance ) )

      V     = V_iiv + V_iov + errorVariance
      V_inv = .safeCholInv( V )

      MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )

      gradMu_T = t( gradientsAdjusted_mu )
      dV_omega = map( seq_len( nOmega ),
                      ~ .outerColPopFim( gradMu_T[ , .x, drop = FALSE ] ) )
      dV_gamma = map( gamma_indices, ~ as.matrix( bdiag(
        .outerColPopFim( A_mu_T[ , .x, drop = FALSE ] ),
        .outerColPopFim( B_mu_T[ , .x, drop = FALSE ] )
      )))

    } else {
      nOmegaRows = nrow( OMEGA )
      tmp        = OMEGA %*% gradientsAdjusted[ seq_len( nOmegaRows ), , drop = FALSE ]
      V          = t( tmp ) %*% gradientsAdjusted[ seq_len( nOmegaRows ), , drop = FALSE ] +
        errorVariance
      V_inv      = .safeCholInv( V )

      MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )

      gradMu_T = t( gradientsAdjusted_mu )
      dV_omega = map( seq_len( nOmega ),
                      ~ .outerColPopFim( gradMu_T[ , .x, drop = FALSE ] ) )
      dV_gamma = if ( has_IOV )
        map( which( gamma_values > 0 ),
             ~ .outerColPopFim( gradMu_T[ , .x, drop = FALSE ] ) )
      else list()
    }

    dV_sigma = map( seq_len( nSigma ), function( i ) {
      sdo = map( seq_len( numberOfOccasions ),
                 ~ varianceResults[[ iter ]]$variances[[ .x ]]$variance$sigmaDerivatives[[ i ]] )
      if ( length( sdo ) == 1L ) as.matrix( sdo[[ 1L ]] )
      else as.matrix( bdiag( sdo ) )
    })

    dV_dlambda = c( dV_omega, dV_gamma, dV_sigma )
    MFVar = .computeMFVar( V_inv, dV_dlambda )

    as.matrix( bdiag( MFbeta_full, MFVar ) ) * proportion
  }

  fisherMatrix = Reduce( "+", lapply( seq_along( allGradientsData ), compute_fisher_one_combination ) )

  nBeta_total = nrow( fisherMatrix ) - nLambda - nOmega
  nMuAndBeta  = nOmega + nBeta_total
  mu_idx_keep = seq_len( nOmega )[ !isFixedMu ]
  beta_idx    = if ( nBeta_total > 0L ) ( nOmega + 1L ):nMuAndBeta else integer( 0L )
  var_idx     = ( nMuAndBeta + 1L ):nrow( fisherMatrix )

  MFbeta     = fisherMatrix[ c( mu_idx_keep, beta_idx ), c( mu_idx_keep, beta_idx ), drop = FALSE ]
  MFVar_full = fisherMatrix[ var_idx, var_idx, drop = FALSE ]
  keepLambda = c( !isFixedOmega, rep( TRUE, nGamma_eff + nSigma ) )

  list(
    MFbeta = MFbeta,
    MFVar  = MFVar_full[ keepLambda, keepLambda, drop = FALSE ]
  )
}

#' Dispatch population variance FIM between simple and covariate/occasion paths.
#' @keywords internal
.evaluateVarianceFIMPop <- function( fim, model, arm,
                                    isFixedMu    = NULL,
                                    isFixedOmega = NULL ) {
  parameters = prop( model, "modelParameters" )
  isFixedMu = isFixedMu %||% map_lgl( parameters,
                                      ~ isTRUE( prop( .x, "fixedMu" ) ) ||
                                        prop( prop( .x, "distribution" ), "mu" ) == 0 )
  isFixedOmega = isFixedOmega %||% map_lgl( parameters,
                                            ~ isTRUE( prop( .x, "fixedOmega" ) ) ||
                                              prop( prop( .x, "distribution" ), "omega" ) == 0 )

  if ( !usesCovariateOccasionStructure( model ) ) {
    .evaluateVarianceFIMPopSimple( model, arm, isFixedMu, isFixedOmega )
  } else {
    .evaluateVarianceFIMPopCovariateOccasion( model, arm, isFixedMu, isFixedOmega )
  }
}
