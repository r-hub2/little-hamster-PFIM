# Population FO variance FIM: simple (flat arm) vs covariatexoccasion mixture.
# Combination kernel lives in computePopFimCombo_Rcpp; this file is layout,
# lambda-dropping, and M = Sum_c pi_c M_c.

#' Chain-rule scale df/deta|₀ per parameter (uses \code{adjustGradient()}).
#'
#' FO pop FIM works on the random-effect (\eqn{\eta}) scale: LogNormal needs
#' \eqn{\times\mu}, Normal is \eqn{1}. Applied to mu-gradient rows only; beta
#' rows stay unscaled.
#' @noRd
#' @keywords internal
.pfimPopMuChainFactors = function( parameters ) {
  vapply( parameters, function( p ) {
    d = prop( p, "distribution" )
    adjustGradient( d, 1, prop( d, "mu" ) )
  }, numeric( 1L ) )
}

#' Outer product of a numeric vector with itself (\eqn{v v^\top}).
#' @noRd
#' @keywords internal
.outerColPopFim = function( v ) tcrossprod( v )

#' Split mu-gradient columns by occasion widths.
#'
#' Multi-occasion IOV: observation columns are concatenated per occasion;
#' \code{occ_col_widths} restores that layout so \eqn{\gamma_i^2} outer products
#' stay block-diagonal across occasions (not pooled into one dense V_IOV).
#' @noRd
#' @keywords internal
.splitGradMuByOccasion = function( G_adj_mu, occ_col_widths ) {
  if ( length( occ_col_widths ) < 2L ) return( list( G_adj_mu ) )
  offsets = c( 0L, cumsum( occ_col_widths ) )
  map( seq_along( occ_col_widths ), function( k ) {
    cols = ( offsets[[ k ]] + 1L ):offsets[[ k + 1L ]]
    G_adj_mu[ , cols, drop = FALSE ]
  } )
}

#' Block-diagonal IOV outer products across occasions for one gamma parameter.
#'
#' For parameter index \code{gamma_idx}, each occasion contributes
#' \eqn{g_{i,k} g_{i,k}^\top}; occasions are independent under FO IOV so the
#' sum is \code{bdiag}, not a full Kronecker coupling.
#' @noRd
#' @keywords internal
.iovOuterBlocks = function( gamma_idx, occ_grads_T ) {
  blocks = map( occ_grads_T, ~ .outerColPopFim( .x[ , gamma_idx, drop = FALSE ] ) )
  if ( length( blocks ) == 1L ) blocks[[ 1L ]] else as.matrix( bdiag( blocks ) )
}

#' One covariate x occasion block (MFbeta + MFVar); R reference for C++ parity tests.
#' @param mu_values Per-parameter chain-rule factors (LogNormal: \eqn{\mu};
#'   Normal: \eqn{1}); from \code{.pfimPopMuChainFactors()}.
#' @noRd
#' @keywords internal
.computePopFimCombo_R = function( gradients,
                                   mu_values,
                                   omega_iiv,
                                   gamma_values,
                                   error_variance,
                                   occ_col_widths,
                                   sigma_derivatives,
                                   has_iov ) {
  n_omega     = length( mu_values )
  n_rows      = nrow( gradients )
  n_occasions = length( occ_col_widths )
  n_sigma     = length( sigma_derivatives )

  # Rows: mu block first, then optional beta rows from covariate effects.
  gradients_mu = gradients[ seq_len( n_omega ), , drop = FALSE ]
  if ( n_rows > n_omega ) {
    gradients_beta = gradients[ ( n_omega + 1L ):n_rows, , drop = FALSE ]
  } else {
    gradients_beta = matrix( 0, nrow = 0L, ncol = ncol( gradients ) )
  }

  # Random-effect cov on eta: multi-occ expands gamma into occasion blocks;
  # single-occ pools omega²+gamma² on the diagonal; no-IOV is diag(omega²) only.
  n_omega_rows = length( omega_iiv )
  OMEGA = if ( has_iov && n_occasions > 1L ) {
    diag( c( omega_iiv, rep( gamma_values^2, n_occasions ) ) )
  } else if ( has_iov ) {
    diag( omega_iiv + gamma_values^2 )
  } else {
    diag( omega_iiv, nrow = n_omega_rows )
  }

  # eta-scale for IIV/IOV terms; MFbeta still uses raw G (mu + beta) vs V⁻¹.
  gradients_adjusted_mu = gradients_mu * mu_values
  gradients_adjusted  = if ( nrow( gradients_beta ) > 0L )
    rbind( gradients_adjusted_mu, gradients_beta )
  else
    gradients_adjusted_mu

  if ( has_iov && n_occasions > 1L ) {
    # Multi-occ FO: V = G_etaᵀ diag(omega²) G_eta + Sum_i gamma_i² bdiag_k(g_{i,k} g_{i,k}ᵀ) + R.
    # (Cannot use a single Gᵀ Ω G with Ω = diag(omega, gamma⊗I_occ) without occasion splits.)
    scaled_mu     = gradients_adjusted_mu * sqrt( omega_iiv )
    V_iiv         = as.matrix( crossprod( scaled_mu ) )
    occ_grads_T   = map( .splitGradMuByOccasion( gradients_adjusted_mu, occ_col_widths ), t )
    gamma_indices = which( gamma_values > 0 )

    V_iov = if ( length( gamma_indices ) > 0L ) {
      reduce( map( gamma_indices, function( i )
        gamma_values[ i ]^2 * .iovOuterBlocks( i, occ_grads_T )
      ), `+` )
    } else {
      matrix( 0, nrow = nrow( error_variance ), ncol = ncol( error_variance ) )
    }

    V     = as.matrix( V_iiv + V_iov + error_variance )
    V_inv = .safeCholInv( V )

    # MFbeta = G V⁻¹ Gᵀ; G includes unscaled beta rows when present.
    MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )
    grad_mu_T   = t( gradients_adjusted_mu )

    # MFVar needs dV/dgamma_i (= occasion outer blocks) then dV/dσ_k.
    dV_gamma = map( gamma_indices, function( i )
      .iovOuterBlocks( i, occ_grads_T )
    )
    dV_full  = c( dV_gamma, sigma_derivatives )
  } else {
    # One occasion (or no IOV): classic V = G_etaᵀ Ω G_eta + R.
    G_mu  = gradients_adjusted[ seq_len( n_omega_rows ), , drop = FALSE ]
    tmp   = OMEGA %*% G_mu
    V     = as.matrix( crossprod( tmp, G_mu ) + error_variance )
    V_inv = .safeCholInv( V )

    MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )
    grad_mu_T   = t( gradients_adjusted_mu )

    gamma_indices = if ( has_iov ) which( gamma_values > 0 ) else integer( 0L )
    dV_full       = sigma_derivatives
  }

  # W = eta-scaled mu grads (and pooled gamma cols when n_occ ≤ 1); dV_full carries σ / multi-occ gamma.
  W = grad_mu_T[ , seq_len( n_omega ), drop = FALSE ]
  if ( has_iov && n_occasions <= 1L && length( gamma_indices ) > 0L )
    W = cbind( W, grad_mu_T[ , gamma_indices, drop = FALSE ] )

  MFVar = computeMFVar_mixed_Rcpp( V_inv, W, map( dV_full, as.matrix ) )

  # Combination FIM is block-diagonal: FE first, then variance (omega / gamma / σ).
  n_fixed = nrow( MFbeta_full )
  n_var   = nrow( MFVar )
  out     = matrix( 0, n_fixed + n_var, n_fixed + n_var )
  out[ seq_len( n_fixed ), seq_len( n_fixed ) ] = MFbeta_full
  out[ n_fixed + seq_len( n_var ), n_fixed + seq_len( n_var ) ] = MFVar
  out
}

#' Standard population path: fixed mu, diagonal IIV, no occasion structure.
#'
#' Builds \eqn{V = G_\eta^\top \Omega G_\eta + R}, then the fixed-effect block
#' \eqn{G V^{-1} G^\top} and the mixed variance block via Rcpp.
#' @param model A \code{Model} object.
#' @param arm An \code{Arm} with flat gradients and variance.
#' @param isFixedMu Logical vector of fixed mu flags.
#' @param isFixedOmega Logical vector of fixed omega flags (unused here; MFVar
#'   already drops fixed sigmas upstream).
#' @return List with \code{MFbeta} and \code{MFVar} matrices.
#' @noRd
#' @keywords internal
.evaluateVarianceFIMPopSimple = function( model, arm, isFixedMu, isFixedOmega ) {
  parameters     = prop( model, "modelParameters" )
  parameterNames = map_chr( parameters, ~ prop( .x, "name" ) )

  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  outputNames      = prop( model, "outputNames" )

  omega_IIV     = .modelOmegaIIVVariance( model )
  muChain       = .pfimPopMuChainFactors( parameters )
  nOmega        = length( parameterNames )
  errorVariance = as.matrix( varianceResults$errorVariance )

  # Stack outputs: rows = parameters, cols = observations (combo kernel layout).
  gradients = if ( length( outputNames ) == 1L ) {
    t( allGradientsData[[ outputNames[[ 1L ]] ]] )
  } else {
    do.call( cbind, lapply( outputNames, function( nm )
      t( allGradientsData[[ nm ]] ) ) )
  }

  # Simple path = one occasion, no IOV - same Armadillo kernel as the cov path.
  fisher = as.matrix( computePopFimCombo_Rcpp(
    gradients,
    unname( muChain ),
    unname( omega_IIV ),
    rep( 0, nOmega ),
    errorVariance,
    as.integer( ncol( gradients ) ),
    lapply( varianceResults$sigmaDerivatives, as.matrix ),
    has_iov = FALSE
  ) )

  # Kernel returns full mu block; drop fixed-mu here. Lambda drop is deferred to
  # PopulationFim::evaluateFim (simple path) so sigma columns stay until then.
  mu_idx_keep = seq_len( nOmega )[ !isFixedMu ]
  var_idx     = if ( nrow( fisher ) > nOmega )
    ( nOmega + 1L ):nrow( fisher ) else integer( 0L )
  list(
    MFbeta = fisher[ mu_idx_keep, mu_idx_keep, drop = FALSE ],
    MFVar  = if ( length( var_idx ) )
      fisher[ var_idx, var_idx, drop = FALSE ]
    else
      matrix( numeric( 0 ), 0L, 0L )
  )
}

#' Covariate / occasion path: sum Fisher blocks over covariate combinations.
#'
#' Pop FIM averages at the FIM level: \eqn{M = \sum_c \pi_c M_c} (not a single
#' expected gradient). Individual/Bayesian use the same mixture contract via
#' \code{.evaluateIndBayesMixtureFim()}.
#' @noRd
#' @keywords internal
.evaluateVarianceFIMPopCovariateOccasion = function( model, arm,
                                                       isFixedMu,
                                                       isFixedOmega,
                                                       isFixedLambda = NULL ) {
  parameters     = prop( model, "modelParameters" )
  parameterNames = map_chr( parameters, ~ prop( .x, "name" ) )

  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  outputNames      = prop( model, "outputNames" )

  omega_IIV     = .modelOmegaIIVVariance( model )
  gamma_values  = vapply( parameters, function(x) pluck( x, "gamma", .default = 0 ), numeric( 1L ), USE.NAMES = FALSE )
  has_IOV       = any( gamma_values > 0 )
  muChain       = .pfimPopMuChainFactors( parameters )

  nOmega = length( parameterNames )
  nGamma = sum( gamma_values > 0 )
  nSigma = length(
    pluck( varianceResults, 1L, "variances", 1L, "variance", "sigmaDerivatives" )
  )
  numberOfOccasions = pluck( varianceResults, 1L, "variances" ) |>
    map_chr( "occasion" ) |> unique() |> length()
  # IOV gamma block only when at least one gamma > 0
  nGamma_eff = if ( has_IOV ) nGamma else 0L
  nLambda    = nOmega + nGamma_eff + nSigma

  compute_fisher_one_combination = function( iter ) {
    # Per combo: occasion grads/R are concatenated (bdiag R); kernel gets occ widths.
    gradients_by_occasion = map( seq_len( numberOfOccasions ), function( occ ) {
      gpo = map( outputNames, ~ t( pluck( allGradientsData, iter, "gradients", occ, "gradient", .x ) ) )
      if ( length( gpo ) == 1L ) gpo[[ 1L ]] else do.call( cbind, gpo )
    })

    variance_by_occasion = map(
      seq_len( numberOfOccasions ),
      ~ as.matrix( pluck( varianceResults, iter, "variances", .x, "variance", "errorVariance" ) )
    )

    # dR/dσ_k must match V's occasion block-diag layout.
    sigma_derivatives = map( seq_len( nSigma ), function( i ) {
      sdo = map( seq_len( numberOfOccasions ),
                 ~ pluck( varianceResults, iter, "variances", .x, "variance", "sigmaDerivatives", i ) )
      if ( length( sdo ) == 1L ) as.matrix( sdo[[ 1L ]] ) else as.matrix( bdiag( sdo ) )
    })

    error_variance = if ( length( variance_by_occasion ) == 1L )
      variance_by_occasion[[ 1L ]]
    else
      as.matrix( bdiag( variance_by_occasion ) )

    block = computePopFimCombo_Rcpp(
      if ( length( gradients_by_occasion ) == 1L ) gradients_by_occasion[[ 1L ]]
      else reduce( gradients_by_occasion, cbind ),
      unname( muChain ),
      unname( omega_IIV ),
      unname( gamma_values ),
      error_variance,
      vapply( gradients_by_occasion, ncol, integer( 1L ) ),
      sigma_derivatives,
      has_iov = has_IOV
    )
    as.matrix( block ) * allGradientsData[[ iter ]]$proportion
  }

  # Mixture of FIMs (not of gradients): M = Sum_c pi_c M_c.
  # n_combo is the covariate Cartesian product (small); reduce is the contract.
  fisherMatrix = reduce(
    map( seq_along( allGradientsData ), compute_fisher_one_combination ),
    `+`
  )

  # Kernel layout: [mu x nOmega | beta | lambda = omega(+gamma)+σ]; drop fixed mu from FE,
  # drop omega² only when fixedOmega (not when fixedMu).
  nBeta_total = nrow( fisherMatrix ) - nLambda - nOmega
  nMuAndBeta  = nOmega + nBeta_total
  mu_idx_keep = seq_len( nOmega )[ !isFixedMu ]
  beta_idx    = if ( nBeta_total > 0L ) ( nOmega + 1L ):nMuAndBeta else integer( 0L )
  var_idx     = ( nMuAndBeta + 1L ):nrow( fisherMatrix )

  MFbeta     = fisherMatrix[ c( mu_idx_keep, beta_idx ), c( mu_idx_keep, beta_idx ), drop = FALSE ]
  MFVar_full = fisherMatrix[ var_idx, var_idx, drop = FALSE ]
  isFixedLambda = isFixedLambda %||% isFixedOmega
  if ( length( isFixedLambda ) != nOmega )
    stop(
      "isFixedLambda length (", length( isFixedLambda ),
      ") must match the number of parameters (", nOmega, ").",
      call. = FALSE
    )
  # Gamma/sigma columns are never flagged fixed via isFixedLambda (length = nOmega).
  keepLambda = c( !isFixedLambda, rep( TRUE, nGamma_eff + nSigma ) )

  list(
    MFbeta = MFbeta,
    MFVar  = MFVar_full[ keepLambda, keepLambda, drop = FALSE ]
  )
}

#' Dispatch population variance FIM between simple and covariate/occasion paths.
#'
#' Simple path: flat arm grads, no beta/IOV occasion layout. Complex path:
#' nested comboxoccasion evaluations and lambda-dropping via \code{isFixedLambda}.
#' @noRd
#' @keywords internal
.evaluateVarianceFIMPop = function( fim, model, arm,
                                    isFixedMu     = NULL,
                                    isFixedOmega  = NULL,
                                    isFixedLambda = NULL ) {
  parameters = prop( model, "modelParameters" )
  isFixedMu = isFixedMu %||% map_lgl( parameters, .paramMuFixed )
  isFixedOmega = isFixedOmega %||% map_lgl( parameters, .paramOmegaFixed )
  isFixedLambda = isFixedLambda %||% isFixedOmega

  if ( !usesCovariateOccasionStructure( model ) ) {
    .evaluateVarianceFIMPopSimple( model, arm, isFixedMu, isFixedOmega )
  } else {
    .evaluateVarianceFIMPopCovariateOccasion(
      model, arm, isFixedMu, isFixedOmega, isFixedLambda
    )
  }
}
