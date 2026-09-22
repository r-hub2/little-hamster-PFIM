# Subject-level FIM helpers (individual / Bayesian).
#
# Population FIM adds information across subjects. Individual and Bayesian FIMs
# describe one subject's precision: when strata (arms x covariate combinations)
# have different protocols, average the covariances then invert:
#   C_bar = sum_s w_s M_s^{-1},   M_eff = C_bar^{-1}.
# Arithmetic mixing of FIMs underestimates shrinkage (Jensen).

#' Subject-level covariance from a (possibly singular) FIM.
#'
#' Correlation-scale eigen-decomposition: identifiable directions get the
#' corresponding pseudoinverse; directions in the numerical null space mark
#' touched parameters with \code{Inf} variance and zero cross-covariances
#' (same contract as \code{.fimPinvCovarianceDiagonal}).
#' @noRd
#' @keywords internal
.pfimSubjectCovariance = function( M ) {
  .pfimPsdPseudoInverse( M )
}

#' Invert a covariance that may contain \code{Inf} on non-identifiable parameters.
#'
#' Only the finite principal block is inverted; other rows/columns stay 0 in the
#' information matrix (no information).
#' @noRd
#' @keywords internal
.pfimCovarianceToFim = function( C ) {
  C = as.matrix( C )
  p = nrow( C )
  if ( !p )
    return( matrix( numeric( 0L ), 0L, 0L ) )
  finite = which( is.finite( diag( C ) ) )
  M = matrix( 0, p, p )
  if ( !length( finite ) )
    return( M )
  Cf = C[ finite, finite, drop = FALSE ]
  Cf[ !is.finite( Cf ) ] = 0
  Mf = tryCatch( .safeCholInv( Cf ), error = function( e ) NULL )
  if ( is.null( Mf ) ) {
    # Singular finite block: correlation-scale pseudoinverse; Inf -> 0 information.
    Mf = .pfimPsdPseudoInverse( Cf )
    Mf[ !is.finite( Mf ) ] = 0
  }
  M[ finite, finite ] = Mf
  0.5 * ( M + t( M ) )
}

#' Harmonic-mean (covariance-mixture) design FIM from subject-level blocks.
#'
#' @param fisherBlocks List of square subject FIMs (same dimension).
#' @param weights Non-negative weights summing to 1 (stratum proportions).
#' @return List with \code{fisherMatrix} (= \eqn{\bar C^{-1}}) and \code{covariance}
#'   (= \eqn{\bar C}).
#' @noRd
#' @keywords internal
.pfimHarmonicMeanFim = function( fisherBlocks, weights ) {
  n = length( fisherBlocks )
  if ( !n )
    return( list(
      fisherMatrix = matrix( numeric( 0L ), 0L, 0L ),
      covariance   = matrix( numeric( 0L ), 0L, 0L )
    ) )
  weights = as.numeric( weights )
  if ( length( weights ) != n )
    .pfimStop( "subject-level FIM weights do not match the number of protocols." )
  if ( any( !is.finite( weights ) ) || any( weights < 0 ) )
    .pfimStop( "subject-level FIM weights must be finite and non-negative." )
  wsum = sum( weights )
  if ( !( wsum > 0 ) )
    .pfimStop( "subject-level FIM weights must sum to a positive value." )
  weights = weights / wsum

  # Inf-aware mix: a positive weight on Inf variance must stay Inf. A naive
  # Reduce(`+`) turns Inf-Inf / 0*Inf into NaN; keep the sequential contract.
  p = nrow( fisherBlocks[[ 1L ]] )
  C_bar = matrix( 0, p, p )
  for ( i in seq_len( n ) ) {
    Mi = as.matrix( fisherBlocks[[ i ]] )
    if ( nrow( Mi ) != p || ncol( Mi ) != p )
      .pfimStop( "subject-level FIM blocks must have the same dimensions." )
    if ( !( weights[[ i ]] > 0 ) )
      next
    Ci = .pfimSubjectCovariance( Mi )
    # Inf variance dominates: any positive weight on Inf stays Inf.
    C_bar = C_bar + weights[[ i ]] * Ci
  }
  # Clean NaN from Inf - Inf or 0 * Inf edge cases.
  C_bar[ is.nan( C_bar ) ] = Inf
  C_bar = 0.5 * ( C_bar + t( C_bar ) )
  # Parameters that are Inf in any stratum stay Inf (null cross-cov).
  inf_idx = which( !is.finite( diag( C_bar ) ) )
  C_bar = .pfimZeroCrossInfDiag( C_bar, inf_idx )
  M_eff = .pfimCovarianceToFim( C_bar )
  list( fisherMatrix = M_eff, covariance = C_bar )
}

#' Drop covariate beta columns from fixed-effect names (Ind / Bayes subject FIM).
#' @noRd
#' @keywords internal
.pfimSubjectFeCols = function( feCols ) {
  feCols[ !startsWith( feCols, "beta_" ) ]
}

#' Bayesian FE column names for eta parameters (\eqn{\omega > 0}; beta excluded).
#'
#' Flat models use bare parameter names; covariate/IOV layouts use \code{mu_*}
#' prefixes. Unlike population FIM, \code{fixedMu}/\code{fixedOmega} do not drop
#' an eta when \eqn{\omega > 0}.
#' @noRd
#' @keywords internal
.pfimBayesianMuCols = function( model, arm = NULL ) {
  parameters = prop( model, "modelParameters" )
  eta = .pfimBayesianEtaParameters( parameters )
  if ( !length( eta ) )
    return( character( 0L ) )
  etaBare = map_chr( eta, ~ prop( .x, "name" ) )

  if ( !is.null( arm ) && usesCovariateOccasionStructure( model ) ) {
    cn = colnames( getArmEvaluationGradientsMatrix( model, arm, evalModel = model ) )
    if ( is.null( cn ) ) cn = character( 0L )
    wanted = paste0( "mu_", etaBare )
    hit = wanted[ wanted %in% cn ]
    if ( length( hit ) ) return( hit )
    # Fallback: bare names present in the gradient layout.
    hit = etaBare[ etaBare %in% cn ]
    if ( length( hit ) ) return( hit )
    return( wanted )
  }

  # Flat path: gradient columns are bare parameter names.
  etaBare
}

#' TRUE when FIM is individual or Bayesian (subject-level precision).
#' @noRd
#' @keywords internal
.pfimIsSubjectLevelFim = function( fim ) {
  S7::S7_inherits( fim, IndividualFim ) || S7::S7_inherits( fim, BayesianFim )
}

#' Bayesian subject FIM on the random-effect scale.
#'
#' Data block \eqn{G^\top V^{-1} G} for every \eqn{\omega > 0} parameter, then
#' \eqn{M^\top M_{\mathrm{pk}} M + \Omega^{-1}} with \eqn{M=\mathrm{diag}(\mu)}
#' for LogNormal. Fix flags do not remove eta.
#' @noRd
#' @keywords internal
.pfimBayesianSubjectFim = function( M_data, model ) {
  parameters = .pfimBayesianEtaParameters( prop( model, "modelParameters" ) )
  n = length( parameters )
  if ( !n )
    return( matrix( numeric( 0L ), 0L, 0L ) )
  if ( nrow( M_data ) != n || ncol( M_data ) != n )
    .pfimInternalStop(
      "Bayesian data FIM dimension does not match the number of IIV parameters."
    )
  mu = vapply( parameters, function( p ) {
    d = prop( p, "distribution" )
    adjustGradient( d, 1, prop( d, "mu" ) )
  }, numeric( 1L ) )
  omega = vapply( parameters, function( p ) prop( prop( p, "distribution" ), "omega" ), numeric( 1L ) )
  mu_d  = if ( n == 1L ) mu[[ 1L ]] else diag( mu )
  as.matrix( t( mu_d ) %*% M_data %*% mu_d + .bayesianOmegaInv( omega ) )
}

#' Bayesian SE / RSE from an eta-scale FIM.
#'
#' For LogNormal, the FIM is on \eqn{\eta}; displayed SE on \eqn{\theta} is
#' \eqn{\mu\cdot\mathrm{SE}_\eta} and RSE is \eqn{100\cdot\mathrm{SE}_\eta}.
#' For Normal, SE stays on the natural scale with RSE \eqn{100\cdot\mathrm{SE}/|\mu|}.
#' @noRd
#' @keywords internal
.pfimBayesianSeRse = function( fisherMatrix, parameters, allNames ) {
  M = as.matrix( fisherMatrix )
  p = nrow( M )
  if ( p != length( allNames ) )
    .pfimInternalStop( "Bayesian SE labels do not match the FIM dimension." )

  cov_mat = tryCatch( .safeCholInv( M ), error = function( e ) NULL )
  singular = is.null( cov_mat )
  if ( singular ) {
    se_eta = sqrt( .fimPinvCovarianceDiagonal( M ) )
  } else {
    se_eta = sqrt( pmax( diag( cov_mat ), 0 ) )
  }
  se_eta[ !is.finite( se_eta ) ] = Inf

  # Map greek-prefixed mu labels back to ModelParameter objects.
  greek_mu = .greekConsole[ "mu" ]
  bare = .stripGreekPrefix( allNames, greek_mu )
  paramByName = set_names( parameters, map_chr( parameters, ~ prop( .x, "name" ) ) )

  objs    = lapply( seq_len( p ), function( i ) paramByName[[ bare[[ i ]] ]] )
  missing = vapply( objs, is.null, logical( 1L ) )
  is_logn = vapply( seq_len( p ), function( i ) {
    if ( missing[[ i ]] )
      FALSE
    else
      S7::S7_inherits( prop( objs[[ i ]], "distribution" ), LogNormal )
  }, logical( 1L ) )
  mu = vapply( seq_len( p ), function( i ) {
    if ( missing[[ i ]] )
      NA_real_
    else
      prop( prop( objs[[ i ]], "distribution" ), "mu" )
  }, numeric( 1L ) )

  se_theta = se_eta
  rse      = rep( NA_real_, p )
  rse[ missing ] = ifelse( is.finite( se_eta[ missing ] ), NA_real_, Inf )

  logn = !missing & is_logn
  se_theta[ logn ] = abs( mu[ logn ] ) * se_eta[ logn ]
  rse[ logn ]      = 100 * se_eta[ logn ]

  nat = !missing & !is_logn
  se_theta[ nat ] = se_eta[ nat ]
  rse[ nat ] = ifelse( abs( mu[ nat ] ) > 0, 100 * se_eta[ nat ] / abs( mu[ nat ] ), Inf )

  seDF = data.frame(
    parametersValues = mu,
    SE  = se_theta,
    RSE = rse
  )
  rownames( seDF ) = allNames
  list(
    SE       = seDF[ , c( "parametersValues", "SE"  ), drop = FALSE ],
    RSE      = seDF[ , c( "parametersValues", "RSE" ), drop = FALSE ],
    table    = seDF,
    SEAndRSE = seDF,
    singular = singular
  )
}

#' Index of the best subject-level protocol by D-criterion.
#'
#' For individual / Bayesian designs the D-criterion is convex in the weights of
#' a covariance mixture, so the optimum is a single protocol (vertex).
#' @param fisherMatrices List of candidate subject FIMs.
#' @return List with \code{index} (1-based) and \code{Dcriterion}.
#' @noRd
#' @keywords internal
.pfimBestSubjectProtocol = function( fisherMatrices ) {
  if ( !length( fisherMatrices ) )
    .pfimInternalStop( "no candidate protocols for subject-level D-criterion." )
  dvals = vapply(
    fisherMatrices,
    function( M ) .fimDcriterionFromMatrix( as.matrix( M ) ),
    numeric( 1L )
  )
  idx = which.max( dvals )
  list( index = as.integer( idx ), Dcriterion = dvals[[ idx ]], dvals = dvals )
}
