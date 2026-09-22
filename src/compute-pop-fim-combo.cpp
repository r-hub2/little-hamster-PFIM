// Population FIM for one covariate × occasion combination.
// Must stay numerically identical to R helper .computePopFimCombo_R
// (cross-tested in testthat) — any change here needs a matching R update.
//
// Pipeline:
//   1. Scale the mu-block of the gradient by mu (LogNormal chain rule).
//   2. Build residual variance V = G Omega G' (+ IOV) + error.
//   3. Assemble dV / d(theta) as rank-1 columns W plus optional full blocks
//      (IOV gamma, residual sigma).
//   4. Return block-diagonal FIM: [ MFbeta  0 ; 0  MFVar ] with
//      MFbeta = G V^{-1} G' and MFVar = mfvar_mixed(V^{-1}, W, full).
//
// Individual / Bayesian FIMs are assembled elsewhere; this kernel does not
// add a Bayesian prior block.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <pfim/pfim-linalg.hpp>
#include <vector>

using pfim::chol_inv;
using pfim::mfvar_mixed;
using pfim::outer_col;

namespace {

// Row i of G_mu <- row i * mu_i (in-place LogNormal adjustment).
void scale_mu_rows( arma::mat& G_mu, const arma::vec& mu ) {
  for ( arma::uword i = 0; i < G_mu.n_rows; ++i )
    G_mu.row( i ) *= mu( i );
}

// Indices of parameters with positive IOV (gamma_i > 0).
std::vector<int> active_gamma( const arma::vec& gamma ) {
  std::vector<int> out;
  for ( arma::uword i = 0; i < gamma.n_elem; ++i )
    if ( gamma( i ) > 0.0 )
      out.push_back( static_cast<int>( i ) );
  return out;
}

// Per-occasion slices of the mu-adjusted gradient (multi-occasion IOV path).
struct OccasionBlocks {
  std::vector<arma::mat> T_blocks;   // occasion k: (n_obs_k x n_omega)
  std::vector<arma::uword> offsets;  // column start of occasion k in full V
};

// Split G_adj_mu columns into occasions using occ_col_widths (must sum to ncol).
OccasionBlocks split_by_occasions( const arma::mat& G_adj_mu,
                                   const Rcpp::IntegerVector& occ_col_widths ) {
  int width_sum = 0;
  for ( int k = 0; k < occ_col_widths.size(); ++k ) {
    if ( occ_col_widths[ k ] <= 0 )
      Rcpp::stop( "computePopFimCombo_Rcpp: occ_col_widths entries must be positive." );
    width_sum += occ_col_widths[ k ];
  }
  if ( width_sum != static_cast<int>( G_adj_mu.n_cols ) )
    Rcpp::stop( "computePopFimCombo_Rcpp: occ_col_widths must sum to gradient columns." );

  OccasionBlocks occ;
  arma::uword col = 0;
  for ( int k = 0; k < occ_col_widths.size(); ++k ) {
    const int w = occ_col_widths[ k ];
    occ.T_blocks.push_back( G_adj_mu.cols( col, col + w - 1 ).t() );
    occ.offsets.push_back( col );
    col += static_cast<arma::uword>( w );
  }
  return occ;
}

// Add gamma_idx contribution: block-diagonal outer products per occasion.
void add_iov_blocks( arma::mat& target, const OccasionBlocks& occ,
                     int gamma_idx, double scale ) {
  for ( size_t k = 0; k < occ.T_blocks.size(); ++k ) {
    const arma::vec col = occ.T_blocks[ k ].col( gamma_idx );
    const arma::uword off = occ.offsets[ k ];
    const arma::uword sz  = col.n_elem;
    target.submat( off, off, off + sz - 1, off + sz - 1 ) += scale * outer_col( col );
  }
}

// Full p x p dV/d(gamma_idx) for mfvar_mixed "full" list (multi-occasion IOV).
arma::mat iov_derivative_block( const OccasionBlocks& occ, int gamma_idx ) {
  const arma::uword p = occ.T_blocks.empty()
    ? 0
    : occ.offsets.back() + occ.T_blocks.back().n_rows;
  arma::mat block( p, p, arma::fill::zeros );
  add_iov_blocks( block, occ, gamma_idx, 1.0 );
  return block;
}

// Marginal variance of observations:
//   multi-occasion IOV: V = (G sqrt(Omega))^T (...) + sum_k gamma_k^2 G_k G_k' + Sigma
//   single-occasion IOV: fold gamma^2 into the IIV scale per parameter
//   no IOV:              V = (G sqrt(Omega))^T (...) + Sigma
// When multi-occasion, `occ` must already be split (shared with build_dV_parts).
arma::mat build_V( const arma::mat& G_adj_mu,
                   const arma::vec& omega_iiv,
                   const arma::vec& gamma,
                   const arma::mat& error_variance,
                   const Rcpp::IntegerVector& occ_col_widths,
                   bool has_iov,
                   const OccasionBlocks* occ ) {
  const int n_omega     = static_cast<int>( G_adj_mu.n_rows );
  const int n_occasions = occ_col_widths.size();
  const std::vector<int> g_idx = active_gamma( gamma );

  if ( has_iov && n_occasions > 1 ) {
    arma::mat scaled = G_adj_mu;
    for ( int i = 0; i < n_omega; ++i )
      scaled.row( i ) *= std::sqrt( omega_iiv( i ) );

    arma::mat V_iov( error_variance.n_rows, error_variance.n_cols, arma::fill::zeros );
    for ( const int idx : g_idx )
      add_iov_blocks( V_iov, *occ, idx, gamma( idx ) * gamma( idx ) );
    return scaled.t() * scaled + V_iov + error_variance;
  }

  arma::mat V = error_variance;
  if ( has_iov ) {
    for ( int i = 0; i < n_omega; ++i ) {
      const arma::vec g  = G_adj_mu.row( i ).t();
      const double w     = omega_iiv( i ) + gamma( i ) * gamma( i );
      V += outer_col( g * std::sqrt( w ) );
    }
  } else {
    arma::mat weighted = G_adj_mu;
    for ( int i = 0; i < n_omega; ++i )
      weighted.row( i ) *= std::sqrt( omega_iiv( i ) );
    V += weighted.t() * weighted;
  }
  return V;
}

// Split variance derivatives into:
//   W    — rank-1 columns (IIV omega; single-occasion gamma reuses G columns)
//   full — dense p x p blocks (multi-occasion gamma, residual sigma)
struct DvParts {
  arma::mat W;
  Rcpp::List full;
};

DvParts build_dV_parts( const arma::mat& G_adj_mu,
                        const arma::vec& gamma,
                        const Rcpp::IntegerVector& occ_col_widths,
                        const Rcpp::List& sigma_derivatives,
                        bool has_iov,
                        const OccasionBlocks* occ ) {
  const int n_omega     = static_cast<int>( G_adj_mu.n_rows );
  const int n_occasions = occ_col_widths.size();
  const std::vector<int> g_idx = active_gamma( gamma );
  const arma::mat G_T = G_adj_mu.t();

  DvParts parts;
  parts.full = Rcpp::List::create();

  if ( has_iov && n_occasions > 1 ) {
    // Multi-occasion: W = omega cols only; gamma derivatives are dense blocks.
    parts.W = G_T.cols( 0, static_cast<arma::uword>( n_omega - 1 ) );
    for ( const int idx : g_idx )
      parts.full.push_back( iov_derivative_block( *occ, idx ) );
  } else if ( has_iov ) {
    // Single occasion: preallocate W = [omega cols | active gamma cols].
    const int n_g = static_cast<int>( g_idx.size() );
    parts.W.set_size( G_T.n_rows, static_cast<arma::uword>( n_omega + n_g ) );
    parts.W.cols( 0, static_cast<arma::uword>( n_omega - 1 ) ) =
      G_T.cols( 0, static_cast<arma::uword>( n_omega - 1 ) );
    for ( int i = 0; i < n_g; ++i )
      parts.W.col( static_cast<arma::uword>( n_omega + i ) ) =
        G_T.col( static_cast<arma::uword>( g_idx[ static_cast<size_t>( i ) ] ) );
  } else {
    parts.W = G_T.cols( 0, static_cast<arma::uword>( n_omega - 1 ) );
  }

  for ( int k = 0; k < sigma_derivatives.size(); ++k )
    parts.full.push_back( Rcpp::as<arma::mat>( sigma_derivatives[ k ] ) );
  return parts;
}

}  // namespace

// Entry point from R: one covariate combination, optional IOV occasions.
// [[Rcpp::export(name = "computePopFimCombo_Rcpp")]]
arma::mat computePopFimCombo_Rcpp( const arma::mat& gradients,
                                   const arma::vec& mu_values,
                                   const arma::vec& omega_iiv,
                                   const arma::vec& gamma,
                                   const arma::mat& error_variance,
                                   const Rcpp::IntegerVector& occ_col_widths,
                                   const Rcpp::List& sigma_derivatives,
                                   bool has_iov ) {
  const int n_omega = static_cast<int>( mu_values.n_elem );

  if ( n_omega <= 0 )
    Rcpp::stop( "computePopFimCombo_Rcpp: at least one parameter is required." );
  if ( static_cast<int>( omega_iiv.n_elem ) != n_omega )
    Rcpp::stop( "computePopFimCombo_Rcpp: omega_iiv length mismatch." );
  if ( static_cast<int>( gamma.n_elem ) != n_omega )
    Rcpp::stop( "computePopFimCombo_Rcpp: gamma length mismatch." );

  if ( static_cast<int>( gradients.n_rows ) < n_omega )
    Rcpp::stop( "computePopFimCombo_Rcpp: gradient rows < parameters." );
  if ( error_variance.n_rows != error_variance.n_cols )
    Rcpp::stop( "computePopFimCombo_Rcpp: error_variance must be square." );
  if ( static_cast<arma::uword>( error_variance.n_rows ) != gradients.n_cols )
    Rcpp::stop( "computePopFimCombo_Rcpp: error variance size mismatch." );

  // Fixed-effect rows only for V / dV (LogNormal-scaled); full gradients
  // (including any extra rows R stacked) are used for MFbeta below.
  arma::mat G_adj_mu = gradients.rows( 0, n_omega - 1 );
  scale_mu_rows( G_adj_mu, mu_values );

  // Split occasion columns once when needed (shared by V and dV assembly).
  OccasionBlocks occ_storage;
  const OccasionBlocks* occ_ptr = nullptr;
  if ( has_iov && occ_col_widths.size() > 1 ) {
    occ_storage = split_by_occasions( G_adj_mu, occ_col_widths );
    occ_ptr = &occ_storage;
  }

  const arma::mat V      = build_V( G_adj_mu, omega_iiv, gamma, error_variance,
                                    occ_col_widths, has_iov, occ_ptr );
  const DvParts dV       = build_dV_parts( G_adj_mu, gamma, occ_col_widths,
                                           sigma_derivatives, has_iov, occ_ptr );
  // Singular V → chol_inv stops (same hard failure as R .safeCholInv path).
  const arma::mat V_inv  = chol_inv( V );
  const arma::mat MFbeta = gradients * V_inv * gradients.t();
  const arma::mat MFVar  = mfvar_mixed( V_inv, dV.W, dV.full );

  // Block-diagonal population FIM (fixed effects | variance parameters).
  const int n_fixed = static_cast<int>( MFbeta.n_rows );
  const int n_var   = static_cast<int>( MFVar.n_rows );
  arma::mat out( n_fixed + n_var, n_fixed + n_var, arma::fill::zeros );
  out.submat( 0, 0, n_fixed - 1, n_fixed - 1 ) = MFbeta;
  out.submat( n_fixed, n_fixed, n_fixed + n_var - 1, n_fixed + n_var - 1 ) = MFVar;
  return out;
}
