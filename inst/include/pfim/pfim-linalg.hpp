// PFIM linear-algebra kernels shared by FIM assembly and discrete optimizers.
//
// Header-only Armadillo helpers: Cholesky inverse, safe solve, and the
// variance-parameter FIM blocks mfvar / mfvar_mixed (see R .computeMFVar).
//
// Numerical contract: try_* helpers return empty / false (non-throwing);
// chol_inv / safe_solve stop on failure — callers that must degrade gracefully
// (FW / Mult mixture FIMs) use try_chol_inv_logdet. No Bayesian prior here.
#pragma once
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

namespace pfim {

namespace {

// Invert A via solve(A, I); return empty matrix on failure or poor residual.
inline arma::mat solve_invert( const arma::mat& A ) {
  arma::mat out;
  const arma::mat I = arma::eye( A.n_rows, A.n_cols );
  if ( !arma::solve( out, A, I, arma::solve_opts::no_approx ) )
    return arma::mat();
  // Relative residual: ||A A^{-1} - I||_inf / (||A||_inf * ||A^{-1}||_inf + eps).
  const double abs_res = arma::norm( A * out - I, "inf" );
  const double scale =
      arma::norm( A, "inf" ) * arma::norm( out, "inf" ) +
      std::numeric_limits<double>::epsilon();
  if ( abs_res > 1e-8 * scale )
    return arma::mat();
  return out;
}

}  // namespace

// Outer product v * v^T (rank-1 contribution to dV).
inline arma::mat outer_col( const arma::vec& v ) {
  return v * v.t();
}

// SPD inverse and log-det via one Cholesky factorization (non-throwing).
inline bool try_chol_inv_logdet( const arma::mat& M,
                                 arma::mat& inv_out,
                                 double& log_det_out ) {
  // NA/Inf fail Armadillo's symmetry check (NA != NA) even after Msym.
  if ( !M.is_finite() )
    return false;

  const arma::mat Msym = 0.5 * ( M + M.t() );
  arma::mat L;
  if ( !arma::chol( L, Msym, "lower" ) )
    return false;

  const arma::vec d = arma::diagvec( L );
  if ( arma::any( d <= 0.0 ) )
    return false;

  log_det_out = 2.0 * arma::as_scalar( arma::sum( arma::log( d ) ) );
  if ( !R_finite( log_det_out ) )
    return false;

  // inv = (L^{-T}) (L^{-1}) from the upper-triangular factor.
  const arma::mat Ui = arma::inv( arma::trimatu( L.t() ) );
  inv_out            = Ui * Ui.t();
  return true;
}

// chol2inv(chol(M)); empty matrix when M is not SPD (non-throwing).
inline arma::mat try_chol_inv( const arma::mat& M ) {
  arma::mat inv;
  double log_det = 0.0;
  if ( !try_chol_inv_logdet( M, inv, log_det ) )
    return arma::mat();
  return inv;
}

// chol2inv(chol(V)); stops when Cholesky fails (.safeCholInv).
// An all-zero or non-PD matrix is an error (no silent eps·I inverse).
inline arma::mat chol_inv( const arma::mat& V ) {
  if ( !V.is_finite() )
    Rcpp::stop( "chol_inv: matrix has non-finite entries." );
  arma::mat out = try_chol_inv( V );
  if ( out.n_elem > 0 )
    return out;
  Rcpp::stop( "chol_inv: Cholesky failed (matrix is zero or not positive definite)." );
}

// solve(M); fails when M is singular (.safeSolve).
inline arma::mat safe_solve( const arma::mat& M, bool* regularized = nullptr ) {
  if ( M.n_rows != M.n_cols )
    Rcpp::stop( "safe_solve: matrix must be square." );

  const arma::mat out = solve_invert( M );
  if ( out.n_elem == 0 )
    Rcpp::stop( "safe_solve: inversion failed (matrix is singular or ill-conditioned)." );
  if ( regularized )
    *regularized = false;
  return out;
}

// 1/2 Tr(V^{-1} dV_i V^{-1} dV_j) (.computeMFVar) — full dV list path.
inline arma::mat mfvar( const arma::mat& V_inv, const Rcpp::List& dV_list ) {
  const int n = dV_list.size();
  if ( n == 0 )
    return arma::mat( 0, 0 );

  const arma::uword p = V_inv.n_rows;
  // Convert each List entry once; reuse for T_k and the Frobenius products.
  std::vector<arma::mat> dV( static_cast<size_t>( n ) );
  std::vector<arma::mat> T( static_cast<size_t>( n ) );

  for ( int k = 0; k < n; ++k ) {
    dV[ static_cast<size_t>( k ) ] = Rcpp::as<arma::mat>( dV_list[ k ] );
    const arma::mat& dVk = dV[ static_cast<size_t>( k ) ];
    if ( dVk.n_rows != p || dVk.n_cols != p )
      Rcpp::stop( "mfvar: dV_list[[%d]] is not %dx%d.", k + 1, p, p );
    T[ static_cast<size_t>( k ) ] = V_inv * dVk * V_inv;
  }

  arma::mat out( n, n, arma::fill::zeros );
  for ( int j = 0; j < n; ++j ) {
    const arma::mat& dVj = dV[ static_cast<size_t>( j ) ];
    for ( int i = 0; i <= j; ++i ) {
      const double val = 0.5 * arma::accu( T[ static_cast<size_t>( i ) ] % dVj );
      out( i, j ) = out( j, i ) = val;
    }
  }
  return out;
}

// Rank-1 columns W (dV_k = w_k w_k^T) plus optional full p x p derivatives (IOV blocks, sigma).
// Tr(V^{-1} w_i w_i^T V^{-1} w_j w_j^T) = (w_i^T V^{-1} w_j)^2.
inline arma::mat mfvar_mixed( const arma::mat& V_inv,
                              const arma::mat& W,
                              const Rcpp::List& dV_full ) {
  const int n1 = static_cast<int>( W.n_cols );
  const int n2 = dV_full.size();
  const int n  = n1 + n2;
  if ( n == 0 )
    return arma::mat( 0, 0 );

  const arma::uword p = V_inv.n_rows;
  if ( W.n_rows != p )
    Rcpp::stop( "mfvar_mixed: W has %d rows, expected %d.", W.n_rows, p );

  arma::mat out( n, n, arma::fill::zeros );
  arma::mat VinW;
  // Rank-1 block: out_ij = 1/2 (w_i^T V^{-1} w_j)^2.
  if ( n1 > 0 ) {
    VinW = V_inv * W;
    const arma::mat A = W.t() * VinW;
    for ( int j = 0; j < n1; ++j ) {
      for ( int i = 0; i <= j; ++i ) {
        const double val = 0.5 * A( i, j ) * A( i, j );
        out( i, j ) = out( j, i ) = val;
      }
    }
  }

  if ( n2 == 0 )
    return out;

  // Convert full dV entries once for T_k, dense block, and rank-1 cross terms.
  std::vector<arma::mat> dV( static_cast<size_t>( n2 ) );
  std::vector<arma::mat> T( static_cast<size_t>( n2 ) );
  for ( int k = 0; k < n2; ++k ) {
    dV[ static_cast<size_t>( k ) ] = Rcpp::as<arma::mat>( dV_full[ k ] );
    const arma::mat& dVk = dV[ static_cast<size_t>( k ) ];
    if ( dVk.n_rows != p || dVk.n_cols != p )
      Rcpp::stop( "mfvar_mixed: dV_full[[%d]] is not %dx%d.", k + 1, p, p );
    T[ static_cast<size_t>( k ) ] = V_inv * dVk * V_inv;
  }

  for ( int j = 0; j < n2; ++j ) {
    const arma::mat& dVj = dV[ static_cast<size_t>( j ) ];
    for ( int i = 0; i <= j; ++i ) {
      const double val = 0.5 * arma::accu( T[ static_cast<size_t>( i ) ] % dVj );
      out( n1 + i, n1 + j ) = out( n1 + j, n1 + i ) = val;
    }
    if ( n1 > 0 ) {
      for ( int i = 0; i < n1; ++i ) {
        const arma::vec wi = VinW.col( static_cast<arma::uword>( i ) );
        const double val = 0.5 * arma::as_scalar( wi.t() * dVj * wi );
        out( i, n1 + j ) = val;
        out( n1 + j, i ) = val;
      }
    }
  }

  return out;
}

}  // namespace pfim
