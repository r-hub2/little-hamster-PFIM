// Multiplicative algorithm for D-optimal weights on a fixed candidate set
// (Pukelsheim / classical continuous design theory).
//
// Contract with R:
//   - fisherMatrices[[i]] are dense p×p elementary FIMs (already assembled in R)
//   - weights live on the simplex; returned weights are the last certified point
//   - Bayesian prior is not added here (R may fold it into F_i before the call)
//
// Iteration (one pass):
//   1. Mixture FIM:     F(w) = sum_i w_i F_i
//   2. D-criterion:     D = det(F)^{1/p} = exp(log_det / p)   // stable scale
//   3. Sensitivity:     phi' = (D/p) F^{-1}
//                       m_i  = tr( phi' F_i )                 // directional gains
//   4. Optimality test: stop if max(m) < (1+delta) * sum(w_i m_i)
//                       certify on the current weights before updating; otherwise
//                       the returned w would no longer be the certified point
//   5. Weight update:   w = w * m^lambda / ||.||_1            // stay on the simplex
//
// Uses exp(log_det/p), never raw det(F), both for the formula and to avoid Inf.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <pfim/pfim-linalg.hpp>
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using pfim::try_chol_inv_logdet;

// [[Rcpp::export]]
Rcpp::List MultiplicativeAlgorithm_Rcpp(
    Rcpp::List fisherMatrices,
    int        n_fim,
    arma::vec  weights,
    int        p,
    double     lambda,
    double     delta,
    int        iteration_init,
    bool       show_process = false )
{
  if ( n_fim <= 0 || p <= 0 || iteration_init < 0 )
    Rcpp::stop( "MultiplicativeAlgorithm: non-positive dimensions." );
  if ( static_cast<int>( fisherMatrices.size() ) < n_fim )
    Rcpp::stop( "MultiplicativeAlgorithm: fisherMatrices length mismatch." );
  if ( static_cast<int>( weights.n_elem ) != n_fim )
    Rcpp::stop( "MultiplicativeAlgorithm: weights length mismatch." );

  std::vector<arma::mat> fim_mats( static_cast<size_t>( n_fim ) );
  for ( int i = 0; i < n_fim; ++i ) {
    fim_mats[ static_cast<size_t>( i ) ] = Rcpp::as<arma::mat>( fisherMatrices[i] );
    const arma::mat& M = fim_mats[ static_cast<size_t>( i ) ];
    if ( M.n_rows != static_cast<arma::uword>( p ) || M.n_cols != static_cast<arma::uword>( p ) )
      Rcpp::stop( "MultiplicativeAlgorithm: FIM %d is not %dx%d.", i + 1, p, p );
  }

  if ( show_process )
    Rcpp::Rcout << "Multiplicative algorithm (up to " << iteration_init << " iterations)\n";

  arma::mat sum_weighted_fims( p, p, arma::fill::zeros );
  arma::mat derivative_phi( p, p );
  arma::mat matmult( p, p );
  arma::vec vector_of_multiplier( n_fim );

  int iter = 0;
  bool converged = false;
  bool singular_fim = false;
  for ( ; iter < iteration_init; ++iter ) {
    Rcpp::checkUserInterrupt();

    // F(w) = sum_i w_i F_i
    sum_weighted_fims.zeros();
    for ( int i = 0; i < n_fim; ++i )
      sum_weighted_fims += fim_mats[ static_cast<size_t>( i ) ] * weights[i];

    arma::mat inv_fim;
    double log_det = 0.0;
    if ( !try_chol_inv_logdet( sum_weighted_fims, inv_fim, log_det ) ) {
      singular_fim = true;
      break;
    }

    // D = exp(log_det / p); phi' = (D/p) F^{-1}; m_i = tr(phi' F_i)
    const double Dcrit = std::exp( log_det / static_cast<double>( p ) );
    derivative_phi     = Dcrit * inv_fim / static_cast<double>( p );

    for ( int i = 0; i < n_fim; ++i ) {
      matmult                 = derivative_phi * fim_mats[ static_cast<size_t>( i ) ];
      vector_of_multiplier[i] = arma::sum( matmult.diag() );
    }

    // Near-optimality on current w (m_i from F(w)): check before the update so
    // returned weights are exactly the certified simplex point.
    if ( vector_of_multiplier.max() <
         ( 1.0 + delta ) * arma::dot( weights, vector_of_multiplier ) ) {
      converged = true;
      break;
    }

    // Multiplicative update on the simplex: w = w * m^lambda / ||.||_1
    arma::vec wm = weights % arma::pow( vector_of_multiplier, lambda );
    const double wm_sum = arma::sum( wm );
    if ( wm_sum <= 0.0 || !R_finite( wm_sum ) ) {
      singular_fim = true;
      break;
    }
    weights = wm / wm_sum;
  }

  if ( show_process )
    Rcpp::Rcout << "Multiplicative algorithm: " << iter
                << " iteration(s), converged = " << ( converged ? "yes" : "no" ) << "\n";

  return Rcpp::List::create(
    Rcpp::Named( "weights" )     = weights,
    Rcpp::Named( "iterations" )  = iter,
    Rcpp::Named( "converged" )   = converged,
    Rcpp::Named( "singularFim" ) = singular_fim
  );
}
