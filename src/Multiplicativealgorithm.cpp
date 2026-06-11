// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace arma;

// Multiplicative weights on candidate FIMs (Seurat et al. 2021).
// Stops when max multiplier is within delta of the weighted mean.

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
  (void) show_process;
  std::vector<arma::mat> fim_mats( static_cast<size_t>( n_fim ) );
  for ( int i = 0; i < n_fim; ++i )
    fim_mats[ static_cast<size_t>( i ) ] = Rcpp::as<arma::mat>( fisherMatrices[i] );

  arma::mat sum_weighted_fims( p, p, arma::fill::zeros );
  arma::mat derivative_phi( p, p );
  arma::mat matmult( p, p );
  arma::vec vector_of_multiplier( n_fim );

  int iter = 0;
  for ( ; iter < iteration_init; ++iter ) {

    sum_weighted_fims.zeros();
    for ( int i = 0; i < n_fim; ++i )
      sum_weighted_fims += fim_mats[ static_cast<size_t>( i ) ] * weights[i];

    double det_val = arma::det( sum_weighted_fims );
    if ( !R_finite( det_val ) || det_val <= 0.0 )
      break;

    arma::mat inv_fim;
    const bool inv_ok = arma::inv_sympd( inv_fim, sum_weighted_fims );
    if ( !inv_ok )
      break;

    const double Dcrit = std::pow( det_val, 1.0 / static_cast<double>( p ) );
    derivative_phi     = Dcrit * inv_fim / static_cast<double>( p );

    for ( int i = 0; i < n_fim; ++i ) {
      matmult                  = derivative_phi * fim_mats[ static_cast<size_t>( i ) ];
      vector_of_multiplier[i]  = arma::sum( matmult.diag() );
    }

    arma::vec wm = weights % arma::pow( vector_of_multiplier, lambda );
    weights      = wm / arma::sum( wm );

    if ( vector_of_multiplier.max() <
         ( 1.0 + delta ) * arma::dot( weights, vector_of_multiplier ) )
      break;
  }

  return Rcpp::List::create(
    Rcpp::Named( "weights" )      = weights,
    Rcpp::Named( "iterationEnd" ) = iter
  );
}
