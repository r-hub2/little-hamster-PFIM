// Rcpp exports for safe matrix inverses used by FIM / FD helpers.
// Implementations live in pfim/pfim-linalg.hpp (.safeCholInv / .safeSolve).
// [[Rcpp::depends(RcppArmadillo)]]
#include <pfim/pfim-linalg.hpp>

// Cholesky-based inverse of an SPD matrix (stops if not positive definite).
// [[Rcpp::export(name = "chol_inv_Rcpp")]]
arma::mat chol_inv_Rcpp( const arma::mat& V ) {
  return pfim::chol_inv( V );
}

// General square solve / inverse (stops if singular or ill-conditioned).
// [[Rcpp::export(name = "safe_solve_Rcpp")]]
arma::mat safe_solve_Rcpp( const arma::mat& M ) {
  return pfim::safe_solve( M );
}
