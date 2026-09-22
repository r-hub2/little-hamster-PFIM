// Rcpp exports for variance-parameter FIM blocks (mfvar / mfvar_mixed).
// Implementations live in pfim/pfim-linalg.hpp; mirrors R `.computeMFVar`.
// [[Rcpp::depends(RcppArmadillo)]]
#include <pfim/pfim-linalg.hpp>

// Full-matrix path: dV_list[[k]] are p x p derivatives of V.
// [[Rcpp::export(name = "computeMFVar_Rcpp")]]
arma::mat computeMFVar_Rcpp( const arma::mat& V_inv, const Rcpp::List& dV_list ) {
  return pfim::mfvar( V_inv, dV_list );
}

// Mixed path: rank-1 columns in W plus optional full dV blocks (IOV / sigma).
// [[Rcpp::export(name = "computeMFVar_mixed_Rcpp")]]
arma::mat computeMFVar_mixed_Rcpp( const arma::mat& V_inv,
                                   const arma::mat& W,
                                   const Rcpp::List& dV_full ) {
  return pfim::mfvar_mixed( V_inv, W, dV_full );
}
