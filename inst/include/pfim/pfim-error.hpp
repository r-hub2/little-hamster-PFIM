// Residual-error variance and dV/d(sigma) diagonals (Combined1 / Combined2).
//
// form: 1 = Combined1  V = (a + b f^c)^2
//       2 = Combined2  V = a^2 + (b f^c)^2   (PopED additive + proportional SDs)
//
// Returns diagonal vectors only. R wraps them in Matrix::Diagonal and embeds
// multi-output blocks sparsely. Empty d_inter / d_slope means "not estimable"
// (value == 0 or fixed) — must match R helper .pfimSigmaIsEstimable().
// Used by pop-FIM combo assembly; does not touch Bayesian priors.
#pragma once
#include <RcppArmadillo.h>
#include <cmath>

namespace pfim {

inline Rcpp::List residual_error_derivatives( double sigma_inter,
                                              double sigma_slope,
                                              bool inter_fixed,
                                              bool slope_fixed,
                                              double c_error,
                                              const arma::vec& f,
                                              int form ) {
  const arma::uword n = f.n_elem;
  if ( n == 0u )
    Rcpp::stop( "residual_error_derivatives: f has length 0." );
  if ( !( c_error > 0.0 ) )
    Rcpp::stop( "residual_error_derivatives: cError must be positive." );
  if ( form != 1 && form != 2 )
    Rcpp::stop( "residual_error_derivatives: form must be 1 or 2." );

  // f^c once; skip pow when c == 1 (dominant case).
  arma::vec f_term = f;
  if ( c_error != 1.0 )
    f_term = arma::pow( f, c_error );

  arma::vec variance( n );
  arma::vec d_inter( n );
  arma::vec d_slope( n );

  if ( form == 2 ) {
    // Combined2: V = a^2 + (b f^c)^2.
    const arma::vec prop_sd = sigma_slope * f_term;
    variance = sigma_inter * sigma_inter + arma::square( prop_sd );
    d_inter.fill( 2.0 * sigma_inter );
    d_slope = 2.0 * prop_sd % f_term;
  } else {
    // Combined1: V = (a + b f^c)^2.
    const arma::vec base = sigma_inter + sigma_slope * f_term;
    variance = arma::square( base );
    d_inter = 2.0 * base;
    d_slope = 2.0 * base % f_term;
  }

  // Estimable iff non-zero and not fixed (same rule as .pfimSigmaIsEstimable).
  const bool keep_inter = ( sigma_inter != 0.0 ) && !inter_fixed;
  const bool keep_slope = ( sigma_slope != 0.0 ) && !slope_fixed;

  return Rcpp::List::create(
      Rcpp::Named( "variance" ) = variance,
      Rcpp::Named( "d_inter" )  = keep_inter ? d_inter : arma::vec(),
      Rcpp::Named( "d_slope" )  = keep_slope ? d_slope : arma::vec() );
}

}  // namespace pfim
