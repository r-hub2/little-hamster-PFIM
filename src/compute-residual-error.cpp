// Rcpp export for Combined1 / Combined2 residual variance diagonals.
//
// form: 1 = Combined1  V=(a+b f^c)^2 ; 2 = Combined2  V=a^2+(b f^c)^2.
// Estimable sigma rule matches ModelError.R: value != 0 && !fixed
// (non-zero fixed sigmas still enter V, never FIM columns).
// Implementation: pfim/pfim-error.hpp.
// [[Rcpp::depends(RcppArmadillo)]]
#include <pfim/pfim-error.hpp>

// [[Rcpp::export(name = "residualErrorDerivatives_Rcpp")]]
Rcpp::List residualErrorDerivatives_Rcpp( double sigmaInter,
                                          double sigmaSlope,
                                          bool sigmaInterFixed,
                                          bool sigmaSlopeFixed,
                                          double cError,
                                          const arma::vec& f,
                                          int form ) {
  return pfim::residual_error_derivatives(
      sigmaInter, sigmaSlope, sigmaInterFixed, sigmaSlopeFixed, cError, f, form );
}
