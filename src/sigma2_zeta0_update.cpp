#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_zeta0_update(int n,
                           double a_sigma2_zeta0,
                           double b_sigma2_zeta0,
                           arma::vec zeta0){
  
double a_sigma2_zeta0_update = 0.50*n +
                               a_sigma2_zeta0;
double b_sigma2_zeta0_update = 0.50*arma::dot(zeta0, zeta0) +
                               b_sigma2_zeta0;
double sigma2_zeta0 = 1.00/R::rgamma(a_sigma2_zeta0_update,
                                     (1.00/b_sigma2_zeta0_update));
  
return(sigma2_zeta0);
  
}