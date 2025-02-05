#include "RcppArmadillo.h"
#include "BPTR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_zeta1_update(int sum_r,
                           double a_sigma2_zeta1,
                           double b_sigma2_zeta1,
                           arma::vec zeta1){
  
double a_sigma2_zeta1_update = 0.50*sum_r +
                               a_sigma2_zeta1;
double b_sigma2_zeta1_update = 0.50*arma::dot(zeta1, zeta1) +
                               b_sigma2_zeta1;
double sigma2_zeta1 = 1.00/R::rgamma(a_sigma2_zeta1_update,
                                     (1.00/b_sigma2_zeta1_update));
  
return(sigma2_zeta1);
  
}