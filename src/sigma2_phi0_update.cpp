#include "RcppArmadillo.h"
#include "BPTR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_phi0_update(int n,
                          double a_sigma2_phi0,
                          double b_sigma2_phi0,
                          arma::vec phi0){
  
double a_sigma2_phi0_update = 0.50*n + 
                              a_sigma2_phi0;
double b_sigma2_phi0_update = 0.50*arma::dot(phi0, phi0) +
                              b_sigma2_phi0;
double sigma2_phi0 = 1.00/R::rgamma(a_sigma2_phi0_update,
                                    (1.00/b_sigma2_phi0_update));
  
return(sigma2_phi0);

}