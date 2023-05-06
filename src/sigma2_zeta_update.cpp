#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_zeta_update(int n,
                          double a_sigma2_zeta,
                          double b_sigma2_zeta,
                          arma::vec zeta){
  
double a_sigma2_zeta_update = 0.50*n +
                              a_sigma2_zeta;
double b_sigma2_zeta_update = 0.50*arma::dot(zeta, zeta) +
                              b_sigma2_zeta;
double sigma2_zeta = 1.00/R::rgamma(a_sigma2_zeta_update,
                                    (1.00/b_sigma2_zeta_update));
  
return(sigma2_zeta);
  
}