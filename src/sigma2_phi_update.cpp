#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_phi_update(arma::mat z,
                         arma::vec delta_trans,
                         arma::vec eta,
                         double a_sigma2_phi,
                         double b_sigma2_phi,
                         int n){
  
  double a_sigma2_phi_update = 0.50*n+
    a_sigma2_phi;
  double b_sigma2_phi_update = 0.50*arma::dot(delta_trans-(z*eta),delta_trans-(z*eta))+
    b_sigma2_phi;
  double sigma2_phi = 1.00/R::rgamma(a_sigma2_phi_update,
                                      (1.00/b_sigma2_phi_update));
  
  return(sigma2_phi);

}