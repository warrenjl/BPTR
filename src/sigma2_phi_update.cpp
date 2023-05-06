#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_phi_update(int n,
                         arma::mat z,
                         double a_sigma2_phi,
                         double b_sigma2_phi,
                         arma::vec delta_trans,
                         arma::vec eta){
  
double a_sigma2_phi_update = 0.50*n + 
                             a_sigma2_phi;
double b_sigma2_phi_update = 0.50*arma::dot((delta_trans - z*eta), (delta_trans - z*eta)) +
                             b_sigma2_phi;
double sigma2_phi = 1.00/R::rgamma(a_sigma2_phi_update,
                                   (1.00/b_sigma2_phi_update));
  
return(sigma2_phi);

}