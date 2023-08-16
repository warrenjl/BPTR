#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_phi1_update(int sum_r,
                          arma::mat z,
                          arma::mat one_star2,
                          double a_sigma2_phi1,
                          double b_sigma2_phi1,
                          arma::vec delta_trans,
                          arma::vec eta,
                          arma::vec phi0){
  
double a_sigma2_phi1_update = 0.50*sum_r + 
                              a_sigma2_phi1;
double b_sigma2_phi1_update = 0.50*arma::dot((delta_trans - z*eta - one_star2*phi0), (delta_trans - z*eta - one_star2*phi0)) +
                              b_sigma2_phi1;
double sigma2_phi1 = 1.00/R::rgamma(a_sigma2_phi1_update,
                                   (1.00/b_sigma2_phi1_update));
  
return(sigma2_phi1);

}