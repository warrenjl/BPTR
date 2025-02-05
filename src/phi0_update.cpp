#include "RcppArmadillo.h"
#include "BPTR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec phi0_update(arma::mat z,
                      arma::mat one_star2_trans,
                      arma::mat one_star2_t_one_star2,
                      int n,
                      arma::vec delta_trans,
                      arma::vec eta,
                      double sigma2_phi0_old,
                      double sigma2_phi1_old){
  
arma::mat cov_phi0 = arma::inv_sympd(one_star2_t_one_star2/sigma2_phi1_old + arma::eye(n, n)/sigma2_phi0_old); 
arma::vec mean_phi0 = cov_phi0*(one_star2_trans*(delta_trans - z*eta))/sigma2_phi1_old;
arma::mat ind_norms = arma::randn(1, 
                                  n);
arma::vec phi0 = mean_phi0 +
                 trans(ind_norms*arma::chol(cov_phi0));

return phi0;

}
