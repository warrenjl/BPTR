#include "RcppArmadillo.h"
#include "BPTR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_update(arma::mat z_trans,
                     arma::mat ztz,
                     arma::mat one_star2,
                     int p_z,
                     double sigma2_eta,
                     arma::vec delta_trans,
                     arma::vec phi0_old,
                     double sigma2_phi1_old){
  
arma::mat cov_eta = arma::inv_sympd(ztz/sigma2_phi1_old + arma::eye(p_z, p_z)/sigma2_eta); 
arma::vec mean_eta = cov_eta*(z_trans*(delta_trans - one_star2*phi0_old))/sigma2_phi1_old;
arma::mat ind_norms = arma::randn(1, 
                                  p_z);
arma::vec eta = mean_eta +
                trans(ind_norms*arma::chol(cov_eta));

return eta;

}
