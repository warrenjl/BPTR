#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_update(arma::vec delta_trans,
                     double sigma2_phi,
                     double sigma2_regress,
                     arma::mat z_trans,
                     arma::mat ztz,
                     int p_z){
  
  arma::mat cov_eta=arma::inv_sympd(ztz/sigma2_phi+arma::eye(p_z,p_z)/sigma2_regress); 
  arma::vec mean_eta=cov_eta*((z_trans*delta_trans)/sigma2_phi);
  arma::mat ind_norms = arma::randn(1,p_z);
  arma::vec eta = mean_eta+
    trans(ind_norms*arma::chol(cov_eta));
  return(eta);
}