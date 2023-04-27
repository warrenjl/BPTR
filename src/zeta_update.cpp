#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec zeta_update(arma::vec y,
                      arma::vec mu_y,
                      arma::vec zeta_old,
                      double sigma2_epsilon,
                      double sigma2_zeta,
                      arma::mat one_star,
                      arma::mat onestar_trans,
                      arma::mat onetone,
                      int n){
  
 arma::vec mean_zeta_piece=y-(mu_y-one_star*zeta_old);
  arma::mat cov_zeta=arma::inv_sympd(onetone/sigma2_epsilon+arma::eye(n,n)/sigma2_zeta); 
 arma::vec mean_zeta=cov_zeta*((onestar_trans*mean_zeta_piece)/sigma2_epsilon);
 arma::mat ind_norms = arma::randn(1, n);
 arma::vec zeta = mean_zeta+
      trans(ind_norms*arma::chol(cov_zeta));
  return(zeta);
}