#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec gamma_update(arma::vec y,
                    arma::vec mu_y,
                    arma::vec gamma_old,
                    double sigma2_epsilon,
                    double sigma2_regress,
                    arma::mat x,
                    arma::mat x_trans,
                    arma::mat xtx,
                    int p_x){
  
  
  
  arma::vec mean_gamma_piece=y-(mu_y-x*gamma_old);
  arma::mat cov_gamma=arma::inv_sympd(xtx/sigma2_epsilon+arma::eye(p_x,p_x)/sigma2_regress); 
  arma::vec mean_gamma=cov_gamma*((x_trans*mean_gamma_piece)/sigma2_epsilon);
  arma::mat ind_norms = arma::randn(1, p_x);
  arma::vec gamma = mean_gamma+
                    trans(ind_norms*arma::chol(cov_gamma));
  return(gamma);
}