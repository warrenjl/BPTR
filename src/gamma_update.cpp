#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List gamma_update(arma::vec y_trans,
                        arma::mat x,
                        arma::mat x_trans,
                        arma::mat xtx,
                        int p_x,
                        double sigma2_gamma,
                        arma::vec mu_y_trans,
                        arma::vec gamma_old,
                        double sigma2_epsilon_old){
  
arma::vec mean_gamma_piece = y_trans -
                             (mu_y_trans - x*gamma_old);
arma::mat cov_gamma = arma::inv_sympd(xtx/sigma2_epsilon_old + arma::eye(p_x, p_x)/sigma2_gamma);
arma::vec mean_gamma = cov_gamma*(x_trans*mean_gamma_piece)/sigma2_epsilon_old;
arma::mat ind_norms = arma::randn(1, 
                                  p_x);
arma::vec gamma = mean_gamma +
                  trans(ind_norms*arma::chol(cov_gamma));
mu_y_trans = mu_y_trans - 
             x*gamma_old +
             x*gamma;
             
return Rcpp::List::create(Rcpp::Named("gamma") = gamma,
                          Rcpp::Named("mu_y_trans") = mu_y_trans);

}