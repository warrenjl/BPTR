#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List zeta_update(arma::vec y_trans,
                       int n,
                       arma::mat one_star,
                       arma::mat one_star_trans,
                       arma::mat one_star_t_one_star,
                       arma::vec mu_y_trans,
                       arma::vec zeta_old,
                       double sigma2_zeta_old,
                       double sigma2_epsilon_old){
  
arma::vec mean_zeta_piece = y_trans - 
                            (mu_y_trans - one_star*zeta_old);
arma::mat cov_zeta = arma::inv_sympd(one_star_t_one_star/sigma2_epsilon_old + arma::eye(n, n)/sigma2_zeta_old); 
arma::vec mean_zeta = cov_zeta*(one_star_trans*mean_zeta_piece)/sigma2_epsilon_old;
arma::mat ind_norms = arma::randn(1, 
                                  n);
arma::vec zeta = mean_zeta +
                 trans(ind_norms*arma::chol(cov_zeta));

mu_y_trans = mu_y_trans - 
             one_star*zeta_old +
             one_star*zeta;

return Rcpp::List::create(Rcpp::Named("zeta") = zeta,
                          Rcpp::Named("mu_y_trans") = mu_y_trans);

}



