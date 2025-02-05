#include "RcppArmadillo.h"
#include "BPTR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List zeta0_update(arma::vec y_trans,
                        int n,
                        arma::mat one_star0,
                        arma::mat one_star0_trans,
                        arma::mat one_star0_t_one_star0,
                        arma::vec mu_y_trans,
                        arma::vec zeta0_old,
                        double sigma2_zeta0_old,
                        double sigma2_epsilon_old){
  
arma::vec mean_zeta0_piece = y_trans - 
                             (mu_y_trans - one_star0*zeta0_old);
arma::mat cov_zeta0 = arma::inv_sympd(one_star0_t_one_star0/sigma2_epsilon_old + arma::eye(n, n)/sigma2_zeta0_old); 
arma::vec mean_zeta0 = cov_zeta0*(one_star0_trans*mean_zeta0_piece)/sigma2_epsilon_old;
arma::mat ind_norms = arma::randn(1, 
                                  n);
arma::vec zeta0 = mean_zeta0 +
                  trans(ind_norms*arma::chol(cov_zeta0));

mu_y_trans = mu_y_trans - 
             one_star0*zeta0_old +
             one_star0*zeta0;

return Rcpp::List::create(Rcpp::Named("zeta0") = zeta0,
                          Rcpp::Named("mu_y_trans") = mu_y_trans);

}



