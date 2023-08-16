#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List zeta1_update(arma::vec y_trans,
                        int sum_r,
                        arma::mat one_star1,
                        arma::mat one_star1_trans,
                        arma::mat one_star1_t_one_star1,
                        arma::vec mu_y_trans,
                        arma::vec zeta1_old,
                        double sigma2_zeta1_old,
                        double sigma2_epsilon_old){
  
arma::vec mean_zeta1_piece = y_trans - 
                             (mu_y_trans - one_star1*zeta1_old);
arma::mat cov_zeta1 = arma::inv_sympd(one_star1_t_one_star1/sigma2_epsilon_old + arma::eye(sum_r, sum_r)/sigma2_zeta1_old); 
arma::vec mean_zeta1 = cov_zeta1*(one_star1_trans*mean_zeta1_piece)/sigma2_epsilon_old;
arma::mat ind_norms = arma::randn(1, 
                                  sum_r);
arma::vec zeta1 = mean_zeta1 +
                  trans(ind_norms*arma::chol(cov_zeta1));

mu_y_trans = mu_y_trans - 
             one_star1*zeta1_old +
             one_star1*zeta1;

return Rcpp::List::create(Rcpp::Named("zeta1") = zeta1,
                          Rcpp::Named("mu_y_trans") = mu_y_trans);

}



