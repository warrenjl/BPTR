#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update(arma::vec y_trans,
                              arma::vec m,
                              arma::vec mu_y_trans,
                              double sigma2_epsilon){

arma::vec dens(sum(m)); dens.fill(0.00);
for(int j = 0; j < sum(m); ++j){
  dens(j) = R::dnorm(y_trans(j),
                     mu_y_trans(j),
                     sqrt(sigma2_epsilon),
                     TRUE);
  }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}

























































