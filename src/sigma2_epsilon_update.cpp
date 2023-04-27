#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(arma::vec y,
                             arma::vec mu_y,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec m){
  
  double a_sigma2_epsilon_update = 0.50*arma::sum(m)+
                                   a_sigma2_epsilon;
  double b_sigma2_epsilon_update = 0.50*arma::dot((y-mu_y),(y-mu_y))+
                                   b_sigma2_epsilon;
  double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                         (1.00/b_sigma2_epsilon_update));
  
  return(sigma2_epsilon);
    
}