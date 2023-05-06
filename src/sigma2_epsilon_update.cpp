#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(arma::vec y_trans,
                             arma::vec m,
                             arma::vec mu_y_trans,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon){
  
double a_sigma2_epsilon_update = 0.50*sum(m) +
                                 a_sigma2_epsilon;
double b_sigma2_epsilon_update = 0.50*arma::dot((y_trans - mu_y_trans), (y_trans - mu_y_trans)) +
                                 b_sigma2_epsilon;
double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                       (1.00/b_sigma2_epsilon_update));
  
return(sigma2_epsilon);
    
}