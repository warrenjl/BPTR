#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double alpha_update(int d,
                    double a_alpha,
                    double b_alpha,
                    arma::vec V){
  
double alpha_a_update = a_alpha + 
                        d - 
                        1.00;
double alpha_b_update = b_alpha - 
                        sum(log(1.00 - V));
double alpha = R::rgamma(alpha_a_update,
                         (1.00/alpha_b_update));
  
return(alpha);
  
}