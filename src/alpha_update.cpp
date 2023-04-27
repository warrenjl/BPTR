#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double alpha_update(arma::vec V,
                    double alpha_a,
                    double alpha_b,
                    int d){
  
  double alpha_a_update = alpha_a + d-1.00;
  double alpha_b_update = (alpha_b-arma::sum(log(1.00-V)));
  double alpha = R::rgamma(alpha_a_update,
                           1.00/alpha_b_update);
  
  return(alpha);
  
}