#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List v_update(arma::mat spatial_dists,
                         int m,
                         double a_phi,
                         double b_phi,
                         arma::vec eta,
                         double tau2,
                         double phi_old,
                         Rcpp::List spatial_corr_info,
                         double metrop_var_phi_trans,
                         int acctot_phi_trans){
 
/*Second*/
arma::vec mu_y_trans_old = mu_y_trans;
arma::vec psi_old = psi;
arma::vec theta_old = theta;

double second = 0.50*log_deter_inv_old - 
                0.50*dot(eta, (corr_inv_old*eta))/tau2 +
                a_phi*phi_trans_old -
                b_phi*exp(phi_trans_old);

/*First*/
double phi_trans = R::rnorm(phi_trans_old, 
                            sqrt(metrop_var_phi_trans));
double phi = exp(phi_trans);
spatial_corr_info = spatial_corr_fun(m,
                                     spatial_dists,
                                     phi);
arma::mat corr_inv = spatial_corr_info[0];
double log_deter_inv = spatial_corr_info[1];

double first = 0.50*log_deter_inv - 
               0.50*dot(eta, (corr_inv*eta))/tau2 +
               a_phi*phi_trans -
               b_phi*exp(phi_trans);

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  
  phi = phi_old;
  spatial_corr_info = spatial_corr_info_old;
  acc = 0;
  
  }
acctot_phi_trans = acctot_phi_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("phi") = phi,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                          Rcpp::Named("spatial_corr_info") = spatial_corr_info);

}



