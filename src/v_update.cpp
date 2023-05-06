#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List V_update(arma::vec y_trans,
                    arma::vec m,
                    arma::mat x,
                    arma::mat one_star,
                    int d,
                    arma::vec mu_y_trans_old,
                    arma::vec gamma,
                    arma::vec zeta,
                    arma::vec V_old,
                    arma::vec theta_old,
                    arma::vec psi_less_old,
                    arma::vec psi_old,
                    arma::vec full_theta_old,
                    double alpha_old,
                    arma::mat z_delta_old,
                    double sigma2_epsilon_old,
                    arma::vec metrop_V,
                    arma::vec acctot_V){
  
arma::vec dens(sum(m)); dens.fill(0.00);
double numer = 0.00;
double denom = 0.00;

arma::vec mu_y_trans = mu_y_trans_old;
arma::vec V = V_old;
arma::vec theta = theta_old;
arma::vec full_theta = full_theta_old;
arma::vec psi_less = psi_less_old;
arma::vec psi = psi_old;

for(int j = 0; j < (d-1); ++j){
  
   mu_y_trans_old = mu_y_trans;
   psi_less_old = psi_less;
   psi_old = psi;
   theta_old = theta;
   full_theta_old = full_theta;
   
   /*Second*/
   for(int k = 0; k < sum(m); ++k){
      dens(k) = R::dnorm(y_trans(k),
                         mu_y_trans_old(k),
                         sqrt(sigma2_epsilon_old),
                         TRUE);
      }
   denom = sum(dens) +
           R::dbeta(V_old(j),
                    1.00,
                    alpha_old,
                    TRUE);

   /*First*/
   V(j) = R::runif((V_old(j) - metrop_V(j)),
                   (V_old(j) + metrop_V(j)));
   if(V(j) < 0.00){
     V(j) = abs(V(j));
     }
   if(V(j) > 1.00){
     V(j) = 2.00 -
            V(j);
     } 
   
   psi_less(0) = 1.00;
   if(d > 2){
     psi_less.subvec(1, (d-2)) = arma::cumprod(1.00 - V.subvec(0, (d-3)));
     }
   psi_less = psi_less%V;
   psi.subvec(0, (d-2)) = psi_less;
   psi(d-1) = 1.00 - 
              sum(psi_less);
   theta(0) = psi(0);
   for(int k = 1; k < d; ++k){
      theta(k) = theta(k-1) +
                 psi(k);
      }

   full_theta.subvec(1, d) = theta;
   
   mu_y_trans = mu_y_trans - 
                log(z_delta_old*full_theta_old) +
                log(z_delta_old*full_theta);
   
   for(int k = 0; k < sum(m); ++k){
      dens(k) = R::dnorm(y_trans(k),
                         mu_y_trans(k),
                         sqrt(sigma2_epsilon_old),
                         TRUE);
      }
   numer = sum(dens) +
           R::dbeta(V(j),
                    1.00,
                    alpha_old,
                    TRUE);

   /*Decision*/
   double ratio = exp(numer - denom);   
   double acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
  
     V(j) = V_old(j);
     theta = theta_old;
     full_theta = full_theta_old;
     psi_less = psi_less_old;
     psi = psi_old;
     mu_y_trans = mu_y_trans_old;
     acc = 0;
  
     }
   acctot_V(j) = acctot_V(j) + 
                 acc;
   
   }

return Rcpp::List::create(Rcpp::Named("V") = V,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("psi_less") = psi_less,
                          Rcpp::Named("psi") = psi,
                          Rcpp::Named("full_theta") = full_theta,
                          Rcpp::Named("mu_y_trans") = mu_y_trans,
                          Rcpp::Named("acctot_V") = acctot_V);

}



