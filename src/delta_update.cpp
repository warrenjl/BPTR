#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update(arma::vec y_trans,
                        int sum_r,
                        arma::vec m,
                        arma::mat t,
                        arma::mat z,
                        arma::vec c,
                        int d,
                        arma::vec a0,
                        arma::vec max_a1_a2,
                        arma::vec choose_vec,
                        Rcpp::IntegerVector power1,
                        Rcpp::IntegerVector power2,
                        Rcpp::IntegerVector phi0_pointer,
                        arma::vec mu_y_trans_old,
                        arma::vec theta,
                        arma::vec full_theta,
                        arma::vec delta_old,
                        arma::vec delta_trans_old,
                        arma::mat z_delta_old,
                        arma::vec eta_old,
                        arma::vec phi0_old,
                        double sigma2_phi1_old,
                        double sigma2_epsilon_old,
                        arma::vec metrop_var_delta,
                        arma::vec acctot_delta){
  
arma::vec dens(sum(m)); dens.fill(0.00);
double numer = 0.00;
double denom = 0.00;
arma::vec vec1(d+1); vec1.fill(0.00);
arma::vec vec2(d+1); vec2.fill(0.00);

arma::vec mu_y_trans = mu_y_trans_old;
arma::vec delta = delta_old;
arma::vec delta_trans = delta_trans_old;
arma::mat z_delta = z_delta_old;

int counter0; counter0 = 0;
int counter1; counter1 = 0;
for(int j = 0; j < sum_r; ++j){
  
   if(j == 0){
     counter0 = 0;
     }
   if(j > 0){
     counter0 = sum(m.subvec(0, (j-1)));
     }
   
   delta_trans_old(j) = delta_trans(j);
   z_delta_old = z_delta;
   mu_y_trans_old = mu_y_trans;
  
   /*Second*/
   dens.fill(0.00);
   for(int k = counter0; k < sum(m.subvec(0,j)); ++k){
      dens(k) = R::dnorm(y_trans(k),
                         mu_y_trans_old(k),
                         sqrt(sigma2_epsilon_old),
                         TRUE);
      }
   denom = sum(dens) +
           R::dnorm(delta_trans_old(j),
                    (dot(z.row(j), eta_old) + phi0_old(phi0_pointer(j) - 1)),
                    sqrt(sigma2_phi1_old),
                    TRUE);
   
   /*First*/
   delta_trans(j) = R::rnorm(delta_trans_old(j),
                             sqrt(metrop_var_delta(j)));
   delta(j) = (a0(j) - max_a1_a2(j))/(1.00 + exp(delta_trans(j)));
  
   counter1 = counter0;
   
   for(int k = 0; k < m(j); ++k){
    
      counter1 = counter1 +
                 1;
      for(int l = 0; l < (d+1); ++l){
                                 
         vec1(l) = pow((t(j,k) + delta(j))/max(c), power1(l));
         vec2(l) = pow((1.00 - (t(j,k) + delta(j))/max(c)), power2(l));
      
         }
      z_delta.row(counter1 - 1) = trans(choose_vec%vec1%vec2);
      
      }
   
   mu_y_trans.subvec(counter0, (sum(m.subvec(0,j)) - 1)) = (mu_y_trans.subvec(counter0, (sum(m.subvec(0,j)) - 1)) - log(z_delta_old.rows(counter0, (sum(m.subvec(0,j)) - 1))*full_theta)) +
                                                           log(z_delta.rows(counter0, (sum(m.subvec(0,j)) - 1))*full_theta);
   
   dens.fill(0.00);
   for(int k = counter0; k < sum(m.subvec(0,j)); ++k){
      dens(k) = R::dnorm(y_trans(k),
                         mu_y_trans(k),
                         sqrt(sigma2_epsilon_old),
                         TRUE);
      }
   numer = sum(dens) +
           R::dnorm(delta_trans(j),
                    (dot(z.row(j), eta_old) + phi0_old(phi0_pointer(j) - 1)),
                    sqrt(sigma2_phi1_old),
                    TRUE);
  
   /*Decision*/
   double ratio = exp(numer - denom);   
   double acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
    
     delta_trans(j) = delta_trans_old(j);
     delta(j) = delta_old(j);
     z_delta = z_delta_old;
     mu_y_trans = mu_y_trans_old;
     acc = 0;
    
     }

   acctot_delta(j) = acctot_delta(j) + 
                     acc;
   
   }

return Rcpp::List::create(Rcpp::Named("delta") = delta,
                          Rcpp::Named("delta_trans") = delta_trans,
                          Rcpp::Named("z_delta") = z_delta,
                          Rcpp::Named("mu_y_trans") = mu_y_trans,
                          Rcpp::Named("acctot_delta") = acctot_delta);

}
