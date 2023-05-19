#include "RcppArmadillo.h"
#include "BETR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List BETR(int mcmc_samples,
                arma::vec y_trans,
                arma::vec m,
                arma::mat t,
                arma::mat x,
                arma::mat z,
                arma::vec a0,
                double a2,
                int d,
                arma::vec metrop_V,
                arma::vec metrop_var_delta,
                Rcpp::Nullable<Rcpp::NumericVector> a1_opt = R_NilValue,
                Rcpp::Nullable<double> sigma2_gamma_prior = R_NilValue,
                Rcpp::Nullable<double> a_sigma2_zeta_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_zeta_prior = R_NilValue,
                Rcpp::Nullable<double> a_alpha_prior = R_NilValue,
                Rcpp::Nullable<double> b_alpha_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_eta_prior = R_NilValue,
                Rcpp::Nullable<double> a_sigma2_phi_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_phi_prior = R_NilValue, 
                Rcpp::Nullable<double> a_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> gamma_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> zeta_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_zeta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> V_init = R_NilValue,
                Rcpp::Nullable<double> alpha_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> delta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> eta_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_phi_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_epsilon_init = R_NilValue){
  
//Defining Parameters and Quantities of Interest
int n = t.n_rows;
arma::mat x_trans = trans(x);
arma::mat xtx = x_trans*x;
arma::mat z_trans = trans(z);
arma::mat ztz = z_trans*z;
int p_x = x.n_cols;
int p_z = z.n_cols;
double theta_0 = 0.00;
arma::vec choose_vec(d + 1); choose_vec.fill(0.00);
for(int j = 0; j < (d + 1); ++j){
   choose_vec(j) = R::choose(d, j);
   }
arma::vec a1(n); a1.fill(0.00);
if(a1_opt.isNotNull()){
  a1 = Rcpp::as<arma::vec>(a1_opt);
  }
arma::vec c(n); c.fill(0.00);
arma::vec temp(2); temp.fill(0.00);
temp(1) = a2;
arma::vec max_a1_a2(n); max_a1_a2.fill(0.00);
for(int j = 0; j < n; ++j){
  
   temp(0) = a1(j);
   max_a1_a2(j) = max(temp);
   c(j) = t(j, (m(j) - 1)) +
          a0(j) -
          max_a1_a2(j);
  
   }
arma::mat one_star(sum(m), n); one_star.fill(0.00);
for(int j = 0; j < m(0); ++j){
  one_star(j,0) = 1;
  }
for(int j = 1; j < n; ++j){
   for(int k = sum(m.subvec(0, (j-1))); k < sum(m.subvec(0,j)); ++k){
      one_star(k,j) = 1;
      }
   }
arma::mat one_star_trans = trans(one_star);
arma::mat one_star_t_one_star = one_star_trans*one_star;

arma::mat gamma(p_x, mcmc_samples); gamma.fill(0.00);
arma::mat zeta(n, mcmc_samples); zeta.fill(0.00);
arma::vec sigma2_zeta(mcmc_samples); sigma2_zeta.fill(0.00);
arma::mat theta(d, mcmc_samples); theta.fill(0.00);
theta.row(d-1).fill(1.00);
arma::mat V((d-1), mcmc_samples); V.fill(0.00);
arma::vec alpha(mcmc_samples); alpha.fill(0.00);
arma::mat delta(n, mcmc_samples); delta.fill(0.00);
arma::mat eta(p_z, mcmc_samples); eta.fill(0.00);
arma::vec sigma2_phi(mcmc_samples); sigma2_phi.fill(0.00);
arma::vec sigma2_epsilon(mcmc_samples); sigma2_epsilon.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_gamma = 10000.00;
if(sigma2_gamma_prior.isNotNull()){
  sigma2_gamma = Rcpp::as<double>(sigma2_gamma_prior);
  }

double a_sigma2_zeta = 0.01;
if(a_sigma2_zeta_prior.isNotNull()){
  a_sigma2_zeta = Rcpp::as<double>(a_sigma2_zeta_prior);
  }

double b_sigma2_zeta = 0.01;
if(b_sigma2_zeta_prior.isNotNull()){
  b_sigma2_zeta = Rcpp::as<double>(b_sigma2_zeta_prior);
  }

double a_alpha = 0.01;
if(a_alpha_prior.isNotNull()){
  a_alpha = Rcpp::as<double>(a_alpha_prior);
  }

double b_alpha = 0.01;
if(b_alpha_prior.isNotNull()){
  b_alpha = Rcpp::as<double>(b_alpha_prior);
  }

double sigma2_eta = 10000.00;
if(sigma2_eta_prior.isNotNull()){
  sigma2_eta = Rcpp::as<double>(sigma2_eta_prior);
  }

double a_sigma2_phi = 0.01;
if(a_sigma2_phi_prior.isNotNull()){
  a_sigma2_phi = Rcpp::as<double>(a_sigma2_phi_prior);
  }

double b_sigma2_phi = 0.01;
if(b_sigma2_phi_prior.isNotNull()){
  b_sigma2_phi = Rcpp::as<double>(b_sigma2_phi_prior);
  }

double a_sigma2_epsilon = 0.01;
if(a_sigma2_epsilon_prior.isNotNull()){
  a_sigma2_epsilon = Rcpp::as<double>(a_sigma2_epsilon_prior);
  }

double b_sigma2_epsilon = 0.01;
if(b_sigma2_epsilon_prior.isNotNull()){
  b_sigma2_epsilon = Rcpp::as<double>(b_sigma2_epsilon_prior);
  }

//Initial Values
gamma.col(0).fill(0.00);
if(gamma_init.isNotNull()){
  gamma.col(0) = Rcpp::as<arma::vec>(gamma_init);
  }

zeta.col(0).fill(0.00);
if(zeta_init.isNotNull()){
  zeta.col(0) = Rcpp::as<arma::vec>(zeta_init);
  }

sigma2_zeta(0) = 0.01;
if(sigma2_zeta_init.isNotNull()){
  sigma2_zeta(0) = Rcpp::as<double>(sigma2_zeta_init);
  }

arma::vec psi_less(d-1); psi_less.fill(1.00);
arma::vec psi(d); psi.fill(0.00);
if(d > 1){
  
  alpha(0) = 1.00;
  if(alpha_init.isNotNull()){
    alpha(0) = Rcpp::as<double>(alpha_init);
    }

  V.col(0).fill(0.50);
  if(V_init.isNotNull()){
    V.col(0) = Rcpp::as<arma::vec>(V_init);
    }

  arma::vec psi_less(d-1); psi_less.fill(1.00); 
  if(d > 2){
    psi_less.subvec(1, (d-2)) = arma::cumprod(1.00 - V.col(0).subvec(0, (d-3)));
    }
  psi_less = psi_less%V.col(0);
  arma::vec psi(d); psi.fill(0.00);
  psi.subvec(0, (d-2)) = psi_less;
  psi(d-1) = 1.00 - 
             sum(psi_less);

  theta(0,0) = psi(0);
  for(int j = 1; j < d; ++j){
     theta(j,0) = theta((j-1), 0) +
                  psi(j);
     }
  
  }
if(theta_init.isNotNull()){
  
  theta.col(0) = Rcpp::as<arma::vec>(theta_init);
  psi(0) = theta(0,0);
  for(int j = 1; j < d; ++j){
     psi(j) = theta(j,0) -
              theta((j-1), 0);
     }
  
  }

eta.col(0).fill(0.00);
if(eta_init.isNotNull()){
  eta.col(0) = Rcpp::as<arma::vec>(eta_init);
  }

sigma2_phi(0) = 0.01;
if(sigma2_phi_init.isNotNull()){
  sigma2_phi(0) = Rcpp::as<double>(sigma2_phi_init);
  }

arma::vec delta_trans(n); delta_trans.fill(0.00);
delta.col(0) = (a0 - max_a1_a2)/(1.00 + exp(delta_trans));  
if(delta_init.isNotNull()){
  
  delta.col(0) = Rcpp::as<arma::vec>(delta_init);
  delta_trans = log((a0 - max_a1_a2)/delta.col(0) - 1.00);
  
  }
Rcpp::IntegerVector power1 = seq(0, d);
Rcpp::IntegerVector power2 = d - power1;
arma::mat z_delta(sum(m), (d+1)); z_delta.fill(0.00);
arma::vec vec1(d+1); vec1.fill(0.00);
arma::vec vec2(d+1); vec2.fill(0.00);
int counter = -1;
for(int j = 0; j < n; ++j){
   for(int k = 0; k < m(j); ++k){ 
    
      counter = counter + 
                1;
     
      for(int l = 0; l < (d+1); ++l){
        
         vec1(l) = pow((t(j,k) + delta(j,0))/max(c), power1(l));
         vec2(l) = pow((1.00 - (t(j,k) + delta(j,0))/max(c)), power2(l));
       
         }
      z_delta.row(counter) = trans(choose_vec%vec1%vec2);
    
      } 
   } 

sigma2_epsilon(0) = 0.01;
if(sigma2_epsilon_init.isNotNull()){
  sigma2_epsilon(0) = Rcpp::as<double>(sigma2_epsilon_init);
  }

arma::vec full_theta(d+1); full_theta.fill(0.00);
full_theta(0) = theta_0;
full_theta.subvec(1, d) = theta.col(0);
arma::vec mu_y_trans = x*gamma.col(0) + 
                       one_star*zeta.col(0) + 
                       log(z_delta*full_theta);

neg_two_loglike(0) = neg_two_loglike_update(y_trans,
                                            m,
                                            mu_y_trans,
                                            sigma2_epsilon(0));

//Metropolis Settings
arma::vec acctot_V(d-1); acctot_V.fill(1.00);
arma::vec acctot_delta(n); acctot_delta.fill(1.00);

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){
  
   //gamma Update
   Rcpp::List gamma_output = gamma_update(y_trans,
                                          x,
                                          x_trans,
                                          xtx,
                                          p_x,
                                          sigma2_gamma,
                                          mu_y_trans,
                                          gamma.col(j-1),
                                          sigma2_epsilon(j-1));
  
   gamma.col(j) = Rcpp::as<arma::vec>(gamma_output[0]);
   mu_y_trans = Rcpp::as<arma::vec>(gamma_output[1]);
   
   //zeta Update
   Rcpp::List zeta_output = zeta_update(y_trans,
                                        n,
                                        one_star,
                                        one_star_trans,
                                        one_star_t_one_star,
                                        mu_y_trans,
                                        zeta.col(j-1),
                                        sigma2_zeta(j-1),
                                        sigma2_epsilon(j-1));
   
   zeta.col(j) = Rcpp::as<arma::vec>(zeta_output[0]);
   mu_y_trans = Rcpp::as<arma::vec>(zeta_output[1]);
   
   //sigma2_zeta Update
   sigma2_zeta(j) = sigma2_zeta_update(n,
                                       a_sigma2_zeta,
                                       b_sigma2_zeta,
                                       zeta.col(j));
   
   if(d > 1){
     
     //V (theta) Update
     Rcpp::List V_output = V_update(y_trans,
                                    m,
                                    x,
                                    one_star,
                                    d,
                                    mu_y_trans,
                                    gamma.col(j),
                                    zeta.col(j),
                                    V.col(j-1),
                                    theta.col(j-1),
                                    psi_less,
                                    psi,
                                    full_theta,
                                    alpha(j-1),
                                    z_delta,
                                    sigma2_epsilon(j-1),
                                    metrop_V,
                                    acctot_V);
   
     V.col(j) = Rcpp::as<arma::vec>(V_output[0]);
     theta.col(j) = Rcpp::as<arma::vec>(V_output[1]);
     psi_less = Rcpp::as<arma::vec>(V_output[2]);
     psi = Rcpp::as<arma::vec>(V_output[3]);
     full_theta = Rcpp::as<arma::vec>(V_output[4]);
     mu_y_trans = Rcpp::as<arma::vec>(V_output[5]);
     acctot_V = Rcpp::as<arma::vec>(V_output[6]);
   
     //alpha Update
     alpha(j) = alpha_update(d,
                             a_alpha,
                             b_alpha,
                             V.col(j));
     
     }
   
   //delta Update 
   Rcpp::List delta_output = delta_update(y_trans,
                                          n,
                                          m,
                                          t,
                                          z,
                                          c,
                                          d,
                                          a0,
                                          max_a1_a2,
                                          choose_vec,
                                          power1,
                                          power2,
                                          mu_y_trans,
                                          theta.col(j),
                                          full_theta,
                                          delta.col(j-1),
                                          delta_trans,
                                          z_delta,
                                          eta.col(j-1),
                                          sigma2_phi(j-1),
                                          sigma2_epsilon(j-1),
                                          metrop_var_delta,
                                          acctot_delta);
   
   delta.col(j) = Rcpp::as<arma::vec>(delta_output[0]);
   delta_trans = Rcpp::as<arma::vec>(delta_output[1]);
   z_delta = Rcpp::as<arma::mat>(delta_output[2]);
   mu_y_trans = Rcpp::as<arma::vec>(delta_output[3]);
   acctot_delta = Rcpp::as<arma::vec>(delta_output[4]);
   
   //eta Update
   eta.col(j) = eta_update(z_trans,
                           ztz,
                           p_z,
                           sigma2_eta,
                           delta_trans,
                           sigma2_phi(j-1));
   
   //sigma2_phi Update
   sigma2_phi(j) = sigma2_phi_update(n,
                                     z,
                                     a_sigma2_phi,
                                     b_sigma2_phi,
                                     delta_trans,
                                     eta.col(j));
   
   //sigma2_epsilon Update
   sigma2_epsilon(j) = sigma2_epsilon_update(y_trans,
                                             m,
                                             mu_y_trans,
                                             a_sigma2_epsilon,
                                             b_sigma2_epsilon);
   
   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(y_trans,
                                               m,
                                               mu_y_trans,
                                               sigma2_epsilon(j));
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     
     if(d > 1){
       
       double accrate_V_min = round(100*(min(acctot_V)/(double)j));
       Rcpp::Rcout << "V Acceptance (min): " << accrate_V_min << "%" << std::endl;
     
       double accrate_V_max = round(100*(max(acctot_V)/(double)j));
       Rcpp::Rcout << "V Acceptance (max): " << accrate_V_max << "%" << std::endl;
     
       }
     
     double accrate_delta_min = round(100*(min(acctot_delta)/(double)j));
     Rcpp::Rcout << "delta Acceptance (min): " << accrate_delta_min << "%" << std::endl;
     
     double accrate_delta_max = round(100*(max(acctot_delta)/(double)j));
     Rcpp::Rcout << "delta Acceptance (max): " << accrate_delta_max << "%" << std::endl;
     
     
     Rcpp::Rcout << "***************************" << std::endl;
    
     }
  
   }
                                  
return Rcpp::List::create(Rcpp::Named("gamma") = gamma,
                          Rcpp::Named("zeta") = zeta,
                          Rcpp::Named("sigma2_zeta") = sigma2_zeta,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("V") = V,
                          Rcpp::Named("alpha") = alpha,
                          Rcpp::Named("delta") = delta,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("sigma2_phi") = sigma2_phi,
                          Rcpp::Named("sigma2_epsilon") = sigma2_epsilon,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_V") = acctot_V,
                          Rcpp::Named("acctot_delta") = acctot_delta);

}

