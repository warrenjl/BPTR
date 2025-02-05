#include "RcppArmadillo.h"
#include "BPTR.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List BPTR(int mcmc_samples,
                arma::vec y_trans,
                arma::vec r,
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
                Rcpp::Nullable<double> a_sigma2_zeta0_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_zeta0_prior = R_NilValue,
                Rcpp::Nullable<double> a_sigma2_zeta1_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_zeta1_prior = R_NilValue,
                Rcpp::Nullable<double> a_alpha_prior = R_NilValue,
                Rcpp::Nullable<double> b_alpha_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_eta_prior = R_NilValue,
                Rcpp::Nullable<double> a_sigma2_phi0_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_phi0_prior = R_NilValue, 
                Rcpp::Nullable<double> a_sigma2_phi1_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_phi1_prior = R_NilValue, 
                Rcpp::Nullable<double> a_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> gamma_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> zeta0_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_zeta0_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> zeta1_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_zeta1_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> V_init = R_NilValue,
                Rcpp::Nullable<double> alpha_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> delta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> eta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> phi0_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_phi0_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_phi1_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_epsilon_init = R_NilValue){
  
//Defining Parameters and Quantities of Interest
int n = r.size();
int sum_r = sum(r);
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
arma::vec a1(sum_r); a1.fill(0.00);
if(a1_opt.isNotNull()){
  a1 = Rcpp::as<arma::vec>(a1_opt);
  }
arma::vec c(sum_r); c.fill(0.00);
arma::vec temp(2); temp.fill(0.00);
temp(1) = a2;
arma::vec max_a1_a2(sum_r); max_a1_a2.fill(0.00);
for(int j = 0; j < sum_r; ++j){
  
   temp(0) = a1(j);
   max_a1_a2(j) = max(temp);
   c(j) = t(j, (m(j) - 1)) +
          a0(j) -
          max_a1_a2(j);
  
   }

arma::mat one_star0(sum(m), n); one_star0.fill(0.00);
double val = 0.00;
for(int j = 0; j < sum(m.subvec(0, (r(0) - 1))); ++j){
   one_star0(j,0) = 1;
   }
val = sum(m.subvec(0, (r(0) - 1)));
for(int j = 1; j < n; ++j){
  
   for(int k = val; k < (val + sum(m.subvec(sum(r.subvec(0, (j-1))), (sum(r.subvec(0,j)) - 1)))); ++k){
      one_star0(k,j) = 1;
      }
   val = val +
         sum(m.subvec(sum(r.subvec(0, (j-1))), (sum(r.subvec(0,j)) - 1)));
  
   }

arma::mat one_star0_trans = trans(one_star0);
arma::mat one_star0_t_one_star0 = one_star0_trans*one_star0;
arma::mat one_star1(sum(m), sum_r); one_star1.fill(0.00);
for(int j = 0; j < m(0); ++j){
   one_star1(j,0) = 1;
   }
for(int j = 1; j < sum_r; ++j){
   for(int k = sum(m.subvec(0, (j-1))); k < sum(m.subvec(0,j)); ++k){
      one_star1(k,j) = 1;
      }
   }
arma::mat one_star1_trans = trans(one_star1);
arma::mat one_star1_t_one_star1 = one_star1_trans*one_star1;

arma::mat one_star2(sum_r, n); one_star2.fill(0.00);
for(int j = 0; j < r(0); ++j){
  one_star2(j,0) = 1;
  }
for(int j = 1; j < n; ++j){
   for(int k = sum(r.subvec(0, (j-1))); k < sum(r.subvec(0,j)); ++k){
      one_star2(k,j) = 1;
      }
   }
arma::mat one_star2_trans = trans(one_star2);
arma::mat one_star2_t_one_star2 = one_star2_trans*one_star2;

Rcpp::IntegerVector phi0_pointer(sum_r); phi0_pointer.fill(0);
for(int j = 0; j < r(0); ++j){
   phi0_pointer(j) = 1;
   }
for(int j = 1; j < n; ++j){
   for(int k = sum(r.subvec(0, (j-1))); k < sum(r.subvec(0,j)); ++k){
      phi0_pointer(k) = (j + 1);
      }
   }

arma::mat gamma(p_x, mcmc_samples); gamma.fill(0.00);
arma::mat zeta0(n, mcmc_samples); zeta0.fill(0.00);
arma::vec sigma2_zeta0(mcmc_samples); sigma2_zeta0.fill(0.00);
arma::mat zeta1(sum_r, mcmc_samples); zeta1.fill(0.00);
arma::vec sigma2_zeta1(mcmc_samples); sigma2_zeta1.fill(0.00);
arma::mat theta(d, mcmc_samples); theta.fill(0.00);
theta.row(d-1).fill(1.00);
arma::mat V((d-1), mcmc_samples); V.fill(0.00);
arma::vec alpha(mcmc_samples); alpha.fill(0.00);
arma::mat delta(sum_r, mcmc_samples); delta.fill(0.00);
arma::mat eta(p_z, mcmc_samples); eta.fill(0.00);
arma::mat phi0(n, mcmc_samples); phi0.fill(0.00);
arma::vec sigma2_phi0(mcmc_samples); sigma2_phi0.fill(0.00);
arma::vec sigma2_phi1(mcmc_samples); sigma2_phi1.fill(0.00);
arma::vec sigma2_epsilon(mcmc_samples); sigma2_epsilon.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_gamma = 10000.00;
if(sigma2_gamma_prior.isNotNull()){
  sigma2_gamma = Rcpp::as<double>(sigma2_gamma_prior);
  }

double a_sigma2_zeta0 = 0.01;
if(a_sigma2_zeta0_prior.isNotNull()){
  a_sigma2_zeta0 = Rcpp::as<double>(a_sigma2_zeta0_prior);
  }

double b_sigma2_zeta0 = 0.01;
if(b_sigma2_zeta0_prior.isNotNull()){
  b_sigma2_zeta0 = Rcpp::as<double>(b_sigma2_zeta0_prior);
  }

double a_sigma2_zeta1 = 0.01;
if(a_sigma2_zeta1_prior.isNotNull()){
  a_sigma2_zeta1 = Rcpp::as<double>(a_sigma2_zeta1_prior);
  }

double b_sigma2_zeta1 = 0.01;
if(b_sigma2_zeta1_prior.isNotNull()){
  b_sigma2_zeta1 = Rcpp::as<double>(b_sigma2_zeta1_prior);
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

double a_sigma2_phi0 = 0.01;
if(a_sigma2_phi0_prior.isNotNull()){
  a_sigma2_phi0 = Rcpp::as<double>(a_sigma2_phi0_prior);
  }

double b_sigma2_phi0 = 0.01;
if(b_sigma2_phi0_prior.isNotNull()){
  b_sigma2_phi0 = Rcpp::as<double>(b_sigma2_phi0_prior);
  }

double a_sigma2_phi1 = 0.01;
if(a_sigma2_phi1_prior.isNotNull()){
  a_sigma2_phi1 = Rcpp::as<double>(a_sigma2_phi1_prior);
  }

double b_sigma2_phi1 = 0.01;
if(b_sigma2_phi1_prior.isNotNull()){
  b_sigma2_phi1 = Rcpp::as<double>(b_sigma2_phi1_prior);
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

zeta0.col(0).fill(0.00);
if(zeta0_init.isNotNull()){
  zeta0.col(0) = Rcpp::as<arma::vec>(zeta0_init);
  }

sigma2_zeta0(0) = 0.01;
if(sigma2_zeta0_init.isNotNull()){
  sigma2_zeta0(0) = Rcpp::as<double>(sigma2_zeta0_init);
  }

zeta1.col(0).fill(0.00);
if(zeta1_init.isNotNull()){
  zeta1.col(0) = Rcpp::as<arma::vec>(zeta1_init);
  }

sigma2_zeta1(0) = 0.01;
if(sigma2_zeta1_init.isNotNull()){
  sigma2_zeta1(0) = Rcpp::as<double>(sigma2_zeta1_init);
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

phi0.col(0).fill(0.00);
if(phi0_init.isNotNull()){
  phi0.col(0) = Rcpp::as<arma::vec>(phi0_init);
  }

sigma2_phi0(0) = 0.01;
if(sigma2_phi0_init.isNotNull()){
  sigma2_phi0(0) = Rcpp::as<double>(sigma2_phi0_init);
  }

sigma2_phi1(0) = 0.01;
if(sigma2_phi1_init.isNotNull()){
  sigma2_phi1(0) = Rcpp::as<double>(sigma2_phi1_init);
  }

arma::vec delta_trans(sum_r); delta_trans.fill(0.00);
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
for(int j = 0; j < sum_r; ++j){
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
                       one_star0*zeta0.col(0) +
                       one_star1*zeta1.col(0) +
                       log(z_delta*full_theta);

neg_two_loglike(0) = neg_two_loglike_update(y_trans,
                                            m,
                                            mu_y_trans,
                                            sigma2_epsilon(0));

//Metropolis Settings
arma::vec acctot_V(d-1); acctot_V.fill(1.00);
arma::vec acctot_delta(sum_r); acctot_delta.fill(1.00);

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
   
   //zeta0 Update
   Rcpp::List zeta0_output = zeta0_update(y_trans,
                                          n,
                                          one_star0,
                                          one_star0_trans,
                                          one_star0_t_one_star0,
                                          mu_y_trans,
                                          zeta0.col(j-1),
                                          sigma2_zeta0(j-1),
                                          sigma2_epsilon(j-1));
   
   zeta0.col(j) = Rcpp::as<arma::vec>(zeta0_output[0]);
   mu_y_trans = Rcpp::as<arma::vec>(zeta0_output[1]);
   
   //sigma2_zeta0 Update
   sigma2_zeta0(j) = sigma2_zeta0_update(n,
                                         a_sigma2_zeta0,
                                         b_sigma2_zeta0,
                                         zeta0.col(j));
   
   //zeta1 Update
   Rcpp::List zeta1_output = zeta1_update(y_trans,
                                          sum_r,
                                          one_star1,
                                          one_star1_trans,
                                          one_star1_t_one_star1,
                                          mu_y_trans,
                                          zeta1.col(j-1),
                                          sigma2_zeta1(j-1),
                                          sigma2_epsilon(j-1));
   
   zeta1.col(j) = Rcpp::as<arma::vec>(zeta1_output[0]);
   mu_y_trans = Rcpp::as<arma::vec>(zeta1_output[1]);
   
   //sigma2_zeta1 Update
   sigma2_zeta1(j) = sigma2_zeta1_update(sum_r,
                                         a_sigma2_zeta1,
                                         b_sigma2_zeta1,
                                         zeta1.col(j));
   
   if(d > 1){
     
     //V (theta) Update
     Rcpp::List V_output = V_update(y_trans,
                                    m,
                                    x,
                                    d,
                                    mu_y_trans,
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
                                          sum_r,
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
                                          phi0_pointer,
                                          mu_y_trans,
                                          theta.col(j),
                                          full_theta,
                                          delta.col(j-1),
                                          delta_trans,
                                          z_delta,
                                          eta.col(j-1),
                                          phi0.col(j-1),
                                          sigma2_phi1(j-1),
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
                           one_star2,
                           p_z,
                           sigma2_eta,
                           delta_trans,
                           phi0.col(j-1),
                           sigma2_phi1(j-1));
   
   //phi0 Update
   phi0.col(j) = phi0_update(z,
                             one_star2_trans,
                             one_star2_t_one_star2,
                             n,
                             delta_trans,
                             eta.col(j),
                             sigma2_phi0(j-1),
                             sigma2_phi1(j-1));
   
   //sigma2_phi0 Update
   sigma2_phi0(j) = sigma2_phi0_update(n,
                                       a_sigma2_phi0,
                                       b_sigma2_phi0,
                                       phi0.col(j));
   
   //sigma2_phi1 Update
   sigma2_phi1(j) = sigma2_phi1_update(sum_r,
                                       z,
                                       one_star2,
                                       a_sigma2_phi1,
                                       b_sigma2_phi1,
                                       delta_trans,
                                       eta.col(j),
                                       phi0.col(j));
   
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
                          Rcpp::Named("zeta0") = zeta0,
                          Rcpp::Named("sigma2_zeta0") = sigma2_zeta0,
                          Rcpp::Named("zeta1") = zeta1,
                          Rcpp::Named("sigma2_zeta1") = sigma2_zeta1,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("V") = V,
                          Rcpp::Named("alpha") = alpha,
                          Rcpp::Named("delta") = delta,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("phi0") = phi0,
                          Rcpp::Named("sigma2_phi0") = sigma2_phi0,
                          Rcpp::Named("sigma2_phi1") = sigma2_phi1,
                          Rcpp::Named("sigma2_epsilon") = sigma2_epsilon,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_V") = acctot_V,
                          Rcpp::Named("acctot_delta") = acctot_delta);

}

