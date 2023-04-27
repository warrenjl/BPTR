#include "RcppArmadillo.h"
#include "BETR2.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update(arma::vec delta,
                    arma::mat z,
                    arma::vec delta_trans,
                    arma::vec metrop_sd_delta,
                    arma::vec y_trans,
                    arma::vec mu_y_trans,
                    arma::mat z_delta,
                    arma::vec theta,
                    arma::vec eta,
                    arma::vec a0,
                    double sigma2_epsilon,
                    double sigma2_phi,
                    double alpha,
                    double theta0,
                    int a2,
                    arma::vec m,
                    arma::vec acctot_V,
                    int n){
  
  int counter0; 
  arma::vec dens(n); dens.fill(0.00);
  arma::vec mean_old;
  int etaS = eta.size();
  for(int j = 0; j < (n-1); ++j){
    if(j==1){
      counter0 = 0; 
    }
    if(j > 1){
      counter0 = arma::sum(m(1,(j-1)));
    }

    arma::vec delta_trans_old = delta_trans;
    arma::vec z_delta_old = z_delta; 
    arma::vec mu_y_trans_old = mu_y_trans; 
    
    double second = 0.00;
    double first = 0.00;
    double ratio = 0.00;
    
    /*Second*/
    for(int i = (counter0 + 1); i<arma::sum(m(1,j)); ++i) {
      dens = R::dnorm(y_trans[i],
                      mu_y_trans_old(i),
                      sqrt(sigma2_epsilon), 
                      TRUE); 
    }
    
    //uvec IDX = as<uvec>(seq(0, z.n_cols, k));
    arma::vec zcol = z*eta;
    
    second = arma::sum(dens) +
             R::dnorm(delta_trans_old[j],
                      zcol[j],
                      sqrt(sigma2_phi),
                      TRUE);                                          
    double delta_double = R::rnorm(delta_trans_old[j],
                                  metrop_sd_delta[j]);
    arma::vec delta_vec = as<arma::vec>(wrap(delta_double));
    delta_trans = delta_vec; 
    delta.col(j) = (a0[j] - a2)/(1.00 + exp(delta_trans[j]));  
    int counter1 = counter0; 
  
    for(int k=1; k<m[j]; ++k){ 
      counter1 = counter1 + 1;
      z_delta.row(counter1) = choose(d, c(0:d))*
          (((t[[j]][k] + delta[i,j])/max(c))^c(0:d))*
          ((1.00 - ((t[[j]][k] + delta[i,j])/max(c)))^(d - c(0:d)))
      
    }
    
    mu_y_trans[(counter0 + 1):sum(m[1:j])]<-(mu_y_trans[(counter0 + 1):sum(m[1:j])] - log(z_delta_old[c((counter0 + 1):sum(m[1:j])),]%*%c(theta0, theta[i,]))) +
      log(z_delta[c((counter0 + 1):sum(m[1:j])),]%*%c(theta0, theta[i,]))     
      
      numer<-sum(dnorm(x = y_trans[(counter0 + 1):sum(m[1:j])],
                       mean = mu_y_trans[(counter0 + 1):sum(m[1:j])],
                                        sd = sqrt(sigma2_epsilon[i-1]), 
                                        log = TRUE)) +
                                          dnorm(x = delta_trans[j],
                                                mean = (z[j,]%*%eta[(i-1),]),
                                                sd = sqrt(sigma2_phi[i-1]),
                                                log = TRUE)
      
      ratio<-exp(numer - denom)
      uni_draw<-runif(n = 1,
                      min = 0.00,
                      max = 1.00)                                         		    		     
      acc<-1
    if(ratio < uni_draw){
      
      acc<-0
      delta_trans[j]<-delta_trans_old
      delta[i,j]<-delta[(i-1), j]
      z_delta<-z_delta_old
      mu_y_trans<-mu_y_trans_old
      
    }
    
    acctot_delta[j]<-acctot_delta[j] +
      acc
      
  }
  
  
  
  
  
  
  arma::rowvec one(1); one.fill(1.00);
  arma::rowvec psi(d); //psi.fill(0.00);
  arma::rowvec stick_to_right(d);
  arma::rowvec v_use(d);
  arma::rowvec stick_to_right_temp(d);
  arma::vec theta_old_use(d+1); theta_old_use.fill(0.00);
  arma::vec theta_use(d+1); theta_use.fill(0.00);
  int n = y_trans.size();
  int V_max = V.size();
  arma::vec dens(n); dens.fill(0.00);
  int acc = 0; 
  
  
  for(int j = 0; j < (d-1); ++j){
    arma::vec mu_y_trans_old = mu_y_trans;
    //arma::vec psi_old = psi; 
    arma::vec theta_old = theta; 
    arma::vec V_old = V; 
    
    double second = 0.00;
    double first = 0.00;
    double ratio = 0.00;
    
    /*Second*/
    for(int k = 0; k < n; ++k){
      dens(k) = R::dnorm(y_trans(k),
           mu_y_trans_old(k),
           sqrt(sigma2_epsilon),
           TRUE);
    }
    
    second = arma::sum(dens) +
      R::dbeta(V(j),
               1.00,
               alpha,
               TRUE); 
    
    V(j) = R::runif((V(j) - metrop_V(j)),
      (V(j) + metrop_V(j)));
    
    if(V(j) < 0.00){
      V(j) = abs(V(j));
    }
    if(V(j) > 1.00) {
      V(j) = 2.00 - V(j);
    }
    
    arma::rowvec v_use = as<arma::rowvec>(wrap(V));
    arma::rowvec stick_to_right_temp = arma::cumprod(1 - v_use);
    arma::rowvec subset = stick_to_right_temp.subvec(0, (V_max-2)); //ask Josh
    arma::rowvec stick_to_right = arma::join_rows(one, subset);
    arma::rowvec psi1 = stick_to_right%v_use;
    //psi.row(d-1) = 1.00-arma::sum(psi); //ask Josh
    double psi2 =  1.00-arma::sum(psi1); 
    arma::rowvec psi2_use = as<arma::rowvec>(wrap(psi2));
    arma::rowvec psi = arma::join_rows(psi1, psi2_use); //ask Josh
    theta(1) = psi(1);
    for(int k = 2; k < d; ++k){
      theta(k) = theta(k-1) + psi(k);
    }
    
    theta_old_use(1) = theta0;  
    int l=1; 
    for(int k = 2; k < (d+1); ++k){
      theta_old_use(k) = theta_old(l);
      l = l+1;
    }
    
    theta_use(1) = theta0;  
    int m=1; 
    for(int k = 2; k < (d+1); ++k){
      theta_use(k) = theta(m);
      m = m+1;
    }
    mu_y_trans = (mu_y_trans - log(z_delta*theta_old_use)) +
      log(z_delta*theta_use);
    
    /*First*/
    for(int k = 0; k < n; ++k){
      dens(k) = R::dnorm(y_trans(k),
           mu_y_trans(k),
           sqrt(sigma2_epsilon),
           TRUE);
    }
    
    first = arma::sum(dens) +
      R::dbeta(V(j),
               1.00,
               alpha,
               TRUE); 
    
    /*Decision*/
    ratio = exp(first - second);   
    acc = 1;
    if(ratio < R::runif(0.00, 1.00)){
      V(j) = V_old(j);
      //psi(j) = psi_old(j);  
      theta(j) = theta_old(j);
      mu_y_trans = mu_y_trans_old; 
      acc = 0;
    }
    acctot_V(j) = acctot_V(j) + acc;   
    
  }
  
  return Rcpp::List::create(Rcpp::Named("V")=V,
                            Rcpp::Named("theta")=theta,
                            Rcpp::Named("psi")=psi,
                            Rcpp::Named("acctot_V")=acctot_V);
}