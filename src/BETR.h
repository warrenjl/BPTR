#ifndef __BETR__
#define __BETR__

Rcpp::List gamma_update(arma::vec y_trans,
                        arma::mat x,
                        arma::mat x_trans,
                        arma::mat xtx,
                        int p_x,
                        double sigma2_gamma,
                        arma::vec mu_y_trans,
                        arma::vec gamma_old,
                        double sigma2_epsilon_old);

Rcpp::List zeta0_update(arma::vec y_trans,
                        int n,
                        arma::mat one_star0,
                        arma::mat one_star0_trans,
                        arma::mat one_star0_t_one_star0,
                        arma::vec mu_y_trans,
                        arma::vec zeta0_old,
                        double sigma2_zeta0_old,
                        double sigma2_epsilon_old);

double sigma2_zeta0_update(int n,
                           double a_sigma2_zeta0,
                           double b_sigma2_zeta0,
                           arma::vec zeta0);

Rcpp::List zeta1_update(arma::vec y_trans,
                        int sum_r,
                        arma::mat one_star1,
                        arma::mat one_star1_trans,
                        arma::mat one_star1_t_one_star1,
                        arma::vec mu_y_trans,
                        arma::vec zeta1_old,
                        double sigma2_zeta1_old,
                        double sigma2_epsilon_old);

double sigma2_zeta1_update(int sum_r,
                           double a_sigma2_zeta1,
                           double b_sigma2_zeta1,
                           arma::vec zeta1);

Rcpp::List V_update(arma::vec y_trans,
                    arma::vec m,
                    arma::mat x,
                    int d,
                    arma::vec mu_y_trans_old,
                    arma::vec V_old,
                    arma::vec theta_old,
                    arma::vec psi_less_old,
                    arma::vec psi_old,
                    arma::vec full_theta_old,
                    double alpha_old,
                    arma::mat z_delta_old,
                    double sigma2_epsilon_old,
                    arma::vec metrop_V,
                    arma::vec acctot_V);

double alpha_update(int d,
                    double a_alpha,
                    double b_alpha,
                    arma::vec V);

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
                        arma::vec acctot_delta);

arma::vec eta_update(arma::mat z_trans,
                     arma::mat ztz,
                     arma::mat one_star2,
                     int p_z,
                     double sigma2_eta,
                     arma::vec delta_trans,
                     arma::vec phi0_old,
                     double sigma2_phi1_old);

arma::vec phi0_update(arma::mat z,
                      arma::mat one_star2_trans,
                      arma::mat one_star2_t_one_star2,
                      int n,
                      arma::vec delta_trans,
                      arma::vec eta,
                      double sigma2_phi0_old,
                      double sigma2_phi1_old);

double sigma2_phi0_update(int n,
                          double a_sigma2_phi0,
                          double b_sigma2_phi0,
                          arma::vec phi0);

double sigma2_phi1_update(int sum_r,
                          arma::mat z,
                          arma::mat one_star2,
                          double a_sigma2_phi1,
                          double b_sigma2_phi1,
                          arma::vec delta_trans,
                          arma::vec eta,
                          arma::vec phi0);

double sigma2_epsilon_update(arma::vec y_trans,
                             arma::vec m,
                             arma::vec mu_y_trans,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon);

double neg_two_loglike_update(arma::vec y_trans,
                              arma::vec m,
                              arma::vec mu_y_trans,
                              double sigma2_epsilon);

Rcpp::List BETR(int mcmc_samples,
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
                Rcpp::Nullable<Rcpp::NumericVector> a1_opt,
                Rcpp::Nullable<double> sigma2_gamma_prior,
                Rcpp::Nullable<double> a_sigma2_zeta0_prior,
                Rcpp::Nullable<double> b_sigma2_zeta0_prior,
                Rcpp::Nullable<double> a_sigma2_zeta1_prior,
                Rcpp::Nullable<double> b_sigma2_zeta1_prior,
                Rcpp::Nullable<double> a_alpha_prior,
                Rcpp::Nullable<double> b_alpha_prior,
                Rcpp::Nullable<double> sigma2_eta_prior,
                Rcpp::Nullable<double> a_sigma2_phi0_prior,
                Rcpp::Nullable<double> b_sigma2_phi0_prior, 
                Rcpp::Nullable<double> a_sigma2_phi1_prior,
                Rcpp::Nullable<double> b_sigma2_phi1_prior, 
                Rcpp::Nullable<double> a_sigma2_epsilon_prior,
                Rcpp::Nullable<double> b_sigma2_epsilon_prior,
                Rcpp::Nullable<Rcpp::NumericVector> gamma_init,
                Rcpp::Nullable<Rcpp::NumericVector> zeta0_init,
                Rcpp::Nullable<double> sigma2_zeta0_init,
                Rcpp::Nullable<Rcpp::NumericVector> zeta1_init,
                Rcpp::Nullable<double> sigma2_zeta1_init,
                Rcpp::Nullable<Rcpp::NumericVector> theta_init,
                Rcpp::Nullable<Rcpp::NumericVector> V_init,
                Rcpp::Nullable<double> alpha_init,
                Rcpp::Nullable<Rcpp::NumericVector> delta_init,
                Rcpp::Nullable<Rcpp::NumericVector> eta_init,
                Rcpp::Nullable<Rcpp::NumericVector> phi0_init,
                Rcpp::Nullable<double> sigma2_phi0_init,
                Rcpp::Nullable<double> sigma2_phi1_init,
                Rcpp::Nullable<double> sigma2_epsilon_init); 

#endif // __BETR__
