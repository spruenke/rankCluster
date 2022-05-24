#ifndef decl_h
#define decl_h

#include <Rcpp.h>
#include <RcppArmadillo.h>

// --- Variance

double Y_abc(arma::vec x_ab, Rcpp::List data, int c);
double kappa_cpp(arma::vec psi, int j);
int g(Rcpp::List data, int unw);

arma::mat sigma_est_cpp(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, int unw);

Rcpp::List ai_est_cpp(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi);

// --- Utils

arma::vec f_psi_cpp(arma::vec x, Rcpp::List data, arma::vec psi);

double f_psi_cpp2(double x, Rcpp::List data, arma::vec psi);

arma::vec f_theta_cpp(arma::vec x, Rcpp::List data, arma::vec theta, Rcpp::List psi);

arma::vec rel_eff_cpp(Rcpp::List data, arma::vec theta, Rcpp::List psi);

// --- Tests

Rcpp::List q_wald_arma(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, arma::mat cmat, int unw);

Rcpp::List q_anova_arma(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, arma::mat cmat, int unw);


#endif