#include <RcppArmadillo.h>
#include "decl.h"
// [[Rcpp::depends(RcppArmadillo)]]

// -- Stats --
// [[Rcpp::export(".q_wald_arma")]]
Rcpp::List q_wald_arma(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, arma::mat cmat, int unw){
  arma::mat res;
  arma::vec p = rel_eff_cpp(data, theta, psi);
  arma::mat sigma = sigma_est_cpp(n, data, theta, psi, unw);
  res = p.t() * cmat.t() * arma::pinv(cmat * sigma * cmat.t()) * cmat * p;
  unsigned int df_1 = arma::rank(cmat * sigma);
  return Rcpp::List::create(res, df_1);
}

// [[Rcpp::export(".q_anova_arma")]]
Rcpp::List q_anova_arma(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, arma::mat cmat, int unw){
  arma::mat res;
  arma::vec p = rel_eff_cpp(data, theta, psi);
  arma::mat sigma = sigma_est_cpp(n, data, theta, psi, unw);
  arma::mat M = cmat.t() * arma::pinv(cmat * cmat.t()) * cmat;
  double nen = arma::trace(M * sigma);
  res = p.t() * M * p / nen;
  double df_1 = pow(arma::trace(M * sigma), 2) / arma::trace(M * sigma * M * sigma);
  return Rcpp::List::create(res, df_1, nen);
}

