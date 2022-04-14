#include <RcppArmadillo.h>
#include "decl.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".f_psi_arma")]]
arma::vec f_psi_cpp(arma::vec x, Rcpp::List data, arma::vec psi){
  unsigned int d = data.length();
  arma::vec m = arma::zeros(d);
  arma::vec res = arma::zeros(x.n_elem);
  for(unsigned int b = 0; b < x.n_elem; b++){
    for(unsigned int j = 0; j < d; j++){
      arma::vec v1 = data[j];
      arma::uvec a1 = arma::find(v1 < x(b));
      arma::uvec a2 = arma::find(v1 == x(b));
      m(j) = (a1.n_elem + 0.5 * a2.n_elem) / v1.n_elem;
    }
    res(b) = arma::sum(psi % m);
  }
  return res;
}

double f_psi_cpp2(double x, Rcpp::List data, arma::vec psi){
  unsigned int d = data.length();
  arma::vec m = arma::zeros(d);
  double res;
  for(unsigned int j = 0; j < d; j++){
    arma::vec v1 = data[j];
    arma::uvec a1 = arma::find(v1 < x);
    arma::uvec a2 = arma::find(v1 == x);
    m(j) = (a1.n_elem + 0.5 * a2.n_elem) / v1.n_elem;
  }
  res = arma::sum(psi % m);
  return res;
}



// [[Rcpp::export(".f_theta_arma")]]
arma::vec f_theta_cpp(arma::vec x, Rcpp::List data, arma::vec theta, Rcpp::List psi){
  arma::vec res = arma::zeros(x.n_elem);
  for(unsigned int j = 0; j < x.n_elem; j++){
    arma::vec ab = arma::zeros(data.length());
    for(unsigned int i = 0; i < ab.n_elem; i++){
      Rcpp::List v2 = data[i];
      arma::vec psi_i = psi[i];
      ab(i) = f_psi_cpp2(x(j), v2, psi_i);
    }
    res[j] = arma::sum(theta % ab);
  }
  return res;
}

//[[Rcpp::export(".rel_eff_arma")]]
arma::vec rel_eff_cpp(Rcpp::List data, arma::vec theta, Rcpp::List psi){
  unsigned int d = data.length();
  arma::vec p = arma::zeros(d);
  
  for(unsigned int i = 0; i < d; i++){
    Rcpp::List vi = data[i];
    unsigned int vi_n = vi.length();
    arma::vec mm = arma::zeros(vi_n);
    for(unsigned int j = 0; j < vi_n; j++){
      arma::vec vij = vi[j];
      mm(j) = arma::sum(f_theta_cpp(vij, data, theta, psi)) / vij.n_elem;
    }
    arma::vec psi_i = psi[i];
    p(i) = arma::sum(psi_i % mm);
  }
  return p;
}

