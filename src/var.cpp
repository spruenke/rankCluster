#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double Y_abc(arma::vec x_ab, Rcpp::List data, int c){
  int m_ab = x_ab.n_elem;
  int d = data.length();
  arma::vec res_1 = arma::zeros(m_ab);
  for(int k = 0; k < m_ab; k++){
    //for(int i = 0; i < d; i++){
    Rcpp::List sub_data = data[c];
    int n_i = sub_data.length();
    arma::vec m_ij(n_i);
    arma::vec vnn(n_i);
    for(int j = 0; j < n_i; j++){
      arma::vec v1 = sub_data[j];
      arma::uvec a1 = arma::find(v1 < x_ab(k));
      arma::uvec a2 = arma::find(v1 == x_ab(k));
      m_ij(j) = (a1.n_elem + 0.5 * a2.n_elem); // v1.n_elem;
      vnn(j)  = v1.n_elem;
    }
    res_1(k) = arma::sum(m_ij) / arma::sum(vnn);
    //}
  }
  return arma::mean(res_1);
}

double kappa_cpp(arma::vec psi, int j){
  double res = 1 - 2 * psi[j] + arma::dot(psi, psi);
  return res;
}

double g(arma::vec n){
  return arma::sum(n);
}

// [[Rcpp::export(".sigma_est_arma")]]
arma::mat sigma_est_cpp(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi){
  int d = n.n_elem;
  arma::vec ind = arma::regspace(0, (d-1));
  Rcpp::List A(d);
  Rcpp::List A_bar(d);
  Rcpp::List sigma_list(d);
  for(int i = 0; i < d; i++){
    int nii = n(i);
    Rcpp::List sublist(nii);
    Rcpp::List subdat = data(i);
    arma::vec psi_i = psi(i);
    // for(int j = 0; j < nii; j++){
    //   arma::vec subvec = arma::zeros(d);
    //   arma::vec subdat_j = subdat[j];
    //   for(int s = 0; s < d; s++){
    //     if(s == i){
    //       arma::vec ind_new = ind(arma::find(ind != s));
    //       for(unsigned int hh = 0; hh < ind_new.n_elem; hh++){
    //         unsigned int h = ind_new(hh);
    //         Rcpp::List dat_h = data[h];
    //         subvec(s) += theta(h) * arma::sum(f_psi_cpp(subdat_j, dat_h, psi_i)) / subdat_j.n_elem;
    //       }
    //     } else {
    //       subvec(s) = (-1) * theta(s) * arma::sum(f_psi_cpp(subdat_j, data[s], psi_i)) / subdat_j.n_elem;
    //     }
    //   }
    //   sublist(j) = subvec;
    // }
    for(int j = 0; j < nii; j++){
      arma::vec subvec = arma::zeros(d);
      arma::vec subdat_j = subdat[j];
      for(int s = 0; s < d; s++){
        if(s == i){
          arma::uvec ind_new = arma::find(ind != s);
          for(unsigned int hh = 0; hh < ind_new.n_elem; hh++){
            unsigned int h = ind_new(hh);
            Rcpp::List dat_h = data[h];
            subvec(s) += theta(h) * Y_abc(subdat_j, data, h);
          }
        } else {
          subvec(s) = (-1) * theta(s) * Y_abc(subdat_j, data, s);
        }
      }
      sublist(j) = subvec;
    }
    A(i) = sublist;
    arma::vec tmp = arma::zeros(d);
    for(int t = 0; t < d; t++){
      for(int j = 0; j < sublist.length(); j++){
        arma::vec A_vec = sublist[j];
        tmp(t) += A_vec(t) * psi_i(j);
      }
    }
    A_bar(i) = tmp;
    Rcpp::List sublist2(nii);
    
    for(int j = 0; j < sublist.length(); j++){
      arma::vec A_ij = sublist[j];
      arma::vec y_lhs = A_ij - tmp;
      sublist2(j) = y_lhs * y_lhs.t(); // sigma_list[[i]][[j]]
    }
    sigma_list(i) = sublist2;
    
  }
  
  // Create sigma matrix
  arma::mat sigma(d, d, arma::fill::zeros);
  for(int i = 0; i < d; i++){
    Rcpp::List siglist = sigma_list[i];
    arma::vec psi_i = psi[i];
    
    for(int j = 0; j < siglist.length(); j++){
      arma::mat subsig = siglist[j];
      sigma += subsig * (pow(psi_i[j], 2)) / kappa_cpp(psi_i, j);
    }
  }
  
  return g(n) * sigma;
  // return sigma_list;
}