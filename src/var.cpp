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
  Rcpp::List A_i_list(d);
  
  for(int i = 0; i < d; i++){
    Rcpp::List data_i = data(i);
    int n_i = n(i);
    Rcpp::List A_ij_list(n_i);
    for(int j = 0; j < n_i; j++){
      arma::vec A_ij = arma::zeros(d);
      arma::vec data_ij = data_i(j);
      for(int h = 0; h < d; h++){
        if(h == i){
          arma::vec ind = arma::regspace(0, d - 1);//
          arma::vec  y = arma::zeros(d);
          for(int s = 0; s < d; s++){
            if(ind(s) != i){
              y(s) = Y_abc(data_ij, data, s) * theta(s);
            }
          }
          A_ij(h) = arma::sum(y);
        } else if(h != i){
          A_ij(h) = (-1) * theta(i) * Y_abc(data_ij, data, h);
          
          
        }  
        
      }
      A_ij_list(j) = A_ij;
    }
    A_i_list(i) = A_ij_list;
  }
  //return A_i_list;
  arma::mat sigma(d, d, arma::fill::zeros);
  for(int i = 0; i < d; i++){
    arma::vec A_ibar = arma::zeros(d);
    Rcpp::List A_i_temp = A_i_list(i);
    arma::vec psi_i = psi(i);
    for(int zz = 0; zz < A_i_temp.length(); zz++){
      arma::vec A_i_z = A_i_temp(zz);
      A_ibar += psi_i(zz) * A_i_z; //A_i_temp(zz);
    }
    Rcpp::List data_i = data(i);
    //     int n_i = data_i.length();
    int n_i = n(i);
    for(int j = 0; j < n_i; j++){
      arma::vec A_ij_temp = A_i_temp(j);
      sigma += ((A_ij_temp - A_ibar) * (A_ij_temp - A_ibar).t()) / (kappa_cpp(psi_i, j)) * pow(psi_i(j), 2);
    }
    //   
  }
  return sigma * g(n);
}


// [[Rcpp::export(".ai_est_arma")]]
Rcpp::List ai_est_cpp(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi){
  int d = n.n_elem;
  Rcpp::List A_i_list(d);
  
  for(int i = 0; i < d; i++){
    Rcpp::List data_i = data(i);
    int n_i = n(i);
    Rcpp::List A_ij_list(n_i);
    for(int j = 0; j < n_i; j++){
      arma::vec A_ij = arma::zeros(d);
      arma::vec data_ij = data_i(j);
      for(int h = 0; h < d; h++){
        if(h == i){
          arma::vec ind = arma::regspace(0, d - 1);//
          arma::vec  y = arma::zeros(d);
          for(int s = 0; s < d; s++){
            if(ind(s) != i){
              y(s) = Y_abc(data_ij, data, s) * theta(s);
            }
          }
          A_ij(h) = arma::sum(y);
        } else if(h != i){
          A_ij(h) = (-1) * theta(i) * Y_abc(data_ij, data, h);
          
          
        }  
        
      }
      A_ij_list(j) = A_ij;
    }
    A_i_list(i) = A_ij_list;
  }
  return A_i_list;
}
