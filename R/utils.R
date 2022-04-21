#' Empirical Distribution function with cluster weights psi
#'
#' @param x A scalar or vector to be checked as F(x)
#' @param i The group/sample index
#' @param data The data, provided as a list of lists
#' @param psi A vector of weights, defaults to unweighted estimator
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return Value(s) of the empirical distribution function 

f_psi = function(x, i, data, psi = NULL, type = NULL){
  #m = numeric(length(data[[i]])) #empty vector for means of clusters
  # If psi is not provided create it
  if(is.null(psi)) {
      if(is.null(type)) psi = weight_fun(data, "unweighted")$psi[[i]]
      psi = weight_fun(data, type)$psi[[i]]
    }#psi = rep(1/length(data[[i]]), length(data[[i]]))
  m = .f_psi_arma(x, data[[i]], psi)
  return(m)
}


#' Empirical Reference Distribution Function with group weights theta and cluster weights psi
#' 
#' @param x A scalar or vector to be checked as F(x)
#' @param data The data, provided as a list of lists
#' @param theta A vector of group weights, defaults to unweighted estimator
#' @param psi A list of vector-weights, defaults to unweighted estimator
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return Value(s) of the empirical distribution function
f_theta = function(x, data, theta = NULL, psi = NULL, type = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  res = .f_theta_arma(x, data, theta, psi)
  return(c(res))
}


#' Computes relative effects with respect to the reference distribution
#' 
#' @param data The data, provided as a list of lists
#' @param theta A vector of group weights, defaults to unweighted estimator
#' @param psi A list of vector-weights, defaults to unweighted estimator
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A vector of the relative effects
rel_eff = function(data, theta = NULL, psi = NULL, type = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  return (c(.rel_eff_arma(data, theta, psi)))
}

#' Inflation-term
#' 
#' @param data The data, provided as a list of lists
#' @param type A string declaring whether a weighted or unweighted estimator will be used, defaults to "unweighted"
#' @return Scalar of inflation
g = function(data, type = "unweighted"){
  ifelse(type == "unweighted", .g_cpp(data, 1), .g_cpp(data, 0))
}

.kappa_r = function(psi, i, j){
  return( 1 - 2*psi[[i]][j] + sum(psi[[i]]^2))
}

#' Computes the Variance-Covariance-Matrix of the relative effects
#' 
#' @param data The data, provided as a list of lists
#' @param theta A vector containing the group weights, defaults to unweighted estimator
#' @param psi A list of vector weights for the clusters, defaults to unweighted estimator
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A Variance-Covariance-Matrix
sigma_est = function(data, theta = NULL, psi = NULL, type = NULL){
  if(is.null(type)) unw = 1
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  if((!is.null(psi)) & !(all(lapply(1:length(psi), FUN = function(x){
    all(psi[[x]] == rep(1 / length(data[[x]]), length(data[[x]])))
  })) == T)){ unw = 0 }
  
  if(!is.null(type)) unw = ifelse(type == "unweighted", 1, 0)
  if(!(type %in% c("unweighted", "weighted"))) unw = 1
  #return( .sigma_est_arma(n, data, theta, psi))
  n = .unsize(data)[[1]]
  return(.sigma_est_arma(n, data, theta, psi, unw))
}



.s_i2 = function(data, i, theta, psi){
  N = .unsize(data)
  n = sum(N[[1]])
  n_i = N[[1]][[i]]
  r_i_bar = 0
  r_i_list = lapply(data[[i]], FUN = function(x){
    n * f_theta(x, data = data, theta = theta, psi = psi) + 0.5}
  ) # Overall Ranks for group i
  M_i  = sum(N[[2]][[i]])
  r_i_int = lapply(data[[i]], FUN = function(x){
    M_i * f_psi(x, i, data, psi = psi[[i]]) + 0.5}
  )
  
  r_i_bar = sum(sapply(r_i_list, mean) * psi[[i]])
  mean_diff = numeric(n_i)
  for(i in 1:n_i){
    mean_diff[i] = mean(r_i_list[[i]] - r_i_int[[i]])
  }
  
  s_i_sq = sum((mean_diff - r_i_bar + (M_i + 1)/2)^2) / (n_i - 1)
  return(s_i_sq)
}

.f_2 = function(data, theta, psi){
  N = .unsize(data)
  n = sum(N[[1]])
  n_i = N[[1]]
  zael = numeric(length(n_i))
  nen = numeric(length(n_i))
  for(i in 1:length(n_i)){
    zael[i] = .s_i2(data, i = i, theta, psi) / (n - n_i[i])
    nen[i]  = (.s_i2(data, i = i, theta, psi) / (n - n_i[i]))^2 / (n_i[i] - 1)
    
  }
  res = sum(zael)^2 / sum(nen)
  return(res)
}

.unsize = function(data){
  N = list()
  N[[1]] = sapply(data, length)
  N[[2]] = lapply(c(1:length(data)), FUN = function(x) sapply(data[[x]], length))
  return(N)
}


.df_sw = function(data, theta, psi, cont){
  n = .unsize(data)[[1]]
  A_i_list = .ai_est_arma(n, data, theta, psi)
  psi_sq = sapply(psi, FUN = function(x) sum(x^2))
  v_list = numeric(nrow(cont))
  for(l in 1:nrow(cont)){
    c_l = cont[l,]
    lambda_list = lapply(A_i_list, FUN = function(x) sapply(x, FUN = function(y) c_l %*% y))
    lambda_bar  = numeric(length(n))
    var_i = numeric(length(n))
    
    for(i in 1:length(n)){
      lambda_bar[i] = psi[[i]] %*% lambda_list[[i]]
      var_i[i] = sum( psi[[i]] * (lambda_list[[i]] - lambda_bar[i])^2 / .kappa_r(psi, i, c(1:length(psi[[i]]))))
    }
    v_list[l] = sum( var_i)^2 / sum(var_i^2 / (n-1) )
    
  }
  return(v_list)
}
