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
#' @param n a vector containing the sample sizes for each group
#' @return Scalar of inflation
g = function(n){
  sum(n)
}



#' Computes the Variance-Covariance-Matrix of the relative effects
#' 
#' @param n A vector of sample sizes for each group
#' @param data The data, provided as a list of lists
#' @param theta A vector containing the group weights, defaults to unweighted estimator
#' @param psi A list of vector weights for the clusters, defaults to unweighted estimator
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A Variance-Covariance-Matrix
sigma_est = function(n, data, theta = NULL, psi = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  return( .sigma_est_arma(n, data, theta, psi))
}