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


.y_abc = function(x_ab, c, data){
  subdata = data[[c]]
  Fcx = numeric(length(x_ab))
  for(k in 1:length(x_ab)){
    x_k = x_ab[k]
    for(i in 1:length(subdata)){
      subsubdata = subdata[[i]]
      Fcx[k] = Fcx[k] + ((length(which(subsubdata < x_k)) + 0.5 * length(which(subsubdata == x_k))) / length(subsubdata))
    }
    Fcx[k] = Fcx[k] / length(subdata)
    
  }
  return(mean(Fcx))
}

.kappa_r = function(psi, i, j){
  return( 1 - 2*psi[[i]][j] + sum(psi[[i]]^2))
}


.sigma_est_r = function(n, data, theta, psi){
  A_i_list = list()
  for(i in 1:length(data)){
    A_ij_list = list()
    
    for(j in 1:length(data[[i]])){
      A_ij = numeric(length(data))
      for(h in 1:length(data)){
        
        if(h == i){
          ### Create vectors of necessary Y's and theta's
          ind = c(1:length(data))[-i]
          y = numeric(length(ind))
          for(s in ind){
            y[s] = .y_abc(data[[i]][[j]], s, data) * theta[s]
          }
          A_ij[h] = sum(y)
          
        } else if (h != i){
          A_ij[h] = -1 * theta[i] * .y_abc(data[[i]][[j]], h, data) 
        }
        
      }
      A_ij_list[[j]] = A_ij
      
    }
    #A_ibar[[i]] = rowMeans(do.call("cbind", A_ij_list)) # Mean of A_ij's <=> Careful: Potentially weighted mean with psi instead
    A_i_list[[i]] = A_ij_list
  }
  
  #####################################
  ##### Berechnung von den Sigmas #####
  
  sigma = matrix(0, length(data), length(data))
  for(i in 1:length(data)){
    A_imat = do.call("cbind", A_i_list[[i]])
    A_ibar = numeric(length(data))
    for(zz in 1:ncol(A_imat)){
      A_ibar = A_ibar + psi[[i]][zz] * A_imat[,zz]
    }
    for(j in 1:length(data[[i]])){
      sigma = sigma + ((A_i_list[[i]][[j]] - A_ibar) %*% t(A_i_list[[i]][[j]] - A_ibar)) /  (.kappa_r(psi, i, j))  * psi[[i]][j]^2 
    }
  }
  return( sigma * g(n))
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
  #return( .sigma_est_arma(n, data, theta, psi))
  return(.sigma_est_r(n, data, theta, psi))
}