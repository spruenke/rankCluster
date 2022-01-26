#' Computes the Wald-Test for clustered data
#' 
#' @param n A vector containing the sample sizes
#' @param data The data, provided as a list of lists
#' @param cont A contrast matrix for the relative effects
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-value and the test decision
q_wald = function(n, data, cont, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  phat = rel_eff(data, theta, psi)
  sigma = sigma_est(n, data, theta, psi)
  stat = (t(phat) %*% t(c_mat) %*% MASS::ginv(c_mat %*% sigma %*% t(c_mat)) %*% c_mat %*% phat)  * g(n)
  df = Matrix::rankMatrix(c_mat %*% sigma)
  pv   = 1 - pchisq(stat, df)
  dec  = pv < alpha
  return(list(Statistic = stat, df = df, p.value = pv, reject = dec))
}


#' Computes the ANOVA-Test for clustered data
#' 
#' @param n A vector containing the sample sizes
#' @param data The data, provided as a list of lists
#' @param cont A contrast matrix for the relative effects
#' @param f_2 A value for the second degree of freedom
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-value and the test decision
q_anova = function(n, data, cont, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  phat = rankCluster::rel_eff(data, theta, psi)
  sigma = sigma_est(n, data, theta, psi)
  M = t(c_mat) %*% MASS::ginv(c_mat %*% t(c_mat)) %*% c_mat
  nen = sum(diag(M %*% sigma))
  stat = t(phat) %*% M %*% phat / nen * g(n)
  df_1  = sum(diag(M %*% sigma))^2 / sum(diag(M %*% sigma %*% M %*% sigma))
  df_2 = .f_2(.unsize(data), data, theta, psi)
  df   = c(df_1, df_2)
  crit = qf(1-alpha, df[1], df[2])
  pv   = 1 - pf(stat, df[1], df[2])
  dec = stat > crit
  #dec  = pv < alpha
  return(list(Statistic = stat, df = df, p.value = pv, reject = dec))
}

#' Computes the MCTP for clustered data
#' 
#' @param n A vector containing the sample sizes
#' @param data The data, provided as a list of lists
#' @param p_null A vector or scalar of relative effects under the null hypothesis, defaults to 0.5
#' @param cont A contrast matrix for the relative effects
#' @param normal Logical, whether normal approximation should be used or not, defaults to FALSE
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom and the test decision
max_T_old  = function(n, data, p_null = 0.5, cont, normal = FALSE, theta = NULL, psi = NULL, alpha, type = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  Sigma = sigma_est(n, data, theta = theta, psi = psi)
  p = rel_eff(data, theta, psi)
  R = cov2cor(Sigma)
  R_c = cov2cor(cont%*%Sigma%*%t(cont))
  stat = sqrt(g(n)) * (p - p_null) / sqrt(diag(Sigma))
  if(normal == TRUE) crit = mvtnorm::qmvnorm(1-alpha/2, tail = "lower.tail", mean = rep(0, length(p)), corr = R)$quantile
  if(normal == FALSE) crit = mvtnorm::qmvt(1-alpha/2, tail = "lower.tail", df = g(n) - 1, corr = R)$quantile
  dec = max(abs(stat)) > crit
  return(list(Statistic = stat, df = g(n) - 1, reject = dec))
}


#' Computes the MCTP for clustered data
#' 
#' @param n A vector containing the sample sizes
#' @param data The data, provided as a list of lists
#' @param p_null A vector or scalar of relative effects under the null hypothesis, defaults to 0.5
#' @param cont A contrast matrix for the relative effects
#' @param normal Logical, whether normal approximation should be used or not, defaults to FALSE
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom and the test decision
max_T  = function(n, data, p_null = 0.5, cont, normal = FALSE, theta = NULL, psi = NULL, alpha, type = NULL){
  # if(is.null(theta)){
  #   if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
  #   theta = weight_fun(data, type)$theta
  # } #theta = rep(1/length(data), length(data))
  # # If psi is not provided create it
  # if(is.null(psi)) {
  #   if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
  #   psi = weight_fun(data, type)$psi
  # }
  # Sigma = sigma_est(n, data, theta = theta, psi = psi)
  # p = rel_eff(data, theta, psi)
  # R = cov2cor(Sigma)
  # R_c = cov2cor(cont%*%Sigma%*%t(cont))
  # stat = numeric(nrow(cont))
  # for(i in 1:nrow(cont)){
  #   stat[i] = sqrt(g(n)) * t(cont[i,]) %*% (p - p_null) * (cont %*% Sigma %*% t(cont))[i,i]^(-0.5)
  #   
  # }
  # #stat =  sqrt(g(n)) * (p - p_null) / sqrt(diag(Sigma))
  # 
  # #if(normal == TRUE) crit = mvtnorm::qmvnorm(1-alpha, tail = "lower.tail", mean = rep(0, length(p)), corr = R)$quantile
  # #if(normal == FALSE) crit = mvtnorm::qmvt(1-alpha, tail = "lower.tail", df = g(n) - 1, corr = R)$quantile
  # if(normal == FALSE) rej = 1-mvtnorm::pmvt(upper = rep(max(abs(stat)), nrow(R_c)), df = g(n) - length(n), corr = R_c, keepAttr = F)
  # if(normal == TRUE) rej = 1 - mvtnorm::pmvnorm(upper = rep(max(abs(stat)), nrow(R)), corr = R)
  # #dec = max(abs(stat)) > crit
  # #return(list(Statistic = stat, df = g(n) - 1, reject = dec))
  # dec = rej <= alpha/2
  # return(list(Statistic = stat, df = g(n) - length(n), reject = dec))
  
  
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  Sigma = sigma_est(n, data, theta = theta, psi = psi)
  p = rel_eff(data, theta, psi)
  R = cov2cor(Sigma)
  R_c = cov2cor(cont%*%Sigma%*%t(cont))
  
  stat = p.adj = numeric(nrow(cont))
  for(i in 1:nrow(cont)){
    stat[i] = sqrt(g(n)) * t(cont[i,]) %*% (p - p_null) * (cont %*% Sigma %*% t(cont))[i,i]^(-0.5)
    if(normal == FALSE) p.adj[i] = 1 - mvtnorm::pmvt(lower = -abs(stat[i]), upper = abs(stat[i]), corr = R_c, df = g(n) - length(n), delta = rep(0, nrow(cont)))
    if(normal == T) p.adj[i] = 1 - mvtnorm::pmvnorm(lower = -abs(stat[i]), upper = abs(stat[i]), corr = R_c, mean = rep(0, nrow(cont)))
  }
  
  # for (pp in 1:nc) {
  #   p.adj[pp] <- 1 - pmvt(lower = -abs(T[pp]), abs(T[pp]), 
  #                         corr = rho.bf, df = df.sw, delta = rep(0, nc))
  # }
  
  
  
  #stat =  sqrt(g(n)) * (p - p_null) / sqrt(diag(Sigma))
  
  #if(normal == TRUE) crit = mvtnorm::qmvnorm(1-alpha, tail = "lower.tail", mean = rep(0, length(p)), corr = R)$quantile
  #if(normal == FALSE) crit = mvtnorm::qmvt(1-alpha, tail = "lower.tail", df = g(n) - 1, corr = R)$quantile
  
  #if(normal == FALSE) rej = 1-mvtnorm::pmvt(upper = rep(max(abs(stat)), nrow(R_c)), df = g(n) - length(n), corr = R_c, keepAttr = F)
  #if(normal == TRUE) rej = 1 - mvtnorm::pmvnorm(upper = rep(max(abs(stat)), nrow(R)), corr = R)
  
  #dec = max(abs(stat)) > crit
  #return(list(Statistic = stat, df = g(n) - 1, reject = dec))
  dec = min(p.adj) <= alpha 
  #dec = rej <= alpha/2
  return(list(Statistic = stat, df = g(n) - length(n), reject = dec))
}




q_comb = function(n, data, p_null = 0.5, cont, normal = FALSE, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  
  
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  Sigma = sigma_est(n, data, theta = theta, psi = psi)
  p = rel_eff(data, theta, psi)
  R = cov2cor(Sigma)
  R_c = cov2cor(cont%*%Sigma%*%t(cont))
  
  g_n = g(n)
  

# Wald --------------------------------------------------------------------

  stat_wald = (t(p) %*% t(cont) %*% MASS::ginv(cont %*% Sigma %*% t(cont)) %*% cont %*% p)  * g_n
  df_wald = Matrix::rankMatrix(cont %*% Sigma)
  pv_wald   = 1 - pchisq(stat_wald, df_wald)
  dec_wald  = pv_wald < alpha
  
  

# ANOVA -------------------------------------------------------------------

  M = t(cont) %*% MASS::ginv(cont %*% t(cont)) %*% cont
  nen = sum(diag(M %*% Sigma))
  stat_anv = t(p) %*% M %*% p / nen * g_n
  f_1  = sum(diag(M %*% Sigma))^2 / sum(diag(M %*% Sigma %*% M %*% Sigma))
  f_2  = .f_2(.unsize(data), data, theta, psi)
  df_anv   = c(f_1, f_2)
  pv_anv   = 1 - pf(stat_anv, df_anv[1], df_anv[2])
  dec_anv = pv_anv < alpha
  

# Max-T -------------------------------------------------------------------

  stat_t = p.adj = numeric(nrow(cont))
  for(i in 1:nrow(cont)){
    stat_t[i] = sqrt(g_n) * t(cont[i,]) %*% (p - p_null) * (cont %*% Sigma %*% t(cont))[i,i]^(-0.5)
    if(normal == FALSE) p.adj[i] = 1 - mvtnorm::pmvt(lower = -abs(stat_t[i]), upper = abs(stat_t[i]), corr = R_c, df = g_n - length(n), delta = rep(0, nrow(cont)))
    if(normal == T) p.adj[i] = 1 - mvtnorm::pmvnorm(lower = -abs(stat_t[i]), upper = abs(stat_t[i]), corr = R_c, mean = rep(0, nrow(cont)))
  }
  
  dec_t = min(p.adj) < alpha 
  

#   -----------------------------------------------------------------------

return(c(dec_wald, dec_anv, dec_t))
  
}
