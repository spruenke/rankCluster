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
  st   = .q_wald_arma(n, data, theta, psi, cont)
  stat = g(n) * st[[1]]
  #df   = Matrix::rankMatrix(cont%*%Sigma)
  df   = st[[2]]
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
q_anova = function(n, data, cont, f_2, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  st   = .q_anova_arma(n, data, theta, psi, cont)
  stat =  st[[1]]  * g(n)
  #df   = c(sum(diag(M%*%Sigma))^2 / sum(diag(M%*%Sigma%*%M%*%Sigma)), f_2)
  df   = c(st[[2]], f_2)
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
  #if(normal == TRUE) crit = mvtnorm::qmvnorm(1-alpha, tail = "lower.tail", mean = rep(0, length(p)), corr = R)$quantile
  #if(normal == FALSE) crit = mvtnorm::qmvt(1-alpha, tail = "lower.tail", df = g(n) - 1, corr = R)$quantile
  if(normal == FALSE) rej = 1-mvtnorm::pmvt(upper = rep(max(abs(stat)), nrow(R)), df = g(n) - length(n), corr = R, keepAttr = F)
  #dec = max(abs(stat)) > crit
  #return(list(Statistic = stat, df = g(n) - 1, reject = dec))
  dec = rej <= alpha/2
  return(list(Statistic = stat, df = g(n) - length(n), reject = dec))
}






