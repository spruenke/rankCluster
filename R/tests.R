#' Computes the Wald-Test for clustered data
#' 
#' @param data The data, provided as a list of lists
#' @param cont A contrast matrix for the relative effects
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-value and the test decision
waldTest = function(data, cont, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  n = .unsize(data)[[1]]
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  c_mat = cont
  p_hat = rel_eff(data, theta, psi)
  sigma = sigma_est(data, theta, psi)
  stat = (t(p_hat) %*% t(c_mat) %*% MASS::ginv(c_mat %*% sigma %*% t(c_mat)) %*% c_mat %*% p_hat)  * g(n)
  df = Matrix::rankMatrix(c_mat %*% sigma)
  pv   = 1 - pchisq(stat, df)
  dec  = pv < alpha
  return(list("Statistic" = stat, "df" = df, "p.value" = pv, "reject" = dec))
}


#' Computes the ANOVA-Test for clustered data
#' 
#' @param data The data, provided as a list of lists
#' @param cont A contrast matrix for the relative effects
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-value and the test decision
anovaTest = function(data, cont, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  n = .unsize(data)[[1]]
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } 
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  c_mat = cont
  p_hat = rel_eff(data, theta, psi)
  sigma = sigma_est(data, theta, psi)
  M = t(c_mat) %*% MASS::ginv(c_mat %*% t(c_mat)) %*% c_mat
  nen = sum(diag(M %*% sigma))
  stat = t(p_hat) %*% M %*% p_hat / nen * g(n)
  df_1  = sum(diag(M %*% sigma))^2 / sum(diag(M %*% sigma %*% M %*% sigma))
  df_2 = .f_2(data, theta, psi)
  df   = c(df_1, df_2)
  crit = qf(1-alpha, df[1], df[2])
  pv   = 1 - pf(stat, df[1], df[2])
  dec = stat > crit
  #dec  = pv < alpha
  return(list("Statistic" = stat, "df" = df, "p.value" = pv, "reject" = dec))
}


#' Computes the MCTP for clustered data
#' 
#' @param data The data, provided as a list of lists
#' @param p_null A vector or scalar of relative effects under the null hypothesis, defaults to 0.5
#' @param cont A contrast matrix for the relative effects
#' @param normal Logical, whether normal approximation should be used or not, defaults to FALSE
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-values and the test decision
mctp  = function(data, p_null = 0.5, cont, normal = FALSE, theta = NULL, psi = NULL, alpha, type = NULL){
  n = .unsize(data)[[1]]
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  Sigma = sigma_est(data, theta = theta, psi = psi)
  p_hat = rel_eff(data, theta, psi)
  R = cov2cor(Sigma)
  R_c = cov2cor(cont%*%Sigma%*%t(cont))
  
  stat_t = sqrt(g(n)) * (cont) %*% (p - p_null) * diag((cont %*% Sigma %*% t(cont))^(-0.5))
  df_t = floor(.df_sw(data, theta, psi, cont))
  
  p_vals_t = numeric(nrow(cont))
  for(i in 1:nrow(cont)){
    p_vals_t[i] = 1 - mvtnorm::pmvt(lower = -abs(stat_t[i]), upper = abs(stat_t[i]), corr = R_c, df = df_t, delta = rep(0, nrow(cont)))
  }
  
  p_vals_mctp = data.frame("p.value" = c(p_vals_t, min(p_vals_t)), row.names = c(row.names(cont), "Overall"))
  pv_t = min(p_vals_t)
  dec_t = pv_t < alpha

  return(list("Statistic" = stat_t, "df" = df_t, "p.value" = p_vals_mctp, "reject" = dec_t))

}


#' Computes comprehensive analysis of the cluster data, containing descriptive statistics, confidence intervals and different test procedures.
#' 
#' @param data The data, provided as a list of lists
#' @param p_null A vector or scalar of relative effects under the null hypothesis, defaults to 0.5
#' @param contrast A contrast statement, either given as a numeric matrix or a character string (e.g. "Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams", "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean"). Defaults to the "GrandMean".
#' @param normal Logical, whether normal approximation should be used or not, defaults to FALSE
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @param type A string indicating whether weighted or unweighted estimator should be used. Only if psi is not provided
clusterTest = function(data, p_null = 0.5, contrast = NULL, normal = FALSE, theta = NULL, psi = NULL, alpha = 0.05, type = NULL){
  n = .unsize(data)[[1]]
  if(is.null(contrast)) cont = multcomp::contrMat(n, "GrandMean")
  if(is.matrix(contrast)){
    if(all(rowSums(contrast) == 0)){
      cont = contrast
    } else {
      stop("Contrasts do not sum up to 0 (rowwise)")
    }
  }
  if(is.character(contrast)) cont = multcomp::contrMat(n, contrast)
  
  if(is.null(theta)){
    if(is.null(type)) theta = weight_fun(data, "unweighted")$theta
    theta = weight_fun(data, type)$theta
  } #theta = rep(1/length(data), length(data))
  # If psi is not provided create it
  if(is.null(psi)) {
    if(is.null(type)) psi = weight_fun(data, "unweighted")$psi
    psi = weight_fun(data, type)$psi
  }
  Sigma = sigma_est(data, theta = theta, psi = psi)
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
  f_2  = .f_2(data, theta, psi)
  df_anv   = c(f_1, f_2)
  pv_anv   = 1 - pf(stat_anv, df_anv[1], df_anv[2])
  dec_anv = pv_anv < alpha
  

# Max-T -------------------------------------------------------------------

  stat_t = sqrt(g(n)) * (cont) %*% (p - p_null) * diag((cont %*% Sigma %*% t(cont))^(-0.5))
  df_t = floor(.df_sw(data, theta, psi, cont))
  
  p_vals_t = numeric(nrow(cont))
  for(i in 1:nrow(cont)){
    p_vals_t[i] = 1 - mvtnorm::pmvt(lower = -abs(stat_t[i]), upper = abs(stat_t[i]), corr = R_c, df = df_t, delta = rep(0, nrow(cont)))
  }
  
  p_vals_mctp = data.frame("p.value" = c(p_vals_t, min(p_vals_t)), row.names = c(row.names(cont), "Overall"))
  pv_t = min(p_vals_t)
  dec_t = pv_t < alpha

# Results  -----------------------------------------------------------------------

  wald_list  = list("Reject" = dec_wald, "Statistic" = stat_wald, "df" = df_wald, "p.value" = pv_wald)
  anova_list = list("Reject" = dec_anv, "Statistic" = stat_anv, "df" = df_anv,  "p.value" = pv_anv)
  mctp_list  = list("Reject" = dec_t, "Statistic" = stat_t, "df" = df_t, "p.value" = pv_t, "p.values" = p_vals_mctp)
  
  p_vals = data.frame("p-value" = c(pv_wald, pv_anv, pv_t), row.names = c("Wald", "ANOVA", "MCTP"))
  
  conf_crit  = mvtnorm::qmvt(alpha, df = df_t, delta = rep(0, length(p)), sigma = Sigma)
  conf_int_l = p - conf_crit / sqrt(g_n) * sqrt(diag(Sigma))
  conf_int_u = p + conf_crit / sqrt(g_n) * sqrt(diag(Sigma))
  
  desc   = data.frame("Estimator" = p, "Std.Dev" = sqrt(diag(Sigma)), "Lower" = conf_int_l, "Upper" = conf_int_u)
  
  return_list = list("Descriptive" = desc, "p.value" = p_vals, wald_list, anova_list, mctp_list)
}
