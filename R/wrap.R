#' Empirical Distribution function with cluster weights psi
#'
#' @param x A scalar or vector to be checked as F(x)
#' @param i The group/sample index
#' @param data The data, provided as a list of lists
#' @param psi A vector of weights, defaults to unweighted estimator
#' @return Value(s) of the empirical distribution function 

f_psi = function(x, i, data, psi = NULL){
  #m = numeric(length(data[[i]])) #empty vector for means of clusters
  # If psi is not provided create it
  if(is.null(psi)) psi = rep(1/length(data[[i]]), length(data[[i]]))
  m = .f_psi_arma(x, data[[i]], psi)
  return(m)
}


#' Empirical Reference Distribution Function with group weights theta and cluster weights psi
#' 
#' @param x A scalar or vector to be checked as F(x)
#' @param data The data, provided as a list of lists
#' @param theta A vector of group weights, defaults to unweighted estimator
#' @param psi A list of vector-weights, defaults to unweighted estimator
#' @return Value(s) of the empirical distribution function
f_theta = function(x, data, theta = NULL, psi = NULL){
  if(is.null(theta)) theta = rep(1/length(data), length(data))
  if(is.null(psi)){
    psi = list()
    for(i in 1:length(data)){
      psi[[i]] = rep(1 / length(data[[i]]), length(data[[i]]))
    }
  }
  res = .f_theta_arma(x, data, theta, psi)
  return(c(res))
}


#' Computes relative effects with respect to the reference distribution
#' 
#' @param data The data, provided as a list of lists
#' @param theta A vector of group weights, defaults to unweighted estimator
#' @param psi A list of vector-weights, defaults to unweighted estimator
#' @return A vector of the relative effects
rel_eff = function(data, theta = NULL, psi = NULL){
  if(is.null(psi)){
    psi = list()
    for(i in 1:length(data)){
      psi[[i]] = rep(1 / length(data[[i]]), length(data[[i]]))
    }
  }
  if(is.null(theta)) theta = rep(1/length(data), length(data))
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
#' @return A Variance-Covariance-Matrix
sigma_est = function(n, data, theta = NULL, psi = NULL){
  if(is.null(psi)){
    psi = list()
    for(i in 1:length(data)){
      psi[[i]] = rep(1 / length(data[[i]]), length(data[[i]]))
    }
  }
  if(is.null(theta)) theta = rep(1/length(data), length(data))
  return( .sigma_est_arma(n, data, theta, psi))
}


#' Computes the Wald-Test for clustered data
#' 
#' @param n A vector containing the sample sizes
#' @param data The data, provided as a list of lists
#' @param cont A contrast matrix for the relative effects
#' @param theta A Vector for the group weights, defaults to unweighted estimator
#' @param psi A list of vectors with the cluster weights, defaults to unweighted estimator
#' @param alpha The significance level, defaults to 0.05
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-value and the test decision
q_wald = function(n, data, cont, theta = NULL, psi = NULL, alpha = 0.05){
  if(is.null(psi)){
    psi = list()
    for(i in 1:length(data)){
      psi[[i]] = rep(1 / length(data[[i]]), length(data[[i]]))
    }
  }
  if(is.null(theta)) theta = rep(1/length(data), length(data))
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
#' @return A list containing the value of the test statistic, the degrees of freedom, the p-value and the test decision
q_anova = function(n, data, cont, f_2, theta = NULL, psi = NULL, alpha = 0.05){
  if(is.null(psi)){
    psi = list()
    for(i in 1:length(data)){
      psi[[i]] = rep(1 / length(data[[i]]), length(data[[i]]))
    }
  }
  if(is.null(theta)) theta = rep(1/length(data), length(data))
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
#' @return A list containing the value of the test statistic, the degrees of freedom and the test decision
max_T  = function(n, data, p_null = 0.5, cont, normal = FALSE, theta = NULL, psi = NULL, alpha){
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
#' @return A list containing the value of the test statistic, the degrees of freedom and the test decision
max_T2  = function(n, data, p_null = 0.5, cont, normal = FALSE, theta = NULL, psi = NULL, alpha){
  Sigma = sigma_est(n, data, theta = theta, psi = psi)
  p = rel_eff(data, theta, psi)
  R = cov2cor(Sigma)
  R_c = cov2cor(cont%*%Sigma%*%t(cont))
  stat = sqrt(g(n)) * (p - p_null) / sqrt(diag(Sigma))
  #if(normal == TRUE) crit = mvtnorm::qmvnorm(1-alpha, tail = "lower.tail", mean = rep(0, length(p)), corr = R)$quantile
  #if(normal == FALSE) crit = mvtnorm::qmvt(1-alpha, tail = "lower.tail", df = g(n) - 1, corr = R)$quantile
  if(normal == FALSE) rej = 1-mvtnorm::pmvt(upper = abs(stat), df = round(sqrt(g(n) - 1)), corr = R)
  #dec = max(abs(stat)) > crit
  #return(list(Statistic = stat, df = g(n) - 1, reject = dec))
  dec = rej <= alpha
  return(list(Statistic = stat, df = g(n) - 1, reject = dec))
}

#' Data-generating function under the null hypothesis
#' 
#' @param n Vector of sample sizes
#' @param m Vector of cluster sizes, defaults to sample size
#' @param dist Distribution under the null, defaults to normal
#' @param corstruct Type of correlation (string), defaults to "independent". Alternatives are "exchangeable" and "wild"
#' @param rho Fixed correlation, defaults to NULL
#' @param rho.max In case of no fixed correlation, specifies the maximal correlation, defaults to NULL
#' @param params List of params necessary for the distribution specified in dist, defaults to standard normal
#' @return Returns a list of lists with the data
h_0_f = function(n = rep(5, 5), m = NULL, dist = "norm", corstruct = "independent", rho = NULL, rho.max = NULL, params = list(mean = 0, sd = 1)){
  if(is.null(m)){
    m = lapply(n, FUN = function(x){
      return(rep(n, 1))
    })
  }
  switch(corstruct,
         independent = {
           l = list()
           for(i in 1:length(n)){
             l[[i]] = list()
             for(j in 1:length(m[[i]])){
               l[[i]][[j]] = do.call(paste0("r", dist), c(m[[i]][j], (params)))
             }
           }
         },
         exchangeable = {
           if(is.null(rho)){
             if(is.null(rho.max)){
               rho = runif(1, 0.05, 0.95)
             } else {
               rho = runif(1, 0, rho.max)
             }
           }
           l = list()
           for(i in 1:length(n)){
             l[[i]] = list()
             for(j in 1:length(m[[i]])){
               R           = matrix(rho, nrow = m[[i]][j], ncol = m[[i]][j])
               diag(R)     = rep(1, nrow(R))
               copu        = copula::mvdc(copula = normalCopula(P2p(R), dim = m[[i]][j], dispstr = "un"), margins = rep(dist, m[[i]][j]), paramMargins = rep(list(params), m[[i]][j]))
               l[[i]][[j]] = copula::rMvdc(1, copu)
             }
           }
         },
         wild = {
           
           l = list()
           for(i in 1:length(n)){
             l[[i]] = list()
             for(j in 1:length(m[[i]])){
               if(is.null(rho)){
                 if(is.null(rho.max)){
                   rho = runif(m[[i]][j], 0.05, 0.95)
                 } else {
                   rho = runif(m[[i]][j], 0, rho.max)
                 }
               }
               R = rho %*% t(rho)
               diag(R)     = rep(1, nrow(R))
               copu        = copula::mvdc(copula = normalCopula(P2p(R), dim = m[[i]][j], dispstr = "un"), margins = rep(dist, m[[i]][j]), paramMargins = rep(list(params), m[[i]][j]))
               l[[i]][[j]] = copula::rMvdc(1, copu)
             }
           }
         }
         
         
  )
  
  return(l)
}

#' Function generating sample- and cluster sizes
#' 
#' @param nn Number of groups/samples
#' @param n_i Typical sample sizes
#' @param m_ij Typical cluster size
#' @param each_s Logical: Does each sample contain a large/small cluster? Default: FALSE
#' @param both_s Logical: Does a sample contain both large/small clusters or only one of them? Default: TRUE
#' @param identical_s Logical: Identical Sample sizes? Default: TRUE
#' @param identical_c Logical: Identical Cluster sizes? Default: TRUE
#' @return List containing the sample- and cluster sizes.
nm_gen = function(nn, n_i, m_ij, each_s = F, both_s = T, identical_s = T, identical_c = T){
  # each_s: Does each sample contain large/small cluster?
  # both_s: Does a sample contain both large/small or only one?
  # identical_s: Identical Sample Sizes?
  # identical_c: Identical Clustersizes?
  if(each_s == T && both_s == F) stop("Invalid conditions")
  if(identical_s == T){
    n = rep(n_i, nn)
    m = list()
    if(identical_c == T){
      for(i in 1:nn){
        m[[i]] = rep(m_ij, n_i)
      }
      
    } else if(identical_c == F){
      m_small = 3
      m_large = 15
      if(each_s == F) other = sample(c(1:nn), 2)
      if(each_s == T) other = list(c(1:nn))
      for(i in 1:nn){
        if(both_s == T){
          if(i %in% other[[1]]){
            m[[i]] = sample(c(m_small, m_large, rep(m_ij, (n_i-2))))
          } else {
            m[[i]] = rep(m_ij, n_i)
          }
        } else if(both_s == F){
          if(i == other[1]){
            m[[i]] = sample(c(m_small, rep(m_ij, (n_i-1))))
          } else if(i == other[2]){
            m[[i]] = sample(c(m_large, rep(m_ij, (n_i-1))))
          } else {
            m[[i]] = rep(m_ij, n_i)
          }
        }
      }
    }
  }
  
  if(identical_s == F){
    n = sample(c(3, 25, rep(n_i, (nn-2))))
    m = list()
    if(identical_c == T){
      for(i in 1:nn){
        m[[i]] = rep(m_ij, n[i])
      }
      
    } else if(identical_c == F){
      m_small = 3
      m_large = 15
      if(each_s == F) other = sample(c(1:nn), 2)
      if(each_s == T) other = list(c(1:nn))
      for(i in 1:nn){
        if(both_s == T){
          if(i %in% other[[1]]){
            m[[i]] = sample(c(m_small, m_large, rep(m_ij, (n[i]-2))))
          } else {
            m[[i]] = rep(m_ij, n[i])
          }
        } else if(both_s == F){
          if(i == other[1]){
            m[[i]] = sample(c(m_small, rep(m_ij, (n[i]-1))))
          } else if(i == other[2]){
            m[[i]] = sample(c(m_large, rep(m_ij, (n[i]-1))))
          } else {
            m[[i]] = rep(m_ij, n[i])
          }
        }
      }
    }
  }
  return(list(n, m))
}




