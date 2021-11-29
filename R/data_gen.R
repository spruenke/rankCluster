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
               copu        = copula::mvdc(copula = copula::normalCopula(copula::P2p(R), dim = m[[i]][j], dispstr = "un"), margins = rep(dist, m[[i]][j]), paramMargins = rep(list(params), m[[i]][j]))
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
               copu        = copula::mvdc(copula = copula::normalCopula(P2p(R), dim = m[[i]][j], dispstr = "un"), margins = rep(dist, m[[i]][j]), paramMargins = rep(list(params), m[[i]][j]))
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
    n = sample(c(8, 25, rep(n_i, (nn-2))))
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

#' Function generating sample- and cluster sizes with only mild unbalancing
#' 
#' @param nn Number of groups/samples
#' @param n_i Typical sample sizes
#' @param m_ij Typical cluster size
#' @param each_s Logical: Does each sample contain a large/small cluster? Default: FALSE
#' @param both_s Logical: Does a sample contain both large/small clusters or only one of them? Default: TRUE
#' @param identical_s Logical: Identical Sample sizes? Default: TRUE
#' @param identical_c Logical: Identical Cluster sizes? Default: TRUE
#' @return List containing the sample- and cluster sizes.
nm_gen2 = function(nn, n_i, m_ij, each_s = F, both_s = T, identical_s = T, identical_c = T){
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
      #if(each_s == F) other = sample(c(1:nn), 2)
      #if(each_s == T) other = list(c(1:nn))
      for(i in 1:nn){
        #if(both_s == T){
          #if(i %in% other[[1]]){
            m[[i]] = sample(c(m_ij:(m_ij - 3), m_ij:(m_ij + 3)), size = n_i, replace = T)#c(m_small, m_large, rep(m_ij, (n_i-2))))
          #} else {
            #m[[i]] = rep(m_ij, n_i)
          #}
       # } else if(both_s == F){
         # if(i == other[1]){
            #m[[i]] = sample(c(m_small, rep(m_ij, (n_i-1))))
         # } else if(i == other[2]){
            #m[[i]] = sample(c(m_large, rep(m_ij, (n_i-1))))
         # } else {
            #m[[i]] = rep(m_ij, n_i)
         # }
        #}
      }
    }
  }
  
  if(identical_s == F){
    n = sample(c(n_i:(n_i - 3), n_i:(n_i + 3)), size = nn, replace = T)
    m = list()
    if(identical_c == T){
      for(i in 1:nn){
        m[[i]] = rep(m_ij, n_i)
      }
      
    } else if(identical_c == F){
      m_small = 3
      m_large = 15
      #if(each_s == F) other = sample(c(1:nn), 2)
      #if(each_s == T) other = list(c(1:nn))
      for(i in 1:nn){
        #if(both_s == T){
        #if(i %in% other[[1]]){
        m[[i]] = sample(c(m_ij:(m_ij - 3), m_ij:(m_ij + 3)), size = n_i, replace = T)#c(m_small, m_large, rep(m_ij, (n_i-2))))
        #} else {
        #m[[i]] = rep(m_ij, n_i)
        #}
        # } else if(both_s == F){
        # if(i == other[1]){
        #m[[i]] = sample(c(m_small, rep(m_ij, (n_i-1))))
        # } else if(i == other[2]){
        #m[[i]] = sample(c(m_large, rep(m_ij, (n_i-1))))
        # } else {
        #m[[i]] = rep(m_ij, n_i)
        # }
        #}
      }
    }
  }
  return(list(n, m))
}