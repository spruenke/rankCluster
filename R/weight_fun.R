#' Empirical Distribution function with cluster weights psi
#'
#' @param data The data, provided as a list of lists
#' @param type Character, specifying whether weighted or unweighted estimator will be used. Defaults to unweighted. The first letter is enough.
#' @return A list of vectors with the weights.
weight_fun = function(data, type = "unweighted"){
  psi = list()
  theta = numeric(length(data))
  if(type == "unweighted" || type == "u"){
    for(i in 1:length(data)){
      psi[[i]] = rep(1 / length(data[[i]]), length(data[[i]]))
    }
    theta = rep(1/length(data), length(data))
  }
  if(type == "weighted" || type == "w"){
    for(i in 1:length(data)){
      psi[[i]] = sapply(data[[i]], FUN = function(x) {
        length(x) / sum(sapply(data[[i]], length))
      }) 
    }
    theta = sapply(data, FUN = function(x){
      length(x) / sum(sapply(data, length))
    })
  }
  return(list("theta" = theta, "psi" = psi))
}