% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weight_fun.R
\name{weight_fun}
\alias{weight_fun}
\title{Empirical Distribution function with cluster weights psi}
\usage{
weight_fun(data, type = "unweighted")
}
\arguments{
\item{data}{The data, provided as a list of lists}

\item{type}{Character, specifying whether weighted or unweighted estimator will be used. Defaults to unweighted. The first letter is enough.}
}
\value{
A list of vectors with the weights.
}
\description{
Empirical Distribution function with cluster weights psi
}
