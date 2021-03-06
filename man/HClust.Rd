% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/misc.R
\name{HClust}
\alias{HClust}
\title{HClust}
\usage{
HClust(data = NULL, method = "average", distance = "euclidean", p = 2,
  sign = c("S", "U"))
}
\arguments{
\item{data}{matrix of beta data}

\item{method}{The agglomeration method to be used. This must be one of "ward", "single", "complete" or "average"}

\item{distance}{The distance measure to be used. This must be one of "euclidiean", "manhattan", "maximum", or "minkowski".}

\item{p}{power of the Minkowski distance.}

\item{sign}{Optional 'U' for unsigned and 'S' for signed dissimilarity}
}
\value{
An object of class *hclust* which describes the tree produced by the clustering process
}
\description{
Hierachical Clustering function based on fast parallel hierachical clustering
}
\author{
Xin Zhou \url{xxz220@miami.edu}
}

