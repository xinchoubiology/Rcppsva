% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/optimization-methods.R
\name{gapStat}
\alias{gapStat}
\title{gapStat}
\usage{
gapStat(data = NULL, dendro = NULL, dendref = NULL, cmax = 0.24,
  cmin = 0.04, mcore = 2, log = FALSE, verbose = TRUE, ...)
}
\arguments{
\item{data}{methylation data matrix}

\item{dendro}{actual dendrogram}

\item{dendref}{reference dendrogram}

\item{cmax}{maximum dissimilarity Default 0.22}

\item{cmin}{minimun dissimilarity Default 0.02}

\item{mcore}{number of threads running}

\item{log}{FALSE(Default). Standardize by log or not}

\item{verbose}{verbose}

\item{...}{optional parameters for \code{distributeRef} and \code{HClust}}
}
\value{
list of clusterDetect
        labels
        number
        Wk  within sum of dissimilarity
        EWk expectation of Wk under null reference distribution
}
\description{
Based on dissmilarity measurement, we define Wk = sum(sum(d_{ij})/2n_r)
}
\details{
grid search for CLuster number detection by gap statistic
}
\author{
Xin Zhou \url{xxz220@miami.edu}
}

