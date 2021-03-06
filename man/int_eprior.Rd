% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{int_eprior}
\alias{int_eprior}
\title{int_eprior}
\usage{
int_eprior(sx, ghat, dhat)
}
\arguments{
\item{sx}{standard data matrix}

\item{ghat}{estimated gamma values}

\item{dhat}{estimated delta values}
}
\value{
list
        g.star
        d.star
}
\description{
Integrating empirical bayesian prior function
by Monte Carlo integration
}
\author{
Xin Zhou \url{xinchoubiology@gmail.com}
}

