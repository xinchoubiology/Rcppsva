% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/dendrogram.R
\name{dendrogram}
\alias{dendrogram}
\title{dendrogram}
\usage{
dendrogram(x, type = c("rectangle", "triangle"), ...)
}
\arguments{
\item{x}{object of hclust. Derived from \code{\link{stats}{dendrogram}}}

\item{type}{type of plot, c("rectangle", "triangle"); Default rectangle}

\item{...}{ignore}
}
\value{
list
        segments data.frame
        labels data.frame
        class
}
\description{
extracting data frame from hclust object for plotting dendrgram
}

