% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/bumphunter.R
\name{bumphuntingEngine}
\alias{bumphuntingEngine}
\title{bumphuntingEngine}
\usage{
bumphuntingEngine(dat.m = NULL, design, sv.m = NULL, chr, pos,
  cluster = NULL, coef = 2, names, cutoff = NULL, qvalue = 0.1,
  maxGap = 500, minDist = 500, robust = FALSE, smooth = FALSE,
  smoothMethod = c("weightedLowess", "loess", "locfit"),
  nullMethod = c("permutation", "bootstrap"), B = 10000, corr = FALSE,
  corFunc = c("spearman", "pearson", "kendall"), combp = FALSE,
  merge = c("single", "complete", "average"), corr.cutoff = 0.8,
  combine = c("stouffer_liptak", "zscore"), verbose = TRUE, ...)
}
\arguments{
\item{dat.m}{450k methylation microarray matrix via preprocessing; m x n matrix}

\item{design}{design matrix; [1, covariate of interested, paired info]}

\item{sv.m}{surrogate variables matrix}

\item{chr}{chromosome vector for probes}

\item{pos}{position vector for probes}

\item{cluster}{The clusters of locations that are to be analyzed together. via correlation
perspevtive, we can use agglomerate cluster index to bump analysis. If cluster is not
available, clusterMaker can be used to cluster nearby locations}

\item{coef}{Integer. covariate of interest's column}

\item{names}{probe names vector ; used in function \link{clusterMaker} and \link{regionSeeker}}

\item{cutoff}{numeric value. Value of estimates coefficient of covariate of interested above cutoff
or below the negative of cutoff will be used as candidate bump regions.}

\item{qvalue}{numeric value. cutoff of pvalue for candidate regions selection ||
qvalue cutoff for combination pvalue method}

\item{maxGap}{if cluster is not availabel. maxGap is used by clusterMaker to define cluster}

\item{minDist}{if clusters are build, the mininal distance between clusters is setted as constraint}

\item{robust}{logic, use robust linear regression or not; Default FALSE}

\item{smooth}{logic. If TRUE then the standard error or correlation of point-wise estimatrs will be used as weigths
in the \link{smoother}}

\item{smoothMethod}{local regression method used in \link{smoother}}

\item{nullMethod}{Method for generating null candidate regions. If ncol(design) > 2. bootstrap method is recommanded}

\item{B}{integer; Denoting the number of resamples to computr null distribution. Default = 0
B is used to generate permutation matrix,
which describes permutations to generate null distribution}

\item{corr}{logical; If TRUE then position correlation matrix is considered as weighted matrix
and clusters modified by correlation constraint; Default FALSE}

\item{corFunc}{Optional; "spearman"(Default) and "pearson"}

\item{combp}{logical; If TRUE then slk p correction will be applied}

\item{merge}{how to merge two sub-clusters; c("single", "complete", "average")}

\item{corr.cutoff}{correlation cutoff in merge algorithm}

\item{combine}{combine pvalue method; c("stouffer_liptak", "zscore")}

\item{verbose}{logical. Optional printing progress message or not}

\item{...}{}
}
\value{
bumps object
}
\description{
bump hunting algorithm for 450k methylation microarray
}
\details{
bumphuntingEngine
}
\examples{
# cluster input format: cg12045430 cg20826792 cg00381604 cg20253340
}

