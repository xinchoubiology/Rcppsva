##' @title print
##' @rdname print-methods
##' @exportMethod print
setMethod("print", signature(x = "bumps"),
          function(x, ...){
            cat(sprintf("'bumps' object with %d bumps\n", ncol(x@table)))
            cat(sprintf("algorithms : ...... \n"))
            for(params in names(x@algorithm)){
              cat(sprintf("\t%s = %s \n", params, x@algorithm[[params]]))
            }
          })


##' bumphuntingEigen
##' 
##' @title bumphuntingEigen
##' @description bump hunting algorithm for 450k methylation microarray
##' @param dat.m 450k methylation microarray matrix via preprocessing; m x n matrix
##' @param design design matrix; [covariate of interested, paired info, surrogate variables]
##' @param chr chromosome vector for probes
##' @param pos position vector for probes
##' @param cluster The clusters of locations that are to be analyzed together. via correlation
##'        perspevtive, we can use agglomerate cluster index to bump analysis. If cluster is not
##'        available, clusterMaker can be used to cluster nearby locations
##' @param coef Integer. covariate of interest's column
##' @param cutoff numeric value. Value of estimates coefficient of covariate of interested above cutoff
##'        or below the negative of cutoff will be used as candidate bump regions. 
##' @param pvalue numeric value. cutoff of pvalue for candidate regions selection
##' @param maxGap if cluster is not availabel. maxGap is used by clusterMaker to define cluster
##' @param minDist if clusters are build, the mininal distance between clusters is setted as constraint
##' @param nullMethod Method for generating null candidate regions. If ncol(design) > 2. bootstrap method is recommanded
##' @param smooth logic. If TRUE then the standard error or correlation of point-wise estimatrs will be used as weigths 
##'        in the \code{loessByCluster}
##' @param smoother smooth function to estimate genome profile.
##' @param B integer; Denoting the number of resamples to computr null distribution. Default = 0
##' @param permutation matrix; define matrix to describe permutations to generate null distribution
##' @param verbose logical. Optional printing progress message or not
##' @export
bumphuntingEigen <- function(dat.m = NULL, design, chr, pos, cluster = NULL, coef = 2,
                             cutoff = NULL, pvalue = 0.01, maxGap = 500, minDist = 500,
                             nullMethod = c("permutation", "bootstrap"),
                             smooth = FALSE, smoother = bumphunter::loessByCluster,
                             B = 10000, permutation = NULL, verbose = TRUE){
  
}