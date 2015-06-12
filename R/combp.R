##' calculation combination pvalue of given regions(pvalues within region are correlated)
##' 
##' @title combine.pvalue
##' @description combine-pvalue accepts a list of pvalues of probes, and correlation matrices
##'              are claculated by input data matrix. Then, it estimate combined pvalue for
##'              candidate regions defined in input region list
##' @param dat.m n x m delta M|beta matrix for n CpG sites across 2*m paired different patient samples
##' @param pvalues pvalue column vector for each probe
##' @param cluster correlated cluster list calculated by \link{corrclusterMaker}
##'        cluster list contains: {single probe cluster, multiple probes cluster}
##' @param chr chromosome vector
##' @param pos position vector
##' @param names probe names vector ; used in function \link{clusterMaker}
##' @param method correlation calculation method; c("spearman", "pearson", "kendall"); Default "spearman"
##' @param combine combine pvalue method; c("stouffer_liptak", "zscore")
##' @param weight NULL; weight provided as \code{1/sigma} of each probe. It means each probe's contribution 
##'        to defined regions
##' @param cutoff cutoff for dmr(differential methylated regions) detection's qvalue(fdr); 
##'        Default 0.1
##' @param cores core number for R backend combp algorithm
##' @return combine-pval
##' @importFrom plyr llply
##' @importFrom qvalue qvalue
##' @importFrom parallel detectCores
##' @importFrom doMC registerDoMC
##' @import data.table
##' @export
combine.pvalue <- function(dat.m = NULL, pvalues, cluster, chr, pos, names,
                           method = c("spearman", "pearson", "kendall"), 
                           combine = c("stouffer_liptak", "zscore"),
                           weight = NULL, cutoff = 0.1,
                           cores = detectCores()){
  combine  <- match.arg(combine)
  method   <- match.arg(method)
  ## regions : cluster with >= 2 probes
  rIndexes <- which(sapply(cluster, length) >= 2)
  multiregions <- cluster[rIndexes]
  
  registerDoMC(cores = cores)
  if(combine == "stouffer_liptak"){
    Cp <- llply(multiregions, 
                .fun = function(ix){
                          sigma <- cov(t(dat.m[ix,]), method = method)
                          combp <- stouffer_liptak.combp(pvalues[ix,], sigma, weight)
                          combp
                       }, 
                .parallel = TRUE)
  }else{
    Cp <- llply(multiregions,
                .fun = function(ix){
                          sigma <- cov(t(dat.m[ix,]), method = method)
                          combp <- zscore.combp(pvalues[ix,], sigma, weight)
                          combp
                       }, 
                .parallel = TRUE)
  }
  ## tabulate 
  Cp <- c(unlist(Cp), pvalues[unlist(cluster[-rIndexes]),])
  Cqval <- qvalue(Cp)$qvalue
  
  cluster <- c(multiregions, cluster[-rIndexes])
  pT <- data.table(names = names, chr = chr, pos = pos)
  setkey(pT, names)
  ## differential
  segments <- list("diff" = which(Cqval <= cutoff), "null" = which(Cqval > cutoff))
  res <- vector("list", 2)
  for(i in 1:2){
    res[[i]] <- data.table(chr    = sapply(segments[[i]], function(ix) pT[cluster[ix]]$chr[1]),
                           start  = sapply(segments[[i]], function(ix) min(pT[cluster[ix]]$pos)),
                           end    = sapply(segments[[i]], function(ix) max(pT[cluster[ix]]$pos) + 1),
                           length = sapply(segments[[i]], function(ix) max(pT[cluster[ix]]$pos) + 1 - min(pT[cluster[ix]]$pos)),
                           L      = sapply(segments[[i]], function(ix) length(cluster[[ix]])),
                           probes = sapply(segments[[i]], function(ix) paste0(cluster[[ix]], collapse = ";")),
                           pvalue = sapply(segments[[i]], function(ix) Cp[ix]),
                           fdr    = sapply(segments[[i]], function(ix) Cqval[ix]))
  }
  names(res) <- names(segments)
  return(res)
}

##' 2 combp functions are borrowed from Combp value function of brentp
##' 
##' Calculate combined p-value by stouffer_liptak method
##' 
##' @title stouffer_liptak.combp
##' @param pvalues vector of pvalues
##' @param sigma covariance of pvalue
##' @param weight use weight for pvalue combination; Default NULL
##' @details if \code{sigma} is not positive definitive, off-diag elements x 0.99
##' @return Cp (Combined pvalue)
##' @export
stouffer_liptak.combp <- function(pvalues, sigma, weight = NULL){
  qvalues <- qnorm(pvalues, mean = 0, sd = 1, lower.tail = TRUE)
  C <- try(chol(sigma), silent = TRUE)
  if(inherits(C, "try-error")){
    sigma <- sigma * 0.9999 + diag(0.0001 * diag(sigma))
    C <- chol(sigma)
  }
  qvalues <- .Internal(backsolve(C, qvalues, ncol(C), TRUE, FALSE))
  if(is.null(weight)){
    Cq <- sum(qvalues) / sqrt(length(qvalues))
  }else{
    Cq <- sum(weight * qvalues) / sqrt(sum(weight^2))
  }
  return(pnorm(Cq, lower.tail = TRUE))
}

##' Calculate the combined pvalue by combined Z score method
##' 
##' @title zscore.combp
##' @param pvalues vector of pvalues
##' @param sigma covariance of pvalue
##' @param weight use weight for pvalue combination; Default NULL
##' @details if \code{sigma} is not positive definitive, off-diag elements x 0.99
##' @return Cp (Combined pvalue)
##' @export
zscore.combp <- function(pvalues, sigma, weight = NULL){
  if(is.null(weight)){
    z <- weighted.mean(qnorm(pvalues, lower.tail = TRUE), weight)
  }else{
    z <- mean(qnorm(pvalues, lower.tail = TRUE))
  }
  sz <- sqrt(length(pvalues) + 2 * sum(sigma[lower.tri(sigma)])) / length(pvalues)
  return(pnorm(z/sz, lower.tail = TRUE))
}
