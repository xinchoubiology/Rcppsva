## Following four find empirical hyper-prior values
aprior <- function(gamma.hat){
  m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2
}

bprior <- function(gamma.hat){
  m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2
}

postmean <- function(g.hat,g.bar,n,d.star,t2){
  (t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)
}

postvar <- function(sum2,n,a,b){
  (.5*sum2+b)/(n/2+a-1)
}

# Pass in entire data set, the design matrix for the entire data, the batch means,
# the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch.
# Uses the EM to find the parametric batch adjustments
it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- apply(!is.na(sdat),1,sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    d.new <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  #cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

#' Monte Carlo integration functions
#' @param sdata standard data matrix
#' @param g.hat estimated gamma values
#' @param d.hat estimated delta values
#' @return adjust data frame
#'         g.star
#'         d.star
int.eprior <- function(sdat,g.hat,d.hat){
  adjust <-  do.call(rbind, int_eprior(sdat, g.hat, d.hat))
  adjust
}

#' Quick calculate f statistic p-values for mod and mod0
#' @title pvalue
#' @param dat.m data matrix for expression/methylation microarray
#' @param mod the model being used to fit the data
#' @param mod0 the null hypothesis model being compared when fitting the data
#' @param verbose Optional; Default FALSE
#' @return p vector of F-statistic p-values for each row of dat
#' @export
#' @author Xin Zhou
pvalue <- function(dat.m = NULL, mod = NULL, mod0 = NULL, ...){
  df1 <- ncol(dat.m) - dim(mod)[2]
  df0 <- ncol(dat.m) - dim(mod0)[2]
  Id <- diag(ncol(dat.m))
  
  ## resid <- dat.m %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid <- dat.m - dat.m %*% mod %*% tcrossprod(solve(crossprod(mod)), mod)
  rss1  <- rowSums(resid*resid)
  
  resid0 <- dat.m - dat.m %*% mod0 %*% tcrossprod(solve(crossprod(mod0)), mod0)
  rss0   <- rowSums(resid0*resid0)
  
  fstats <- ((rss0 - rss1)/(df0 - df1))/(rss1 / df1)
  p <- 1 - pf(fstats, df1 = (df0 - df1), df2 = df1)
  p
}

#' estimate coefficient of interest variables with known Svs, \link{mlm.tstat}
#' 
#' @title mlm.fit
#' @param dat.m n x m matrix of methylation microarray
#' @param design design matrix for expression data matrix(data.m)
#' @param coef covariate of interest; Default = 2
#' @param B permutation number of covariates of interest
#' @param full return full regression object; Optional FALSE
#' @return list
#'         coefficient
#'         stdev_unscale
#'         sigma
#'         df.residule
#' @export
#' @importFrom doMC registerDoMC
#' @examples
#' sd <- 0.3 * sqrt(4/rchisq(100, df = 4))
#' y  <- matrix(rnorm(100*6, sd = sd), 100, 6)  # each row of data is generate by sd[i] ~ invchisq
#' rownames(y) <- paste("cg", 1:100)
#' # introduce 2 cgs which are DMPs 
#' y[1:2, 4:6] <- y[1:2, 4:6] + 2 # have significant differential when we introduce the poi(cancer - normal)
#' pheno <- factor(c(0,0,0,1,1,1))
#' levels(pheno) <- c("normal", "cancer")
#' design <- model.matrix(~ pheno)
#' fit <- mlm.fit(y, design)
#' sig.tab <- mlm.tstat(fit)
#' library(doMC)
#' registerDoMC(2)
mlm.fit <- function(dat.m = NULL, design = NULL, coef = 2, B = NULL, full = FALSE, mcore = 4){
  x0 <- design[,coef, drop = FALSE]
  xb <- design[,-coef, drop = FALSE]
  xperm <- if(is.null(B)){
    x0
  }else{
    ## replicate(B, sample(x0))
    replicate(B, unlist(lapply(1:(nrow(x0)/2), function(x) sample(c(1,0)))))
  }
  result <- beta_regress(M = dat.m, pv = xperm, svs = xb, full = as.numeric(TRUE))
  if(!is.null(rownames(dat.m))){
    rownames(result$coef) <- rownames(result$sigma) <- rownames(dat.m)
  }
  ## normalized by stdev_unscaled
  if(!is.null(B) && B >= 2){
    result$coef <- t(apply(result$coef, 1, '/', t(result$stdev_unscaled)))
  }else{
    result$coef <- result$coef / result$stdev_unscaled[1]
  }
  result$stdev_unscaled <- NULL
  result
}

#' bootstrap wrapper
#' 
#' @title bootstrap.fit
#' @param dat.m n x m matrix of methylation microarray
#' @param design design matrix for expression data matrix(data.m)
#' @param coef covariate of interest; Default = 2
#' @param B permutation number of covariates of interest
#' @return list
#'         coefficient
#'         sigma
#'         df.residule
#' @export
bootstrap.fit <- function(dat.m = NULL, design = NULL, coef = 2, B = NULL){
  mod   <- design
  mod0  <- design[,-coef,drop=FALSE]
  xperm <- if(is.null(B)){
    matrix(1:ncol(dat.m), ncol = 1)
  }else{
    replicate(B, sample(1:ncol(dat.m), replace = TRUE), simplify=TRUE) - 1
  }
  result <- bootstrap_regress(M = dat.m, mod = mod, modn = mod0, B = xperm)
  if(!is.null(rownames(dat.m))){
    rownames(result$coef) <- rownames(result$sigma) <- rownames(dat.m)
  }
  result
}

#' moderest lm T statistic
#' @title mlm.tstat
#' @param fit object of mlm.fit; contains coefficient, sigma
#' @return list of 2 components
#'         post.sigma
#'         df.total
#'         t score
#' @export
mlm.tstat <- function(fit){
  s2 <- apply(fit$sigma^2, 2, limma::squeezeVar, fit$df.residuals)
  B  <- ncol(fit$coef)
  t  <- lapply(1:B, function(i) fit$coef[,i] / sqrt(s2[[i]]$var.post))
  
  df.total <- fit$df.residuals + sapply(s2, "[[", "df.prior")
  post.sigma <- sqrt(do.call(cbind, lapply(s2, "[[", "var.post")))
  ## p.value <- lapply(1:B, function(i) 2 * pt(-abs(t[[i]]), df = df.total[[i]]))
  
  return(list(post.sigma = post.sigma, df.total = df.total))
}

#' If our data is paired, wilcoxon rank test can be used to test the significnat of sites
#' 
#' @title wilcox.fit
#' @description calculate wolcoxon signed rank sum test for each CpG probes of microarray;
#'              support multiple core parallism calculation.
#' @param dat.m n x m delta M|beta matrix for n CpG sites across 2*m paired different patient samples
#' @param alternative a character string specifying the alternative hypothesis; 
#'        c("two.sided", "greater", "less") and "two.sided" is default
#' @param qvalue0 false discovery rate's threshold; Default = 0.1
#' @param mcore Cpu cores can be used in test
#' @return res list
#'             siggene data.frame
#'             null    data.frame
#' @importFrom doMC registerDoMC
#' @importFrom plyr aaply
#' @importFrom qvalue qvalue
#' @export
#' @author Xin Zhou
wilcox.fit <- function(dat.m = NULL, alternative = c("two.sided", "greater", "less"), qvalue0 = 0.1, mcore = 4){
  alternative <- match.arg(alternative)
  registerDoMC(cores = mcore)
  if(ncol(dat.m) >= 10){
    rank.m <- aaply(dat.m, 1, .fun = function(d){
      r <- rep(0, length(d))
      r[which(d != 0)] <- rank(d[which(d != 0)])
    }, .parallel = TRUE
    )
    sgn.m  <- aaply(dat.m, 1, sign, .parallel = TRUE)
    W      <- abs(rowSums(rank.m * sgn.m))
    Nr     <- rowSums(abs(sgn.m))
    sigma  <- sqrt(Nr * (Nr + 1) * (2 * Nr + 1)/6)
    p <- pnorm(W, mean = 0.5, sd = sigma, lower.tail = TRUE)
    pval <- switch(alternative,
                   two.sided = pmin(p, 1-p),
                   greater   = 1-p,
                   less      = p
    )
    qval <- p.adjust(pval, method = "fdr")
    sig.cg <- which(qval <= qvalue0)
  } else {
    pval <- aaply(dat.m, 1, .fun = function(x){
      wilcox.test(x, alternative = alternative)$p.value
    }, .parallel = TRUE
    )
    qval <- p.adjust(pval, method = "fdr")
    sig.cg <- which(qval <= qvalue0)
  }
  
  table <- data.frame(pvalue = pval, qvalue = qval)
  rownames(table) <- rownames(dat.m)
  
  res <- list(siggene = table[sig.cg,], null = table[-sig.cg,])
  res
}

#' clusterMaker function based on spatial distance
#'
#' @title clusterMaker
#' @param chr Chromosome vector
#' @param pos position numeric vector
#' @param maxGap max gap length between 2 probes within a region
#' @param names probe names vector
#' @return cluster data.frame of length m(m probes). Probes within same region have same region identities list.
#' @importFrom plyr llply
#' @import GenomicRanges
#' @details If you define the maxGap 500, so we can guarantee probes in different regions has distance >= maxGap.
#'          As the result, if we extend region upstream & downstream maxGap / 2, our extended regions have no overlap
#' @examples 
#' require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' Location <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@@data$Locations)
#' Probes <- rownames(Location)[rownames(Location) %in% rownames(dat.m)]
#' BED <- cbind(Location[Probes, 1:2])
#' BED <- cbind(BED, dat.m[Probes, 4:3])
#' cluster <- clusterMaker(chr = BED$chr, pos = BED$pos, maxGap = 1000, names = rownames(BED))
#' @export
clusterMaker <- function(chr, pos, maxGap = 500, names){
  genome     <- GRanges(seqnames = chr, ranges = IRanges::IRanges(start = pos, width = 1), names = names)
  seqlevels(genome) <- sort(seqlevels(genome))
  genome     <- sort(genome)
  chrIndexes <- split(genome, seqnames(genome))
  cluster    <- llply(chrIndexes, function(chr) cumsum((diff(c(0, start(chr))) > maxGap) + 0))
  offset     <- 0
  for(i in 1:length(cluster)){
    cluster[[i]] <- cluster[[i]] + offset
    offset <- tail(cluster[[i]], 1)
  }
  cluster <- unlist(cluster)
  names(cluster) <- genome$names
  cluster
}

#' correlated cluster maker. Build the CpG cluster by correlation
#' 
#' @description By \code{clusterMaker}, we can generate probe cluster under constraint of maximum 
#'              distance. Adding the correlation constraint, we will truncate exsited cluster and
#'              calculate their correlation matrix \code{sigma} out
#' @title corrclusterMaker
#' @param dat.m n x m delta M|beta matrix for n CpG sites across 2*m paired different patient samples
#' @param chr Chromosome vector
#' @param pos position numeric vector
#' @param cluster vector of probes clusters By distance(maxGap) constraints
#' @param names probe names vector ; used in function \link{clusterMaker}
#' @param maxGap max gap length between 2 probes within a region ; Default = 500
#'        also used in function \link{clusterMaker}
#' @param cutoff threshold to filter distance based cluster; Default = 0.8
#' @param method correlation calculation method; c("spearman", "pearson", "kendall"); Default "spearman"
#' @param merge how to merge two sub-clusters; c("single", "complete", "average")
#' @param pos CpG position with their probes id
#' @param corrmat correlation matrix for each cluster; Default NULL
#' @param cores number of thread used by \link{corrclusterMaker}
#' @details only clusters contain >= 2 probes can get a corrected p-value
#' @return list of clusterID vector
#' @importFrom plyr llply
#' @importFrom doMC registerDoMC
#' @importFrom parallel detectCores
#' @export
corrclusterMaker <- function(dat.m = NULL, chr, pos, cluster = NULL, 
                             maxGap = 500, names, cutoff = 0.8, 
                             corrmat = NULL, cores = detectCores(),
                             method = c("pearson", "spearman", "kendall"), 
                             merge = c("single", "complete", "average")){
  method <- match.arg(method)
  merge  <- match.arg(merge)
  
  if(is.null(cluster)){
    cluster <- clusterMaker(chr = chr, pos = pos, maxGap = maxGap, names = names)
  }
  cnames      <- names(cluster)
  names(pos)  <- names
  rawcluster  <- split(cnames, cluster)
  combIndex   <- which(sapply(rawcluster, function(c) length(c) > 1))
  multcluster <- rawcluster[combIndex]
  
  ## build correlation matrix within clusters
  if(is.null(corrmat)){
    if(cores >= 2){
      registerDoMC(cores = cores)
      corrmat    <- llply(multcluster, .fun = function(ix){
                                             dist <- abs(outer(pos[ix], pos[ix], "-"))
                                             # colnames(dist) <- rownames(dist) <- ix
                                             cor(t(dat.m[ix,]), method = method) * (dist < maxGap)
                                           }, .parallel = TRUE)
    }else{
      corrmat    <- llply(multcluster, .fun = function(ix){
                                             dist <- abs(outer(pos[ix], pos[ix], "-"))
                                             colnames(dist) <- rownames(dist) <- ix
                                             cor(t(dat.m[ix,]), method = method) * (dist < maxGap)
                                           })
    }
  }
  corrcluster    <- Dbpmerge(corrmat, merge, cutoff)
  cluster        <- c(corrcluster, rawcluster[-combIndex])
  names(cluster) <- seq(1:length(cluster))
  
  cluster
}

#' sub-regions merge and split
#' 
#' @title Dbpmerge
#' @param c.mat list of correlation probes matrices
#' @param merge merge method; c("single", "complete", "average")
#' @param cutoff correlation >= cutoff can be merge into a sub-regions
#' @return list of all sub-clusters
#' @importFrom plyr llply
#' @export
Dbpmerge <- function(c.mat = NULL, merge = c("single", "complete", "average"), cutoff = 0.8){
  merge <- match.arg(merge)
  corrcluster <- llply(c.mat, .fun = function(mx){
                                    pname <- rownames(mx)
                                    diag(mx) <- 0
                                    mx <- switch(merge, single   = single_linkage(mx),
                                                        complete = complete_linkage(mx),
                                                        average  = average_linkage(mx))
                                    mx <- (mx <= cutoff) + 0.0
                                    cIndexes <- cumsum(c(1, clique_merge(mx)))
                                    split(pname, cIndexes)
                                   })
  unlist(corrcluster, recursive = FALSE)
}



#' segment vectors into positive, null and negative
#' 
#' @title segmentsMaker
#' @param cluster vector of probes clusters
#' @param beta beta value for probes
#' @param chr chromosome vector
#' @param pos position vector
#' @param cutoff cutoff value for positive & negative value; 
#'        <= cutoff[1] : hypo 
#'        >= cutoff[4] : hyper
#'        (cutoff[2], cutoff[3]) : null
#' @param permbeta null hypothesis distribution
#' @return cluster list and each cluster contains 3 list
#'         `hyper`| `+`
#'         `hypo` | `-`
#'         `null` | `0`
#' @details bumphunting kernel; and we can use comb-p method in segmentMaker.
#' @import data.table
#' @export
segmentsMaker <- function(cluster, beta, cutoff, chr, pos, permbeta){
  cnames <- names(cluster)
  all <- ((beta >= cutoff[4]) - (beta <= cutoff[1]) + 2 * (beta >= cutoff[2]) * (beta <= cutoff[3]))[cnames,]
  sgn <- as.numeric(c(0, abs(diff(as.numeric(all)))) != 0)
  splitter <- cumsum(sgn) + cluster
  all[all == 1]  <- "hyper"
  all[all == 2]  <- "null"
  all[all == -1] <- "hypo"
  all[all == 0]  <- "unknown"
  perm <- permbeta[cnames, ]
  colnames(perm) <- paste0("perm", 1:ncol(permbeta))
  
  segments <- data.table(names = cnames, cluster = cluster, beta = beta[cnames,], status = all, segs = splitter, chr = chr, pos = pos, perm)
  segments <- segments[segments$status != "unknown", ]
  segments <- plyr:::splitter_d(segments, .(status))
  segments <- list(hyper = plyr:::splitter_d(segments[[1]], .(segs)),
                   hypo  = plyr:::splitter_d(segments[[2]], .(segs)),
                   null  = plyr:::splitter_d(segments[[3]], .(segs)))
  segments
}

#' segmentsCluster clustering correlated probes
#' 
#' @title segmentsCluster
#' @param cluster vector of probes clusters
#' @param beta beta value for probes
#' @param chr chromosome vector
#' @param pos position vector
#' @param permbeta null hypothesis distribution
#' @return data.table of all clusters
#' @export
segmentsCluster <- function(cluster, beta, chr, pos, permbeta){
  cnames <- names(cluster)
  perm   <- permbeta[cnames, ]
  segments <- data.table(names = cnames, cluster = cluster, beta = beta[cnames,], chr = chr, pos = pos, perm)
  segments <- plyr:::splitter_d(segments, .(cluster))
  segments
}

#' bump hunting algorithm one clusters predefined by distance/correlation constraints
#' 
#' @title regionSeeker
#' @param beta vector of different probes
#' @param chr Chromosome vector
#' @param pos position numeric vector
#' @param pvalc pvalue cutoff for bump selection; Default 0.1
#' @param cluster list clusters defined by arguments or \link{clusterMaker}
#' @param maxGap max gap length between 2 probes within a region ; 
#'        used in function \link{clusterMaker}
#' @param names probe names vector ; used in function \link{clusterMaker}
#' @param cutoff cutoff threshold for bump selection. only segements within cluster satisfy cutoff
#'        are returned as predicted region c(sig-cutoff, null-cutoff)
#' @param permbeta Null hypothesis beta distribution 
#'        [added permbeta(H0 distribution) to the end columns of region table from \link{regionSeeker}]
#' @param drop FALSE; create discriminate table (significant | null) or not
#' @param corr logical, whether the cluster is generated under distance constraint or 
#'        under correlation constraint. Default FALSE
#' @param mcores multiple threads used in \link{regionSeeker}
#' @param qvalue cutoff of FDR(false discovery rate)
#' @param verbose FALSE
#' @details If an arbitary threshold is defined, regionSeeker will return a table (within / without)
#'          the threshold. In bump hunting algorithms, these contiguous probes mean bumps.
#' @import GenomicRanges
#' @import data.table
#' @return Table of predict regions
#' @export
regionSeeker <- function(beta, chr, pos, cluster = NULL, maxGap = 500, names, pvalc = 0.1,
                         cutoff = c(quantile(abs(beta), 1-pvalc), quantile(abs(beta), pvalc)),
                         permbeta = NULL, mcores = 2, corr = FALSE, qvalue = 0.1,
                         drop = FALSE, verbose = FALSE){
  if(is.null(cluster)){
    cluster <- clusterMaker(chr = chr, pos = pos, maxGap = maxGap, names = names)
  }
  if(inherits(cluster, "list")){
    cluster <- List2Index(cluster = cluster, chr = chr, pos = pos, names = names)
  }
  if(mcores >= 2){
    registerDoMC(cores = mcores)
  }
  
  genome     <- GRanges(seqnames = chr, ranges = IRanges::IRanges(start = pos, width = 1), names = names, chr = chr)
  seqlevels(genome) <- sort(seqlevels(genome))
  genome     <- sort(genome)
  L          <- ncol(permbeta)
  
  if(corr){
    segments <- segmentsCluster(cluster = cluster, beta = beta,
                                chr = genome$chr,
                                pos = start(genome),
                                permbeta = permbeta)
    
    out      <- llply(segments, function(ix){
                                  data.table(chr     = ix[1,4],
                                             start   = min(ix[,5]), 
                                             end     = max(ix[,5]) + 1,
                                             length  = max(ix[,5]) - min(ix[,5]) + 1,
                                             value   = mean(ix[,3]), 
                                             area    = abs(sum(ix[,3])), 
                                             cluster = ix[1,2], 
                                             L = nrow(ix),
                                             cgnames = paste0(ix[,1],collapse = ";"),
                                             pvalue  = min(sum(mean(ix[,3]) >= colMeans(ix[,-(1:5)])), sum(mean(ix[,3]) <= colMeans(ix[,-(1:5)]))) / L,
                                             parea   = sum(abs(sum(ix[,3])) <= abs(colSums(ix[,-(1:5)])))/ L)
                                }, .parallel = TRUE)
    
    res  <- do.call(rbind, out)
    qval <- p.adjust(res$parea, method = "fdr")
    res$qvalue <- qval
    
    diff <- res[(res$value <= -cutoff[1] || res$value >= cutoff[1]),]
    null <- res[(res$value >= -cutoff[2] && res$value <= cutoff[2]),]
    return(list(diff = diff[diff$qvalue <= qvalue,], null = null[null$qvalue > qvalue,]))
  }else{
    segments <- segmentsMaker(cluster = cluster, beta = beta,
                              cutoff = c(-cutoff[1], -cutoff[2], cutoff[2], cutoff[1]),
                              chr = genome$chr,
                              pos = start(genome), 
                              permbeta = permbeta)
    res <- vector("list", 3)
    
    if(verbose)
      cat("[regionSeeker]\t Estimating Statistics For Regions Hyper, Hypo, NULL ...")
    
    for(i in 1:3){
      out <- llply(segments[[i]], function(ix){
                                    data.table(chr     = ix[1,6], 
                                               start   = min(ix[,7]), 
                                               end     = max(ix[,7]) + 1, 
                                               length  = max(ix[,7]) - min(ix[,7]) + 1,
                                               value   = mean(ix[,3]), 
                                               area    = abs(sum(ix[,3])), 
                                               cluster = ix[1,2], 
                                               L = nrow(ix), 
                                               cgnames = paste0(ix[,1],collapse = ";"),
                                               pvalue  = min(sum(mean(ix[,3]) >= colMeans(ix[,-(1:7)])), sum(mean(ix[,3]) <= colMeans(ix[,-(1:7)]))) / L,
                                               parea   = sum(abs(sum(ix[,3])) <= abs(colSums(ix[,-(1:7)])))/ L,
                                               status  = names(segments)[i])
                                  }, .parallel = TRUE)
      res[[i]] <- do.call(rbind, out)
    }
    
    names(res) <- names(segments)
    if(drop){
      res <- list(diff = rbind(res$hyper, res$hypo), null = res$null)
    }
    return(res)
  }
}

#' fitting L/S model with missing values
#' @description Beta.NA function is borrow from package sva
#' @title Beta.NA
Beta.NA <- function(y,X){
  des <- X[!is.na(y),]
  y1 <- y[!is.na(y)]
  B <- solve(t(des)%*%des)%*%t(des)%*%y1
  B
}

#' get contains NA's row, col
#' @description Index.NA function is borrow from package sva
#' @title Index.NA
Index.NA <- function(mat, by = c("row", "col")){
  by <- match.arg(by)
  if(by == "row"){
    unique(which(is.na(mat)) %% nrow(mat))
  }
  else{
    unique(which(is.na(mat)) %% ncol(mat))
  }
}

#' convert cluster list to Index with increased number
#' @title List2Index
#' @param cluster vector of probes clusters By distance(maxGap) constraints
#' @param chr Chromosome vector
#' @param pos position numeric vector
#' @param names probe names vector ; used in function \link{clusterMaker}
#' @return clusterIndex
List2Index <- function(cluster, chr, pos, names){
  ## prepare IRange
  genome     <- GRanges(seqnames = chr, ranges = IRanges::IRanges(start = pos, width = 1), names = names)
  seqlevels(genome) <- sort(seqlevels(genome))
  genome     <- sort(genome)
  
  ## convert from cluster to index
  ## create cluster vector
  indexes <- unlist(sapply(1:length(cluster), function(ix) rep(ix, length(cluster[[ix]]))))
  names(indexes) <- unlist(cluster)
  cluster <- indexes[genome$names]
  indexes <- seq(1:length(unique(cluster)))
  names(indexes) <- unique(cluster)
  cluster <- indexes[as.character(cluster)]
  names(cluster) <- genome$names
  
  cluster
}

#' execute hierachical clustering by beta value data directly
#' @title hclust_from_data
#' @param data beta value matrix
#' @param link character for linkage type[default "ward" out of 
#'        c("ward", "average", "complete", "single")]
#' @param dist character for distance type[default "euclidean" out of 
#'        c("euclidean", "manhanttan")]
#' @param minkowski soft_power of minkowski distance
#' @return hclust object for our cluster
#' @export
hclust_from_data <- function(data, link, dist, minkowski){
  .Call('Rcppsva_hclust_from_data', 
        data, 
        as.integer(link), 
        as.integer(dist),
        as.numeric(minkowski), DUP = FALSE, NAOK = FALSE, PACKAGE = 'Rcppsva')
}

#' Hierachical Clustering function based on fast parallel hierachical clustering
#' @title HClust
#' @param data matrix of beta data
#' @param method The agglomeration method to be used. This must be one of "ward", "single", "complete" or "average"
#' @param distance The distance measure to be used. This must be one of "euclidiean", "manhattan", "maximum", or "minkowski".
#' @param p power of the Minkowski distance.
#' @param sign Optional 'U' for unsigned and 'S' for signed dissimilarity
#' @return An object of class *hclust* which describes the tree produced by the clustering process
#' @export
#' @author Xin Zhou \url{xxz220@@miami.edu}
HClust <- function(data = NULL, method = "average", distance = "euclidean", p = 2, sign = c("S", "U")){
  options(warn = -1)
  if(method == "ward" && distance != "euclidean"){
    warning("Distance method is forced to (squared) 'euclidean' distance for Ward's method")
    distance = "euclidean"
  }
  if(is.null(data)){
    stop("Beta matrix is needed for Chclust")
  }
  if(distance == "spearman"){
    data <- rankm(as.matrix(data), byrow = TRUE)
    distance <- "pearson"
  }
  
  sign <- match.arg(sign)
  if(distance == "pearson"){
    if(sign == "S"){
      distance <- "scosine"
    } else{
      distance <- "ucosine"
    }
  }
  
  method <- pmatch(method, linkage_kinds())
  if(is.na(method)){
    stop("clustering method is not support by Rcppsva ...")
  }
  if(method == -1){
    stop("ambigous clustering method ...")
  }
  distance <- pmatch(distance, distance_kinds())
  if(is.na(distance)){
    stop("distance metric is not support by Rcppsva ...")
  }
  if(distance == -1){
    stop("ambigous distance metric ...")
  }
  
  ## Debug printf
  ## cat("  start hierachical clustering ... \n")
  hcl <- hclust_from_data(as.matrix(data), link = method, dist = distance, minkowski = p)
  hcl$labels      <- row.names(data)
  hcl$methods     <- linkage_kinds()[method]
  hcl$call        <- match.call()
  hcl$dist.method <- distance_kinds()[distance]
  class(hcl)      <- "hclust"
  hcl
}

# optimal cluster script is a collection for optimal cluster number determination
# optimal height selection, when given an array of candidate height, we can use gap
# statistic for cut-height determination.

#' To calculate the optimal cluster number of our data, we need a reference distribution
#' generation for NULL_DENDROGRAM
#' @title distributeRef
#' @param data  data matrix to sampled
#' @param size  bootstrap sampling size
#' @param by    data will independent by row('R')[Default] or by column
#' @param n     restrict sampled data matrix size
#' @param umethod 'uniform'(Default). Method for sampling reference null distribution.
#'               c('uniform', 'bootstrap')
#' @param mcores parallel threads to be used
#' @param verbose FALSE 
#' @return list of NULL dendrogram objects
#' @export 
#' @author Xin Zhou \url{xxz220@@miami.edu}
distributeRef <- function(data = NULL, size = 10, by = 'R', n = 1000, 
                          umethod = c('uniform', 'bootstrap', 'MCMC'), verbose = FALSE,
                          mcores = 2, ...){
  options(warn = -1)
  if(is.null(data)){
    stop("background data matrix is not supported . ")
  }
  umethod <- match.arg(umethod)
  if(umethod == "MCMC"){
    R <- svd(data)
    X <- data %*% R$v
  }
  NULL_dendrogram <- llply(1:size, function(x){
                                      if(verbose){
                                        cat(sprintf("\t simulating %d round reference distribution ...\n", x))
                                      }
                                      set.seed(x + 1001)
                                      if(mcores >= 2){
                                        registerDoMC(cores = mcores)
                                      }
                                      sdat <- do.call(cbind, 
                                                       llply(1:ncol(data), 
                                                             function(ix)
                                                             {
                                                               if(umethod == "bootstrap"){
                                                                 data[sample(1:nrow(data), n, replace = TRUE), ix, drop = TRUE]
                                                               } else if(umethod == "MCMC"){
                                                                 runif(n, min = min(X[,ix]), max = max(X[,ix]))
                                                               } else{
                                                                 # uniform distribution subtitute bootstrap
                                                                 runif(n, min = min(data[,ix]), max = max(data[,ix]))
                                                               }
                                                             }, .parallel = TRUE)
                                                     )
                                      if(umethod == "MCMC"){
                                        sdat <- tcrossprod(sdat, R$v)
                                      }
                                      if(verbose){
                                        cat(sprintf("\t hierachical clustering ... \n"))
                                      }
                                      list(hclust = HClust(sdat, ...), data = sdat)
                                   })
  NULL_dendrogram
}