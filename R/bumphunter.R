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


##' bumphuntingEngine
##' 
##' @title bumphuntingEngine
##' @description bump hunting algorithm for 450k methylation microarray
##' @param dat.m 450k methylation microarray matrix via preprocessing; m x n matrix
##' @param design design matrix; [1, covariate of interested, paired info]
##' @param sv.m surrogate variables matrix 
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
##' @param robust logic, use robust linear regression or not; Default FALSE
##' @param smooth logic. If TRUE then the standard error or correlation of point-wise estimatrs will be used as weigths 
##'        in the \code{loessByCluster}
##' @param smoother smooth function to estimate genome profile.
##' @param B integer; Denoting the number of resamples to computr null distribution. Default = 0
##'        B is used to generate permutation matrix, 
##'        which describes permutations to generate null distribution
##' @param cor logical; If TRUE then position correlation matrix is considered as weighted matrix
##' @param corFunc Optional; "spearman"(Default) and "pearson"
##' @param combp logical; If TRUE then slk p correction will be applied
##' @param verbose logical. Optional printing progress message or not
##' @param ...
##' @return bumps object
##' @importFrom limma lmFit eBayes 
##' @examples 
##' # cluster input format: cg12045430 cg20826792 cg00381604 cg20253340 
##  #                                3          3          3          4 
##  # -like vector
##' @export
bumphuntingEngine <- function(dat.m = NULL, design, sv.m = NULL, chr, pos, cluster = NULL, coef = 2,
                              cutoff = NULL, pvalue = 0.01, maxGap = 500, minDist = 500,
                              nullMethod = c("permutation", "bootstrap"), robust = FALSE,
                              smooth = FALSE, smoother = bumphunter::loessByCluster,
                              B = 10000, cor = FALSE, 
                              corFunc = c("spearman", "pearson"), combp = FALSE,
                              verbose = TRUE, ...){
  nullMethod <- match.arg(nullMethod)
  mod        <- cbind(design, sv.m)
  # make cluster
  if(is.null(cluster))
    cluster <- clusterMaker(chr = chr, pos = pos, maxGap = maxGap, names = rownames(dat.m))
  
  # estimate each position coefficient profile
  # smooth means need 1/sigma as weight or not
  # robust means squeeze variance or not
  # combp  means calculate p-value or not
  if(smooth){
    if(!robust){
      tmp   <- mlm.fit(dat.m = dat.m, design = mod, coef = 2, full = TRUE)
      beta  <- tmp$coef
      sigma <- tmp$sigma
      if(combp){
        t   <- beta0 / tmp$stdev_unscaled / sigma0
        p   <- pt(t, tmp$df.residuals)
      }
    }else{
      tmp    <- lmFit(dat.m, mod)
      contrasts <- cbind("C-N" = c(0, 1, rep(0, ncol(mod) - 2)))
      tmp    <- contrasts.fit(tmp, contrasts)
      tmp    <- eBayes(tmp)
      beta0  <- tmp$coefficients
      sigma0 <- sqrt(tmp$s2.post)
      if(combp)
        p   <- tmp$p.value
    }
    rm(tmp)
  }else{
    if(!robust){
      tmp <- mlm.fit(dat.m = dat.m, design = mod, coef = 2, full = TRUE)
      beta  <- tmp$coef
      sigma <- NULL
      if(combp){
        t    <- beta0 / tmp$stdev_unscaled / sigma0
        p   <- pt(t, tmp$df.residuals)
      }
    }else{
      tmp    <- lmFit(dat.m, mod)
      contrasts <- cbind("C-N" = c(0, 1, rep(0, ncol(mod) - 2)))
      tmp    <- contrasts.fit(tmp, contrasts)
      tmp    <- eBayes(tmp)
      beta  <- tmp$coefficients
      sigma <- NULL
      if(combp)
        p   <- tmp$p.value
    }
    rm(tmp)
  }
  
  # Region search is based on :
  #           + combination-p method 
  #           + permutation method
  if(!combp && B > 0){
    # permutation calculation
    if(nullMethod == "permutation"){
      tmp <- mlm.fit(dat.m = dat.m, design = mod, coef = 2, B = B, full = TRUE)
      beta0 <- tmp$coef
      sigma0 <- tmp$sigma 
      if(robust){
        tmpr <- mlm.tstat(tmp)
        sigma0 <- tmpr$post.sigma
        rm(tmpr)
      }
      rm(tmp)
    }
    # bootstrap perform 100x speedup than permutation
    else if(nullMethod == "bootstrap"){
      tmp <- bootstrap.fit(dat.m = dat.m, design = mod, coef = 2, B = B)
      beta0 <- tmp$coef
      sigma0 <- tmp$sigma
      if(robust){
        tmpr <- mlm.tstat(tmp)
        sigma0 <- tmpr$post.sigma
        rm(tmpr)
      }
      rm(tmp)
    }
    # smooth processing
    # regionSeeker by a soft threshold and their null hypothesis
    region <- regionSeeker(beta = beta, chr = chr, pos = pos, maxGap = maxGap, names = rownames(dat.m), drop = TRUE, permbeta = beta0)
  }
}