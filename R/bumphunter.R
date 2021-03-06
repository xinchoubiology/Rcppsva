#' @title show
#' @rdname show-methods
#' @exportMethod show
setMethod("show", signature(object = "bumps"),
          function(object){
            cat(sprintf("'bumps' object with %d DMR bumps\n", nrow(object@bumps$diff)))
            cat(sprintf("'bumps' object with %d NULL H bumps for discriminative analysis\n", 
                        nrow(object@bumps$null)))
            cat(sprintf("algorithms : ...... \n"))
            for(params in names(object@algorithm)){
              cat(sprintf("\t%s = %s \n", params, object@algorithm[[params]]))
            }
          })

#' @title get.bumps
#' @rdname get.bumps-methods
#' @exportMethod get.bumps
setMethod("get.bumps", signature(object = "bumps"),
          function(object, ...){
            return(object@bumps)
          })

#' bumphuntingEngine
#' 
#' @title bumphuntingEngine
#' @description bump hunting algorithm for 450k methylation microarray
#' @param dat.m 450k methylation microarray matrix via preprocessing; m x n matrix
#' @param design design matrix; [1, covariate of interested, paired info]
#' @param sv.m surrogate variables matrix 
#' @param chr chromosome vector for probes
#' @param pos position vector for probes
#' @param cluster The clusters of locations that are to be analyzed together. via correlation
#'        perspevtive, we can use agglomerate cluster index to bump analysis. If cluster is not
#'        available, clusterMaker can be used to cluster nearby locations
#' @param coef Integer. covariate of interest's column
#' @param names probe names vector ; used in function \link{clusterMaker} and \link{regionSeeker}
#' @param cutoff numeric value. Value of estimates coefficient of covariate of interested above cutoff
#'        or below the negative of cutoff will be used as candidate bump regions. 
#' @param qvalue numeric value. cutoff of pvalue for candidate regions selection || 
#'                              qvalue cutoff for combination pvalue method
#' @param maxGap if cluster is not availabel. maxGap is used by clusterMaker to define cluster
#' @param minDist if clusters are build, the mininal distance between clusters is setted as constraint
#' @param nullMethod Method for generating null candidate regions. If ncol(design) > 2. bootstrap method is recommanded
#' @param robust logic, use robust linear regression or not; Default FALSE
#' @param smooth logic. If TRUE then the standard error or correlation of point-wise estimatrs will be used as weigths 
#'        in the \link{smoother}
#' @param smoothMethod local regression method used in \link{smoother}
#' @param B integer; Denoting the number of resamples to computr null distribution. Default = 0
#'        B is used to generate permutation matrix, 
#'        which describes permutations to generate null distribution
#' @param corr logical; If TRUE then position correlation matrix is considered as weighted matrix
#'        and clusters modified by correlation constraint; Default FALSE
#' @param corFunc Optional; "spearman"(Default) and "pearson"
#' @param combp logical; If TRUE then slk p correction will be applied
#' @param merge how to merge two sub-clusters; c("single", "complete", "average")
#' @param corr.cutoff correlation cutoff in merge algorithm
#' @param combine combine pvalue method; c("stouffer_liptak", "zscore")
#' @param verbose logical. Optional printing progress message or not
#' @param ...
#' @return bumps object
#' @importFrom limma lmFit eBayes 
#' @examples 
#' # cluster input format: cg12045430 cg20826792 cg00381604 cg20253340 
##  #                                3          3          3          4 
##  # -like vector
#' @export
bumphuntingEngine <- function(dat.m = NULL, design, sv.m = NULL, 
                              chr, pos, cluster = NULL, coef = 2,
                              names, cutoff = NULL, qvalue = 0.1, maxGap = 500, 
                              minDist = 500, robust = FALSE, smooth = FALSE, 
                              smoothMethod = c("weightedLowess", "loess", "locfit"),
                              nullMethod = c("permutation", "bootstrap"), B = 10000, corr = FALSE, 
                              corFunc = c("spearman", "pearson", "kendall"), combp = FALSE,
                              merge = c("single", "complete", "average"), corr.cutoff = 0.8,
                              combine = c("stouffer_liptak", "zscore"),
                              verbose = TRUE, ...){
  
  nullMethod <- match.arg(nullMethod)
  corFunc    <- match.arg(corFunc)
  merge      <- match.arg(merge)
  combine    <- match.arg(combine)
  smoothMethod <- match.arg(smoothMethod)
  ## model matrix
  mod        <- cbind(design, sv.m)
  ## make cluster
  if(verbose){
    cat(sprintf("[Bumphunting]\t Generating clusters by gap constraint = %d \n", maxGap))
  }
  if(is.null(cluster)){
    cluster <- clusterMaker(chr = chr, pos = pos, maxGap = maxGap, names = names)
    # make correlated cluster or not
    if(corr){
      if(verbose){
        cat(sprintf("[Bumphunting]\t Splitting clusters whose correlation >= %f \n", cor.cutoff))
      }
      cluster <- corrclusterMaker(dat.m = dat.m, chr = chr, pos = pos, names = names, cluster = cluster, cutoff = cor.cutoff, maxGap = maxGap, method = corFunc, merge = merge)
    }
  }
  # estimate each position coefficient profile
  # robust means squeeze variance or not
  # combp  means calculate p-value or not
  if(verbose){
    cat(sprintf("[Bumphunting]\t Estimating beta value for each probes \n"))
  }
  if(!robust){
    tmp   <- mlm.fit(dat.m = dat.m, design = mod, coef = 2, full = TRUE)
    beta  <- tmp$coef
    sigma <- tmp$sigma
    if(combp){
      t   <- beta / sigma
      p   <- 2 * pmin(pt(t, tmp$df.residuals), 1 - pt(t, tmp$df.residuals))
    }
    rm(tmp)
  }else{
    tmp   <- lmFit(dat.m, mod)
    contrasts <- cbind("C-N" = c(0, 1, rep(0, ncol(mod) - 2)))
    tmp   <- contrasts.fit(tmp, contrasts)
    tmp   <- eBayes(tmp)
    beta  <- tmp$coefficients
    sigma <- sqrt(tmp$s2.post)
    if(combp)
      p   <- tmp$p.value
    rm(tmp)
  }
  weight <- NULL
  if(smooth){
    weight <- sigma
    if(verbose){
      cat(sprintf("[Bumphunting]\t Smoothing beta values within cluster \n"))
    }
    beta <- smoother(beta = beta, pos = pos, names = names, cluster = cluster, weight = weight, method = smoothMethod)
  }
  # Region search is based on :
  #           + combination-p method 
  #           + permutation method
  if(!combp && B > 0){
    if(verbose){
      cat(sprintf("[Bumphunting]\t Computing null distribution of beta by %s \n", nullMethod))
    }
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
    # regionSeeker by a soft threshold and their null hypothesis
    if(verbose){
      cat(sprintf("[Bumphunting]\t Finding DMRs... \n"))
    }
    region <- regionSeeker(beta = beta, chr = chr, pos = pos, names = names, cluster = cluster, maxGap = maxGap, drop = TRUE, permbeta = beta0, corr = corr, qvalue = qvalue)
  } else if(combp){
    if(verbose){
      cat(sprintf("[Bumphunting]\t Finding DMRs by Comb-p method... \n"))
    }
    ## TODO fixed weight argument error
    weight <- NULL
    region <- combine.pvalue(dat.m = dat.m, pvalues = p, cluster = cluster, chr = chr, pos = pos, names = names, method = method, combine = combine, weight = weight, cutoff = qvalue)
  }
  
  algorithm <- list(method = if(combp) combine else nullMethod, 
                    smooth = if(smooth) "smoothMethod" else "NULL",
                    correlation = if(corr) corFunc else "NULL",
                    corrcutoff = if(corr) corr.cutoff else "NULL",
                    bumpcitoff = comb.cutoff, cluster = merge, qvalue = qvalue)
  
  new("bumps", algorithm = algorithm, bumps = region)
}