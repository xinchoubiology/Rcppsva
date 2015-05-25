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

##' Monte Carlo integration functions
##' @param sdata standard data matrix
##' @param g.hat estimated gamma values
##' @param d.hat estimated delta values
##' @return adjust data frame
##'         g.star
##'         d.star
int.eprior <- function(sdat,g.hat,d.hat){
  adjust <-  do.call(rbind, int_eprior(sdat, g.hat, d.hat))
  adjust
}

##' Quick calculate f statistic p-values for mod and mod0
##' @title pvalue
##' @param dat.m data matrix for expression/methylation microarray
##' @param mod the model being used to fit the data
##' @param mod0 the null hypothesis model being compared when fitting the data
##' @param qval0 threshold of qvalue
##' @param verbose Optional; Default FALSE
##' @return p vector of F-statistic p-values for each row of dat
##' @export
##' @author Xin Zhou
pvalue <- function(dat.m = NULL, mod = NULL, mod0 = NULL, qval0 = 0.1, verbose = FALSE){
  df  <- ncol(dat.m) - dim(mod)[2] + 1
  df0 <- ncol(dat.m) - dim(mod0)[2] + 1
  p   <- rep(0, nrow(dat.m))
  res <- dat.m - dat.m %*% mod %*% tcrossprod(solve(crossprod(mod)), mod)
  rss <- rowSums(res * res)
  rm(res)
  
  res0 <- dat.m - dat.m %*% mod0 %*% tcrossprod(solve(crossprod(mod0)), mod0)
  rss0 <- rowSums(res0 * res0)
  rm(res0)
  
  f.stats <- (rss0 - rss) / (df0 - df) / rss / df
  pv.v    <- 1 - pf(f.stats, df1 = df0 - df, df2 = df)
  pv.s    <- sort(pv.v, decreasing = FALSE, index.return = TRUE)
  qv.v    <- qvalue(pv.s$x)$qvalue
  nsig    <- length(which(qv.v < qvalue0))
  if(verbose){
    cat(sprintf("%d CpG sites are significante differential in cancer ...\n", nsig))
  }
  if(nsig > 0){
    pred.idx <- pv.s$ix[1:nsig]
    lm.o     <- lm(t(dat.m) ~ mod - 1)
    t.stats  <- sapply(summary(lm.o), function(x) x$coefficients[2,3])
    res      <- cbind(t.stats, pv.s[1:nsig], qv.v[1:nsig])
    colnames(res) <- c("t-stat", "fpval", "qval")
    rownames(res) <- rownames(dat.m)[pred.idx]
  }
  else{
    res <- NULL
  }
  
  res
}




