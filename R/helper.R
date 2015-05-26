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
##' @param verbose Optional; Default FALSE
##' @return p vector of F-statistic p-values for each row of dat
##' @export
##' @author Xin Zhou
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




