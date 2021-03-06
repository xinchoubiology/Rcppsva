#' estimate the dimension of Methylation status methylation by RMT
#'
#' @title estDimRMT
#' @param dat.m n x m M|beta matrix for n CpG sites across m different patient samples.
#' @param plot Optional bool value, plot eigenvalue distribution or not. Default = TRUE
#' @return cor Cross correlation matrix across different samples
#' @return dim Estimate dimension of data
#' @return empeigen Empirical distribution of eigenvalue
#' @return theigen Theoretical distribution of eigenvalue
#' @return evalue Eigen values of correlation matrix
#' @return ggplot Options ggplot2 object
#' @import ggplot2
#' @export
#' @author Xin Zhou \url{xinchoubiology@@gmail.com}
estDimRMT <- function(dat.m, plot = TRUE){
  ## transform value to "M"
  ## normalize by each column
  M <- apply(dat.m, 2, function(X){(X-mean(X))/sd(X)})
  ## define random P(lambda)
  ### Q = nrow(M) / ncol(M) > 1
  sigma2 <- var(as.vector(M));
  Q <- nrow(M) / ncol(M)
  lambdaMAX <- (1 + 1/Q + 2*sqrt(1/Q))
  lambdaMIN <- (1 + 1/Q - 2*sqrt(1/Q))
  ## lambdaMAX <- sigma2 * (1 + 2*sqrt(1/Q))
  ## lambdaMIN <- sigma2 * (1 - 2*sqrt(1/Q))

  ### P(lambda)'s sampling
  ns <- ncol(M)
  delta <- lambdaMAX - lambdaMIN
  roundN <- 3
  step <- round(delta/ns, roundN)
  while(step == 0){
    roundN <- roundN + 1;
    step <- round(delta/ns, roundN); ## maximal space
  }

  ### P_{rm}(lambda) = Q/2pi * sqrt(lambdaMAX-lambda)sqrt(lambda-lambdaMIN)/lambda
  lambda.v <- seq(lambdaMIN, lambdaMAX, by = step);
  dens.v <- vector();  ## this is
  ii <- 1
  for(lambda in lambda.v){
    dens.v[ii] <- Q / (2 * pi) * sqrt(lambdaMAX - lambda) * sqrt(lambda - lambdaMIN) / lambda;
    ## dens.v[ii] <- Q / (2 * pi * sigma2) * sqrt(max(0, (lambdaMAX - lambda) * (lambda - lambdaMIN))) / lambda;
    ii <- ii + 1;
  }

  ### so the theigen can be represent by lambda and their dens.v
  theigen <- list(min = lambdaMIN, max = lambdaMAX, step = step, lambda = lambda.v, dens = dens.v)

  ### to represent different samples's cross-correlation, we use the t(M) %*% M
  C <- 1 / nrow(M) * t(M) %*% M;
  eigen.v <- sort(arma_eigen(C), decreasing = TRUE);
  ### eigens <- eigen(C, symmetric = TRUE)
  ### eigen.v <- eigens$values
  ### via eigen.v get empirical density
  empeigen <- density(eigen.v, from = min(eigen.v), to = max(eigen.v), cut = 0);
  dim <- length(which(eigen.v > theigen$max))

  if(plot){
    theigen.df <- data.frame(lambda = theigen$lambda, density = theigen$dens);
    i <- min(which(empeigen$x >= min(eigen.v)))
    f <- max(which(empeigen$x <= max(eigen.v)))
    empeigen.df <- data.frame(lambda = empeigen$x[i:f], density = empeigen$y[i:f])
    obj <- ggplot() + xlim(min(theigen.df$lambda, empeigen.df$lambda),max(theigen.df$lambda, empeigen.df$lambda)) + ylim(0,max(theigen.df$density, empeigen.df$density))
    # draw random matrix
    obj <- obj + geom_line(data = theigen.df, aes(x = lambda, y = density), colour = "darkgreen", size = 1.2)
    # draw C matrix
    obj <- obj + geom_line(data = empeigen.df, aes(x = lambda, y = density), colour = "red", size = 1.2)
    # abline marker $K$ lambdas which are larger than lambdaMAX
    obj <- obj + geom_vline(xintercept = eigen.v[1:dim], alpha = 0.75, colour = "darkblue", lty = "dashed", size = 0.2)
    return(list(cor = C, dim = dim, empeigen = empeigen, theigen = theigen, evalue = eigen.v, ggplot = obj))
  }else{
    return(list(cor = C, dim = dim, empeigen = empeigen, theigen = theigen, evalue = eigen.v))
  }
}

#' isvaFunction
#'
#' @title isvaFn
#' @param dat.m Data matrix whose rows represent different labeling features and columns stand for the different samples
#' @param design design matrix for expression/methylation microarray
#' @param qcutoff qvalue's cutoff, control the FDR(false discovery rate)
#' @return n.isv number of ISVs
#' @return isv matrix of ISV
#' @export
#' @importFrom qvalue qvalue
#' @importFrom fastICA fastICA
#' @author Xin Zhou \url{xinchoubiology@@gmail.com}
isvaFn <- function(dat.m = NULL, design = NULL, type = c("M", "Beta"), qcutoff = 0.05, verbose = FALSE){
  type <- match.arg(type)
  if(type == "Beta"){
    cat(sprintf("M values transform ...\n"))
    dat.m <- log2(dat.m / (1 - dat.m))
  }
  lm.o <- lm(t(dat.m) ~ -1 + design);
  res.m <- t(lm.o$residuals);
  ## estimate dimension for ICA
  rmt.o <- estDimRMT(res.m, plot = FALSE);
  ncomp <- rmt.o$dim;
  cat(sprintf("rmt evaluated ISVs dimension is %d \n", ncomp));

  if(ncomp == 0){
    return(list(n.isv = ncomp, isv = NULL))
  }

  ## perform ICA on res.m
  ICA.o <- fastICA(res.m,n.comp = ncomp);

  tmp.m <- t(ICA.o$A)
  isv.m <- tmp.m
  sd    <- 1/sqrt(ncol(dat.m) - 3)
  for(i in 1:ncol(tmp.m)){
    cor.v <- as.vector(cor(t(dat.m),tmp.m[, i]))                   ## for each component calculate its correlation with each CpG's spectrum
    z.v <- 0.5*log((1+cor.v)/(1-cor.v))
    pv.v <- 2 * pnorm(abs(z.v), 0, sd, lower.tail = FALSE)         ## calculate P-value for each regression of CpG
    tmp.s <- sort(pv.v, decreasing = FALSE, index.return = TRUE)   ## return list(x = ...,ix = ...)
    qv.o <- qvalue(p = pv.v)
    nsig <- length(which(qv.o$qvalue <= qcutoff))
    ## select best feature
    if(nsig < 500){
      nsig <- 500
    }
    ## reductive fastICA
    reduct.m <- dat.m[tmp.s$ix[1:nsig],]
    ICA.o <- fastICA(reduct.m, n.comp = ncomp)
    ## select best correlated vectors A*
    cor.v <- abs(cor(tmp.m[,i], t(ICA.o$A)))
    kmax <- which.max(cor.v)
    isv.m[,i] <- t(ICA.o$A)[,kmax]
    if(verbose){
      cat(sprintf("Build %d SV ...\n", i))
    }
  }
  return(list(n.sv = ncomp, sv = isv.m))
}

#' svaReg Do regression with selected surrogate variables
#'
#' @title svaReg
#' @description when we search all of the factors which means surrogate variables
#'              update our model.matrix and calculate moderate test
#' @param dat.m n x m M|beta matrix for n CpG sites across m different patient samples.
#' @param design phenotype of interested; m-length vector reresent patient samples' phenotype.
#'        Optional: [combine with sample pair information] and mandatory column 2 is covariate
#'        of interest
#' @param sv.m  surrogate variables matrix calculate from \link{isvaFn} other sva methods.
#' @param qvalue0 false discovery rate's threshold; Default = 0.1
#' @param backend backend regression packages; Default = "NULL", switch to limma for moderate statistic
#' @param verbose Optional; Default FALSE
#' @importFrom limma lmFit eBayes contrasts.fit
#' @importFrom qvalue qvalue
#' @return res data.frame
#' @export
#' @author Xin Zhou \url{xinchoubiology@@gmail.com}
svaReg <- function(dat.m = NULL, design = NULL, sv.m = NULL, qvalue0 = 0.1, backend = c("NULL", "limma"), verbose = FALSE){
  backend <- match.arg(backend)
  if(is.null(sv.m)){
    modelSv <- model.matrix(~ -1 + design)
    modelNull <- model.matrix(~ -1 + design[,-2])
  }
  else{
    modelSv <- model.matrix(~ -1 + design + sv.m)
    modelNull <- model.matrix(~ -1 + sv.m + design[,-2])
  }

  if(backend == "NULL"){
    pval   <- pvalue(dat.m, modelSv, modelNull)
    pval.s <- sort(pval, decreasing = FALSE, index.return = TRUE)
    qval   <- qvalue(pval.s$x)$qvalue
    nsig   <- length(which(qval < qvalue0))
    df     <- ncol(dat.m) - ncol(modelSv)
    cat(sprintf("Numeber DMP = %d \n", nsig))
    if(nsig > 0){
      ## p.index <- pval.s$ix[1:nsig]
      p.index <- pval.s$ix
      tmp.lm  <- lm(t(dat.m[p.index,]) ~ modelSv - 1)
      lm.stat <- do.call(rbind, lapply(summary(tmp.lm), function(x) x$coefficients[2,c(1,3)]))
      ## res     <- cbind(lm.stat, pval[p.index], qval[1:nsig])
      res     <- cbind(lm.stat, pval[p.index], qval)
      colnames(res) <- c("Estimate", "t-stats", "f-pvalue", "fdr")
      rownames(res) <- rownames(dat.m[p.index,])
      res <- list(siggene = res[1:nsig,], null = res[(nsig+1):nrow(dat.m),])
    }
    else{
      res <- NULL
    }
    return(res)
  }
  else if(backend == "limma"){
      fit  <- lmFit(dat.m, modelSv)
      nsv  <- ncol(modelSv) - 2
      contrast.matrix <- cbind("cancer-normal" = c(0, 1, rep(0, nsv)))
      fitc <- contrasts.fit(fit,contrast.matrix)
      fitc <- eBayes(fitc)
      ## res  <- topTableF(fitc, number = nsig)
      res  <- topTableF(fitc, number = Inf)
      qval <- qvalue(res$P.Value)$qvalue
      nsig <- length(which(qval < qvalue0))
      cat(sprintf("Numeber DMP = %d \n", nsig))
      res  <- list(siggene = res[1:nsig,], null = res[(nsig+1):nrow(dat.m),])
      return(res)
  }
}
