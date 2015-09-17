#' grid search for CLuster number detection by gap statistic
#' @title gapStat
#' @description Based on dissmilarity measurement, we define Wk = sum(sum(d_{ij})/2n_r)
#' @param data    methylation data matrix
#' @param dendro  actual dendrogram
#' @param dendref reference dendrogram
#' @param cmax    maximum dissimilarity Default 0.22
#' @param cmin    minimun dissimilarity Default 0.02
#' @param by      step size of height serarch
#' @param log     FALSE(Default). Standardize by log or not
#' @param mcore   number of threads running
#' @param verbose verbose 
#' @param ...     optional parameters for \code{distributeRef} and \code{HClust}
#' @return list of clusterDetect
#'         labels
#'         number
#'         Wk  within sum of dissimilarity
#'         EWk expectation of Wk under null reference distribution
#' @export
#' @author Xin Zhou \url{xxz220@@miami.edu}
gapStat <- function(data = NULL, dendro = NULL, dendref = NULL, 
                    cmax = 0.24, cmin = 0.04, by = -0.005, mcore = 2,
                    log = FALSE, verbose = TRUE, ...){
  options(warn = -1)
  if(is.null(data)){
    stop("Methylation beta / M matrix is needed for clustering")
  }
  params <- list(...)
  print(params)
  if(is.null(dendro)){
    dendro <- HClust(data = data, distance = params$distance, sign = params$sign)
  }
  if(is.null(dendref)){
    dendref <- distributeRef(data = data, ...)
  }
  
  if(mcore >= 2){
    registerDoMC(cores = mcore)
  }
  
  K <- NULL
  Wk       <- NULL
  ExpWk    <- NULL
  SWk      <- NULL
  Ik       <- NULL   # Inter clusters dissimilarity
  Hk       <- NULL
  
  # cut clustering by number
  kmax <- sum(dendro$height >= cmax) + 1
  kmin <- sum(dendro$height >= cmin) + 1
  for(knum in kmax:kmin){
    if(verbose){
      cat("  Detecting clustering performance of height cut @ ", height, "...\n")
    }
    Htree  <- cutree(dendro, h = knum)
    Hgroup <- split(1:length(Htree), Htree)
    Hcorr  <- llply(Hgroup, function(x){
                              if(length(x) <= 1){
                                0
                              } else{
                                N <- length(x)
                                W <- cor(t(data[x,]), method = params$distance)
                                if(params$sign == "S"){
                                  sum((1 - W) / 2) / (2 * N)
                                } else{
                                  sum((1 - abs(W))) / (2 * N)
                                }
                              }
                            }, .parallel = TRUE)
    
    ## get cluster number
    K      <- c(K, knum)
    height <- tail(dendro$height, knum-1)[1]
    Hk     <- c(Hk, height)
    
    Wk0    <- NULL
    # specify the data imbalance scale
    # scale <- nrow(data) / nrow(dendref[[1]]$data)
    for(i in 1:length(dendref)){
      Reftree  <- cutree(dendref[[i]]$hclust, k = knum)
      Refgroup <- split(1:length(Reftree), Reftree)
      Refcor   <- llply(Refgroup, function(x){
                                    if(length(x) <= 1){
                                      0
                                    } else{
                                      N <- length(x)
                                      W <- cor(t(dendref[[i]]$data[x,]), method = params$distance)
                                      if(params$sign == "S"){
                                        sum((1 - W) / 2) / (2 * N) # * scale
                                      } else{
                                        sum(1 - abs(W)) / (2 * N)  # * scale
                                      }
                                    }
                                  }, .parallel = TRUE)
      Wk0 <- c(Wk0, do.call(sum, Refcor))
    }
    if(log){
      ExpWk <- c(ExpWk, mean(log(Wk0)))
      SWk   <- sqrt(1+1/length(dendref)) * sd(log(Wk0))
      Wk    <- c(Wk, log(do.call(sum, Hcorr)))
    } else{
      ExpWk <- c(ExpWk, mean(Wk0))
      SWk   <- sqrt(1+1/length(dendref)) * sd(Wk0)
      Wk    <- c(Wk, do.call(sum, Hcorr))
    }
  }
  ret <- data.frame(H = Hk, K = K, Wk = Wk, EWk = ExpWk, Gap = ExpWk - Wk, Sk = SWk)
  
  ret
}

#' optimal cluster number retrieve
#' @title optimModule
#' @param gstat gap statistic data
#' @param dendro dendrogram object of actual clustering
#' @param plot  Optional, plot or not
#' @return module object
#' @export
optimModule <- function(gstat = NULL, dendro = NULL, plot = TRUE, verbose = TRUE){
  options(warn = -1)
  if(is.null(gstat) || is.null(dendro)){
    stop("gap statistic data frame and dendrogram object are needed for concersion ... ")
  }
  cat("Gap statistic for cluster number ... \n")
  ggWk  <- ggplot(gstat, aes(x = K, y = Wk)) + 
    geom_line(lty = 1, lwd = 1.5) + ggtitle("number & cluster") + 
    geom_line(aes(x = K, y = EWk), lty = 1, lwd = 1.5, color = "blue") +
    geom_point()
  
  ggGap <- ggplot(gstat, aes(x = K, y = Gap)) + geom_line(lty = 1, lwd = 1.5) + 
    geom_point() +  geom_errorbar(aes(ymin=Gap-Sk, ymax=Gap+Sk), width = 2, color = "blue") + ggtitle("number & cluster gap statistic")
  
  if(plot){
    ggGap
  }
  
  ixmax <- which.max(gstat$Gap)
  Kmin  <- NULL
  for(i in 1:(nrow(gstat)-1)){
    if(gstat$Gap[i] >= gstat$Gap[i+1] + gstat$Sk[i+1]){
      Kmin <- gstat$K[i]
      break
    }
  }
  trees <- cutree(dendro, h = gstat$H[ixmax])
  ret   <- list(labels = trees, Kmax = gstat$K[ixmax], Kmin = Kmin, height = gstat$H[ixmax], ggWk = ggWk, ggGap = ggGap)
  class(ret) <- "module"
  
  ret
}



#' @rdname print
#' @param x module object
#' @param ... other arguments
#' @return NULL
#' @method print module
#' @export
print.module <- function(x, ...){
  cat(sprintf(" %d clusters based on hierachical clustering gap statistic \n", x$Kmax))
  cat(sprintf(" minimal cluster number is %d", x$Kmin))
  cat(sprintf(" cut height threshold = %f \n", x$height))
  cat(sprintf(" minimal correlation within defined clusters = %f \n", 1 - 2 * x$height))
}