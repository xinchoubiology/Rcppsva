## smoother functions for beta value have uniform callable 
## smoother(beta, chr, pos, cluster, names, weight, verbose = TRUE)

#' robust local weighted regression, fitting a smooth nonparamatic regression
#' and smoothing methylation measurement.
#' 
#' @title smoother
#' @param beta beta value for probes
#' @param pos position vector
#' @param names probe names vector ; used in function \link{clusterMaker} and \link{regionSeeker}
#' @param cluster cluster object for local regression
#' @param weight weight for local regression
#' @param method c("weightedLowess", "loess", "locfit")
#' @param mcores thread used in function
#' @importFrom limma loessFit
#' @importFrom parallel detectCores
#' @importFrom doMC registerDoMC
#' @export
smoother <- function(beta, pos, names, cluster, weight, 
                     method = c("weightedLowess", "loess", "locfit"),
                     mcores = detectCores()){
  method = match.arg(method)
  # List2NumericVector
  if(inherits(cluster, "numeric")){
    cluster <- split(names(cluster), cluster)
  }
  registerDoMC(cores = mcores)
  fit <- llply(cluster, .fun = function(c){
                                if(length(c) >= 2){
                                  y <- beta[c,]
                                  x <- pos[which(names %in% c)]
                                  out <- loessFit(y, x, weights = weight[c,], method = method)
                                  return(out$fitted)
                                }else{
                                  return(beta[c,])
                                }
                               }, .parallel = TRUE)
  
  beta <- matrix(unlist(fit), ncol = 1)
  rownames(beta) <- unlist(cluster)
  beta
}