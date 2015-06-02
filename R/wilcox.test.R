##' If our data is paired, wilcoxon rank test can be used to test the significnat of sites
##' 
##' @title wilcox.fit
##' @description calculate wolcoxon signed rank sum test for each CpG probes of microarray;
##'              support multiple core parallism calculation.
##' @param dat.m n x m delta M|beta matrix for n CpG sites across 2*m paired different patient samples
##' @param alternative a character string specifying the alternative hypothesis; 
##'        c("two.sided", "greater", "less") and "two.sided" is default
##' @param qvalue0 false discovery rate's threshold; Default = 0.1
##' @param mcore Cpu cores can be used in test
##' @return res list
##'             siggene data.frame
##'             null    data.frame
##' @importFrom doMC registerDoMC
##' @importFrom plyr aaply
##' @importFrom qvalue qvalue
##' @export
##' @author Xin Zhou
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