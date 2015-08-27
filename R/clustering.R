# functions for clustering regions; regions has probe(s)[predefined clusters]
# clustering methods contains three candidate methods:
# 1) Hierarchical Clustering(HC) & DynamicTreeCut (WGCNA)
# 2) NMF method (robust)[But NMF must be define the k]
# Using Package WGCNA

#' get eigen-gene from pre-clustered prototype [Methods: BMC2007_Steve_Horvath]
#' and the second-layer network analysis will be based on [eigengene network]
#' @title Eigenbeta
#' @description CpG sites are clustered by bump-hunter or combp. So the correlation(similarity) network
#'              will be represented by eigen unit(which is the 1st PC)
#' @param mset    normalized methylation matrix of 485512 CpG sites
#' @param cluster List; clustered neighboured CpG sites 
#' @param nPC     number of principle components
#' @param verbose Optional; TRUE(default), whether to print processing or not
#' @param softpower 6; For hub beta vector generation
#' @param corr    Correlation method; Optional; "pearson" default; "spearman" and "kendall"
#' @param align   Use the average beta(Methylation) value to control the orientation of eigenunit
#'                align =  "along average"(default, means that eigenunit should be aligned with average beta value)
#'                align =  "" (means letting the eigenunit vector as themselves)
#' @return list
#'         eigenBeta eigen beta vectors' matrix for pre-clustered prototype units
#'         averageBeta average beta vectors for each clustered CpG unit
#' @export
#' @author Xin Zhou \url{xxz220@@miami.edu}
Eigenbeta <- function(mset = NULL, cluster = NULL, nPC = 1, verbose = TRUE,
                      softpower = 6, corr = c("pearson", "spearman", "kendall"),
                      align = c("along average", "")){
  options(warn = -1)
  
  align <- match.arg(align)  
  corr  <- match.arg(corr)
  if(is.null(mset) || is.null(cluster)){
    stop("Methylation beta matrix & cluster definition list are needed by function Eigenbeta ...")
  }
  if(!is.element(align, c("along average", ""))){
    stop("Eigen vector must be oriented by [along average] method or [None] method ... ")
  }
  
  ## TRUE for region should be calculate their Eigenbeta vector
  combRegion  <- (cluster[[4]] != 1)
  eigenBeta   <- matrix(NA, nrow = length(combRegion), ncol = ncol(mset))
  averageBeta <- matrix(NA, nrow = length(combRegion), ncol = ncol(mset))
  
  eigenBeta[!combRegion,]   <- mset[cluster[!combRegion][[5]],]
  averageBeta[!combRegion,] <- mset[cluster[!combRegion][[5]],]
  
  rownames(eigenBeta)   <- 1:length(combRegion)
  colnames(eigenBeta)   <- colnames(mset)
  rownames(averageBeta) <- 1:length(combRegion)
  colnames(averageBeta) <- colnames(mset)
  if(verbose){
    print(paste("  moduleEigen vectors : Calculating ", sum(combRegion), "combined CpG regions' eigen-beta vector."))
  }
  
  for(i in which(combRegion)){
    index   <- str_split(cluster[i][['probes']], pattern = ";")[[1]]
    dat.mat <- mset[index, ,drop = FALSE]
    n       <- nrow(dat.mat)
    p       <- ncol(dat.mat)
    
    if(verbose){
      print(" Calculating SVD ...")
    }
    usv <- svd(dat.mat, nu = min(n, p, nPC), nv = min(n, p, nPC))
    pc  <- usv$v[,1]
    
    ## assignment
    eigenBeta[i,] <- pc
    if(align == "along average"){
      if(verbose){
        print(" aligning eigenBeta with average beta value")
      }
      averageBeta[i,] <- colMeans(dat.mat)
      if(cor(averageBeta[i,], eigenBeta[i,], use = "p") < 0){
        eigenBeta[i,] <- eigenBeta[i,] * -1
      }
    }
  }
  
  list(eigenBeta = eigenBeta, averageBeta = averageBeta)
}


