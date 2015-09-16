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

#' based on WGCNA and Dynamic Tree Cut, clustering the DMRs
#' Dynamic Tree cut method to determine co-methylation modules of different regions
#' @title moduleSearch
#' @description cutreeDynamic is a dynamic clustering method borrowed from package dynamicTreeCut,
#'              and in this method, we will clustered our pre-clustered beta matrix and calculate
#'              out their co-regulation methylation regions
#' @param beta.expr pre-clustered genome region for co-regulation methylation status
#' @param cor.type  correlation type for dissimilarity calculation; "pearson"(Default) and "spearman"
#' @param sim.type  similarity type for hierachical clustering; "signed"(Default) and "unsigned"
#' @param method    hierachical clustering method; c("average", "ward", "single", "complete")
#' @param deepSplit [0,1,2,3,4]. Parameter in cutreeDynamic to control the sensitivity of module detection
#' @param cutHeight h_{max} in DynamicTreeCut method
#' @param minClusterSize c_{min} minimum size of each cluster
#' @param pamStage  Only used for method "hybrid" dynamicTreecut; If TRUE, the (PAM-like)stage will be
#'                  performed
#' @param pamRespectsDendro Logical
#' @param usebigmemory use bigmemory or not 
#' @return A list with following component:
#'         clusters  a vector of cluster labels for all CpG regions
#'         dendrogram dendrogram of Rclusterpp.hclust
#' @export
#' @importFrom dynamicTreeCut cutreeDynamic
#' @author Xin Zhou \url{xxz220@@miami.edu}
#  Notice！This function has somthing wrong on memory management
#  I'll introduce the Rclusterpp module in my package and finally
#  I can directly transfer the beta value matrix in to hclust function
moduleSearch <- function(beta.expr = NULL, cor.type = c("pearson", "spearman"),
                         sim.type = c("S", "U"),
                         method   = c("average", "ward", "single", "complete"),
                         deepSplit = 1, cutHeight = 0.9, 
                         minClusterSize = min(10, ncol(beta.expr)/2), 
                         pamStage = TRUE, pamRespectsDendro = TRUE, 
                         verbose = 2, ...){
  options(warn = -1)
  if(is.null(beta.expr)){
    stop("methylation beta value matrix is needed by moduleSearch function ...")
  }
  cor.type <- match.arg(cor.type)
  sim.type <- match.arg(sim.type)
  method   <- match.arg(method)
  
  
  # clustering by Rclusterpp.hclust
  dendro <- HClust(beta.expr, method = method, distance = cor.type, sign = sim.type)
  
  # reference dendrogram
  dendro.ref <- distributeRef(data = beta.expr, size = 20, by = 'R', n = nrow(beta.expr), distance = "pearson", sign = "S")
  
  dendromodule <- clusterDetect(data = beta.expr, dendro = dendro, dendref = dendro.ref, log =TRUE, distance = cor.type, sign = sim.type)
  
  modules   <- optimModule(dendromodule, dendro = dendro, plot = TRUE)
  
  modules
}

