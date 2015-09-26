## Package Rcppsva : Package of comprehensive differential methylation region motif search engine

package reconstruction path:
generic definition : AllGeneric.R 
class definition   : AllClasses.R
* Flowchart:
  + Raw data reading : `minfi package`
    - minfi::preprocessRaw
    - minfi::preprocessIllumina
    - minfi::getBeta minfi::getM
    - minfi::preprocessSWAN
    
  + Batch effect elimination : `combat.R`
    - ComBat
  
  + Surrogate variable analysis : `sva.R`
    - estDimRMT
    - isva
  
  + Normalized beta/M value expression matrix ~ design matrix(phenotype + surrogate variables) : `sva.R`
    - svaReg
  
  + Sampling null distribution for coefficient of phenotype : `sampling.R`
    - mlm.fit : permutation sampling
    - bootstrap.fit : bootstrap sampling
  
  + Clustering CpG sites by distance : `misc.R`( Given spatial constraint on genome )
    - agglomerateBydist
    - diffbumpFilter (filter Bumps by difference p-value of null distribution)
    - corrbumpFilter (filter Bumps by correlation constraint of null distribution)
  
  + Dbpmerge -- Agglomerative clustering for correlation CpG regions search : `merge.R`
    - Dbpmerge (Pre-clustered region are tested by `combinePval` method)
  
  + Combin p-value method via stouffer-liptik method : `combp.R`
    - combinep (combine pvalue method)
      - stoufferliptakCombp
      - zscoreCombp
  
  + Local weight smooth method : `lowess.R`
    - lowess 
  
  + Hierachical clustering : `hclust.R`
    - eigenM : representing M of DMR(regions) by eigen-value & sample-wise mean
    - hclust
    - referenceDist
    - gapStat
    - optimalHC
    - ultilities: `util.R`:
      - dendrogram dispatch methods
  
  + Auxiliary function : `dendrogram.R`
    - dendrogram plot [not necessary]
  
  + Unirform interface : `bumphunter.R`
    - bumphunterEigen
    
  + TCGA-Assembler Downloading : `get450K.R`
    - downloadMethylURL
    - getClinical
    - getMethyltion
    - get450K (read in methylation 450K data and coefficients)

  + BMarry Class abstraction
    - methods wrapper for all above methods
  
  + 450K preprocessing 
    - getTCGAssemble (load TCGA Assembler Package in /inst)
    - get450Kurl 
    - get450Ksheet
    - get450Kzip
