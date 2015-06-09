load(".RData")
library(Rcppsva)
library(doMC)

##########################################
# Test correlated matrix list calculation
##########################################
rawCluster <- split(names(cluster), cluster)
combIndex <- which(sapply(rawCluster, function(c) length(c) > 1))
multClust <- rawCluster[combIndex]
pos <- BED[,2]
names(pos) <- rownames(BED)

registerDoMC(cores = detectCores())

start <- Sys.time()
corrmatrix  <- llply(multClust, .fun = function(ix){
  dist <- abs(outer(pos[ix], pos[ix], "-"))
  colnames(dist) <- rownames(dist) <- ix
  cor(t(mset.ComBat[ix,]), method = "spearman") * (dist <= 1000)
}, .parallel = TRUE)

end  <- Sys.time()

cat(sprintf("Time elasped = %f mins On %d cores \n", end - start, detectCores()))

## cluster by Dbpmerge
## test on 16 cores ~ 
corregions <- Dbpmerge(corrmatrix, merge = "average", cutoff = 0.8)

save.image(file = ".RData")