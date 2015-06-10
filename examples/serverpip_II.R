load(".RData")
library(Rcppsva)
library(doMC)

system.time(combpregion <<- combine.pvalue(dat.m = mset.ComBat, pvalues = p, cluster = corrcluster, chr = BED$chr, pos = BED$pos, names = rownames(BED), method = "spearman", combine = "stouffer_liptak"))

save.image(file = "n2.RData")
