## pipeline dependence package
library(minfi)
library(minfiData)
library(sva)
library(Rcppsva)
library(stringr)

## minfiData has direction "extdata"
##                           |------ 2 Sentrix ID's Folder
##                           |------ SampleSheet.csv (same peopel has same ID; id1, id2, id3)
## Lvl.1 Data store in baseDir
baseDir <- system.file("extdata", package = "minfiData")

## Read 450K
list.files(baseDir)  ## show different slides
## IDAT file
list.files(file.path(baseDir, "5723646052"))  ## show different arrays

## targets
targets <- read.450k.sheet(baseDir)  ## different sentrix ID means different Slides(Batch)
## different sentrix position means different Array(not Batch)

## read microArray data
RGset   <- read.450k.exp(targets = targets)

# phenotype pd; in phenotye table, 
# the Sample_Groupe is treated as factor of status[normal; cancer]
pd <- pData(RGset)

## QC quality
qcReport(rgSet = RGset, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Group, pdf = "minif_test.pdf")

## Convert from RGset to Methyl and Unmethyl
mset.Raw <- preprocessRaw(RGset)
## Control Probes = 622399 - 485512 = 136887; 450K probes = [Type-I + Type-II]

## Control Probes is useful, add the control process before convert
mset.Norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 1)

## Via the mset Norm it is probability that U & M are all under the Background Signal
## As the result, we use M/(M + U + 100) for beta, and logit(beta) for M
## How to choose beta threshold is important
## Also define a CN value(copy number)
Beta <-  getBeta(mset.Norm,type = "Illumina")
## Zeros in your dataset have been replaced with 0.000001
M    <-  getM(mset.Norm, type = "Beta", offset = 100, betaThreshold = 0.000001)

## MDS plot ==> First sense of relationship between samples
## Multiple Dimension scaling: Visualize the similarity of individual cases of a dataset.
mdsPlot(M, numPositions = 10000, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
## cmdscale is to obtain firsr sense the different between different sample
## They are batch effect

## Another question is that Type I probe and Type II probe are different
## Use SWAN to erase the technique different between 2 type
mset.SWAN <- preprocessSWAN(RGset, mset.Norm)
par(mfrow = c(1,2))
plotBetasByType(mset.Norm[,1], main = "Raw")
plotBetasByType(mset.SWAN[,1], main = "SWAN")

## ComBat for batch effect
## In ChAMP the code told us: batch = pd$Slide.
## Therefore in extdata, Slide = 5723646052 & 5723646053
## Erase Batch effect by ComBat (dat - batch_effect)
M <- getM(mset.SWAN, type = "Beta", offset = 100, betaThreshold = 0.000001)

## + person : Taking person pair information into account
pheno.v  <- factor(pd$status, levels = c("normal", "cancer"))
mod   <- model.matrix(~ pheno.v + person, data = pd)
mod0  <- model.matrix(~ person, data = pd)

mset.ComBat <- ComBat(dat = M, batch = pd$Slide, mod = mod, par.prior = FALSE)

## Before next step, I did the profiler for ComBat function
## library(lineprof)
## tmp <- lineprof(ComBat(dat = M[1:10000,], batch = pd$Slide, mod = mod, par.prior = FALSE), torture = TRUE)
## shine(tmp)

## choose Combat-n cuz Paper: ChaoChen_PlosOne_2011
## after ComBat, we use SVA to learn all the other unmodel factors
## mset.SV <- isvaFn(dat.m = mset.ComBat, design = mod, type = "M", verbose = TRUE)
n.sv <- num.sv(dat = mset.ComBat, mod = mod, method = "leek")
mset.SV <- sva(dat = mset.ComBat, mod = mod, mod0 = mod0, n.sv = n.sv)

## calculate DMP(differential methylation position)
## pheno.v  <- factor(pd$status, levels = c("normal", "cancer"))
## rank "cancer" followed by "normal" so estimated means cancer - normal

## test differential methylation(siggene) by linear regression on M data(link function from beta is logit)
## generally, siggene found by "NULL" method will be less that those found by "limma" eBayes method
mset.DMP <- svaReg(dat.m = mset.ComBat, design = mod, sv.m = mset.SV$sv, qvalue0 = 0.05, backend = "NULL", verbose = TRUE)
mset.DMP.limma <- svaReg(dat.m = mset.ComBat, design = mod, sv.m = mset.SV$sv, qvalue0 = 0.1, backend = "limma", verbose = TRUE)

## plot CpG
cpgs <- rownames(mset.DMP.limma$null)[85680:85683]
par(mfrow = c(2,2))
plotCpg(mset.ComBat, cpg = cpgs, pheno = pd$Sample_Group)

## 3 methods for differential methylation regions search 
## whether paired-limma could improve predictiion
sim.data <- function(){
  sd   <- 0.3 * sqrt(4/rchisq(100, df = 4))
  y  <- matrix(rnorm(100*6, sd = sd), 100, 6)
  rownames(y) <- paste("cg", 1:100)
  for(i in 1:3){
    y[,(i+3)] = y[,i] + rnorm(100, sd =  abs(mean(y[,i])))
    y[1:3,(i+3)] = y[1:3,i] + c(0.1,1,2);
  }
  targets  <- data.frame(pair = factor(c(1,2,3,1,2,3)), type = factor(c(0,0,0,1,1,1)))
  design   <- model.matrix(~pair+type, data = targets)
  design_n <- model.matrix(~type, data = targets)
  fit <- lmFit(y, design)
  fit <- eBayes(fit)
  fit_n <- lmFit(y, design_n)
  fit_n <- eBayes(fit_n)
  res <- topTable(fit, coef = 4)
  res_n <- topTable(fit_n, coef = 2)
}

## if we wish to use the wilcoxon rank sum test
## if pdata has person pairs and status information, we can use split to 
## Wilcoxon fit has bugs
paired.v <- split(pd, pd$person)
paired.v <- lapply(paired.v, function(x) return(x[order(x$status),]))
mset.Delta <- mset.ComBat[,sapply(paired.v, function(x) rownames(x[which(x$status == "cancer"),]))] - mset.ComBat[,sapply(paired.v, function(x) rownames(x[which(x$status == "normal"),]))]
mset.DMP.wilcox <- wilcox.fit(dat.m = mset.Delta, alternative = "two.sided", qvalue0 = 0.1)

## create modified p-value data frame
mset.limma <- list(table = rbind(mset.DMP.limma$siggene, mset.DMP.limma$null))
class(mset.limma) <- "bed"
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
print(dat.m = mset.limma, bed = "siggene.bed", db = IlluminaHumanMethylation450kanno.ilmn12.hg19)

## comb-p pipeline --seed 0.01 --dist 1000 --acf-dist 1:1000:50 --step 1000 -p out.prefix siggene.bed
modsv <- cbind(mod, sv = mset.SV$sv)
tmp <- mlm.fit(dat.m = mset.ComBat[1:10,], design = modsv, coef = 2, B = NULL, full = TRUE)
tmp <- bootstrap.fit(dat.m = mset.ComBat[1:10,], design = modsv, coef = 2, B = 100)

## test bump hunting kernel
Location <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
Probes <- rownames(Location)[rownames(Location) %in% rownames(dat.m)]
BED <- cbind(Location[Probes, 1:2])
BED <- cbind(BED, dat.m[Probes, 4:3])
## elapsed ~ 1.345 secs
cluster <- clusterMaker(chr = BED$chr, pos = BED$pos, maxGap = 1000, names = rownames(BED))

tmp <- mlm.fit(dat.m = mset.ComBat, design = modsv, coef = 2, B = NULL, full = TRUE)
## Segment extract from predefined regions
## Test segmentMaker
B <- tmp$coef
cutoff <- c(-quantile(abs(B), 0.8), -quantile(abs(B), 0.2), quantile(abs(B), 0.2), quantile(B, 0.8))

## CpG information
genome     <- GRanges(seqnames = BED$chr, ranges = IRanges::IRanges(start = BED$pos, width = 1), 
                      names = rownames(BED), chr = BED$chr)

seqlevels(genome) <- sort(seqlevels(genome))
genome     <- sort(genome)

## elapsed ~ 1.160 secs
segment.set <- segmentsMaker(cluster, B, cutoff, genome$chr, start(genome))

## regionSeeker
## elapsed of regionSeeker : 96.543 secs
B <- tmp$coef
cutoff = c(quantile(abs(B), 0.8), quantile(abs(B), 0.2))
regions <- regionSeeker(beta = B, chr = BED$chr, pos = BED$pos, maxGap = 1000, names = rownames(BED), drop = TRUE)

## permutation
## elapsed = 20.058 secs
## tmp <- bootstrap.fit(dat.m = mset.ComBat, design = modsv, coef = 2, B = 100)
## tmp1 <- mlm.fit(dat.m = mset.ComBat[1:10,], design = modsv, coef = 2, B = 10, full = TRUE)
## tmp1$coef <- t(apply(tmp1$coef, 1, '/', t(tmp1$stdev_unscaled)))
## tmp2 <- bootstrap.fit(dat.m = mset.ComBat[1:10,], design = modsv, coef = 2, B = 10)
## tmp <- mlm.fit(dat.m = mset.ComBat, design = modsv, coef = 2, B = 100, full = TRUE)
## tmp$coef <- t(apply(tmp$coef, 1, '/', t(tmp$stdev_unscaled)))
tmp <- bootstrap.fit(dat.m = mset.ComBat, design = modsv, coef = 2, B = 100)
permB   <- tmp$coef

## segments extraction processing
cutoff <- c(-quantile(abs(B), 0.8), -quantile(abs(B), 0.2), quantile(abs(B), 0.2), quantile(B, 0.8))
segment_perm <- segmentsMaker(cluster, B, cutoff, genome$chr, start(genome), permB)

regions_perm <- regionSeeker(beta = B, chr = BED$chr, pos = BED$pos, maxGap = 1000, permbeta = permB, names = rownames(BED), drop = TRUE)

##########################################################
## correlated regions extracted && test Dbpmerge function
##########################################################
rawCluster <- split(names(cluster), cluster)
combIndex <- which(sapply(rawCluster, function(c) length(c) > 1))
multClust <- rawCluster[combIndex]
pos <- BED[,2]
names(pos) <- rownames(BED)

## corrmatrix calculation too slow
## working on 16 cores ~ 450 secs
## Time difference of 8.203694 mins on 16 cores
registerDoMC(cores = detectCores())

start <- Sys.time()
## -----
corrmatrix  <- llply(multClust, .fun = function(ix){
  dist <- abs(outer(pos[ix], pos[ix], "-"))
  cor(t(mset.ComBat[ix,]), method = "spearman") * (dist <= 1000)
}, .parallel = TRUE)
## -----
end <- Sys.time()


## cluster by Dbpmerge
## test on 4 cores ~ 55 secs for 55003 clusters  => 287288 clusters
## But cluster contains more than one probes => 45203 clusters
corregions <- Dbpmerge(corrmatrix, merge = "average", cutoff = 0.8)
names(corregions) <- seq(1, length(corregions))
corrClust  <- corregions[which(sapply(corregions, function(c) length(c) > 1))]

## Visualization correlation
corrHeatmap <- llply(corrClust[1:100], .fun = function(ix){
  dist <- abs(outer(pos[ix], pos[ix], "-"))
  cor(t(mset.ComBat[ix,]), method = "spearman") * (dist <= 1000)
}, .parallel = TRUE)

## correlation visualization
## cormatrix <- list(cor = corrHeatmap$`103.13`)
## class(cormatrix) <- "corr"
## plot(cormatrix, cutoff = 0.80)

## build correlated clusters
## for 1/4 probes ~ 403 secs => all probes 1600 secs
## corrcluster <- corrclusterMaker(dat.m = mset.ComBat, chr = BED$chr, pos = BED$pos, names = rownames(BED), maxGap = 1000, cutoff = 0.8, merge = "average", method = "spearman")
corrcluster <- c(corregions, rawCluster[-combIndex])
names(corrcluster) <- seq(1, length(corrcluster))

## robust estimation and p-value calculation
tmpr <- mlm.fit(dat.m = mset.ComBat, design = modsv, coef = 2, B = NULL, full = TRUE)
betar <- tmpr$coef
sigmar <- tmpr$sigma
fitr <- mlm.tstat(tmpr)
tr   <- betar / fitr$post.sigma
p   <- 2 * pmin(pt(tr, fitr$df.total), 1 - pt(tr, fitr$df.total))

## get multiple-probe regions
multiregions <- corrcluster[which(sapply(corrcluster, length) >= 2)]

## normal cluster
regions_perm <- regionSeeker(beta = B, chr = BED$chr, pos = BED$pos, cluster = cluster, maxGap = 1000, permbeta = permB, names = rownames(BED), drop = TRUE)

## permcorregion
corrcluster <- corrclusterMaker(dat.m = mset.ComBat, chr = BED$chr, pos = BED$pos, names = rownames(BED), corrmat = corrmatrix, maxGap = 1000, cutoff = 0.8, merge = "average", method = "spearman")
corregion_perm <- regionSeeker(beta = B, chr = BED$chr, pos = BED$pos, cluster = corrcluster, maxGap = 1000, permbeta = permB, names = rownames(BED), drop = TRUE, mcores = 4)
## test new data.table create method
## time elapsed :  443.170 secs ~ 16 cores ~ 2 mins
corregion_perm <- regionSeeker(beta = B, chr = BED$chr, pos = BED$pos, cluster = corrcluster, maxGap = 1000, permbeta = permB, names = rownames(BED), drop = TRUE, mcores = 4)

## calculate out combine p-value regions
## time on 16 cores = 70 secs/10,000 regions; ~ Total 3000 secs = 50 mins
## cutoff increase, segments goes down. times -7 mins / -0.1 cutoff
## Final total processing time = 5927.57 secs ~ 1.6 hrs and time elasped in combp = 2798.003
## This implies me that tabulate is also a time-consuming step.
combpregion <- combine.pvalue(dat.m = mset.ComBat, pvalues = p, cluster = corrclusters, chr = BED$chr, pos = BED$pos, names = rownames(BED), method = "spearman", combine = "stouffer_liptak")

## smooth beta
betas <- smoother(beta = betar, pos = BED$pos, names = rownames(BED), cluster = corrcluster, weight = sigmar, method = "weightedLowess")
