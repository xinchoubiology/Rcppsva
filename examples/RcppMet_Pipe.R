## pipeline dependence package
library(minfi)
library(minfiData)
library(Rcppsva)

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
mod   <- model.matrix(~ status + person, data = pd)

mset.ComBat <- ComBat(dat = M, batch = pd$Slide, mod = mod, par.prior = FALSE)

## Before next step, I did the profiler for ComBat function
## library(lineprof)
## tmp <- lineprof(ComBat(dat = M[1:10000,], batch = pd$Slide, mod = mod, par.prior = FALSE), torture = TRUE)
## shine(tmp)

## choose Combat-n cuz Paper: ChaoChen_PlosOne_2011
## after ComBat, we use SVA to learn all the other unmodel factors
mset.SV <- isvaFn(dat.m = mset.ComBat, design = mod, type = "M", verbose = TRUE)

## calculate DMP(differential methylation position)
pheno.v  <- factor(pd$status, levels = c("normal", "cancer"))
## rank "cancer" followed by "normal" so estimated means cancer - normal
mod   <- model.matrix(~ pheno.v + person, data = pd)

## test differential methylation(siggene) by linear regression on M data(link function from beta is logit)
## generally, siggene found by "NULL" method will be less that those found by "limma" eBayes method
mset.DMP <- svaReg(dat.m = mset.ComBat, design = mod, sv.m = mset.SV$isv, qvalue0 = 0.05, backend = "NULL", verbose = TRUE)
mset.DMP.limma <- svaReg(dat.m = mset.ComBat, design = mod, sv.m = mset.SV$isv, qvalue0 = 0.05, backend = "limma", verbose = TRUE)

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
