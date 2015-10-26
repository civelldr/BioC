## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(ShortRead)  # best way to read in FASTQ files

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("ShortRead"))

## ----fastq1--------------------------------------------------------------
fastqDir <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")
fastqDir
fastqPath <- list.files(fastqDir, pattern = ".fastq.gz$", full = TRUE)[1]
reads <- readFastq(fastqPath)
reads

## ----fastq2--------------------------------------------------------------
fqFile <- FastqFile(fastqPath)
fqFile
reads <- readFastq(fqFile)

## ----accessorFastq-------------------------------------------------------
sread(reads)[1:2]
quality(reads)[1:2]
id(reads)[1:2]

## ----convertQual---------------------------------------------------------
qm <- as(quality(reads), "matrix")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(Rsamtools)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Rsamtools")

## ----bamPath-------------------------------------------------------------
bamPath <- system.file("extdata", "ex1.bam", package="Rsamtools")
bamFile <- BamFile(bamPath)
bamFile

## ----bamFileInfo---------------------------------------------------------
seqinfo(bamFile)

## ----scanBam-------------------------------------------------------------
aln <- scanBam(bamFile)
length(aln)
class(aln)

## ----lookAtBam-----------------------------------------------------------
aln <- aln[[1]]
names(aln)
lapply(aln, function(xx) xx[1])  # first element of each name
aln$seq

## ----yieldSize-----------------------------------------------------------
yieldSize(bamFile) <- 1  # means that every time you call scanBam you get one read
open(bamFile)
scanBam(bamFile)[[1]]$seq  # essentially this is iterating over the bamfile
scanBam(bamFile)[[1]]$seq
## Cleanup
close(bamFile)
yieldSize(bamFile) <- NA

## ----ScanBamParams-------------------------------------------------------
gr <- GRanges(seqnames = "seq2",
              ranges = IRanges(start = c(100, 1000), end = c(1500,2000)))
params <- ScanBamParam(which = gr, what = scanBamWhat())
scanBamWhat() # which pieces of the bam file we want to read in
aln <- scanBam(bamFile, param = params)
names(aln)
aln[[1]]$pos
head(aln[[1]]$pos) # reads are long, so the left most overlap the positions of the granges we provided

## ----summary-------------------------------------------------------------
quickBamFlagSummary(bamFile)

## ----BamViews------------------------------------------------------------
bamView <- BamViews(bamPath)
aln <- scanBam(bamView)
names(aln)

## ----BamViews2-----------------------------------------------------------
bamRanges(bamView) <- gr
aln <- scanBam(bamView)
names(aln)
names(aln[[1]])

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(oligo) # affy and nimblegen arrays GE and snp arrays
library(GEOquery)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("oligo", "GEOquery))

## ----getData-------------------------------------------------------------
library(GEOquery)
getGEOSuppFiles("GSE38792") # remember, raw data is in supplementary

# probset group of probes that all measure the same target

list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL") # samples from control to samples with sleep apnea

## ----readData, message=FALSE---------------------------------------------
library(oligo)
celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)

## ----show----------------------------------------------------------------
rawData  # pd.hugene.1.0.st.v1  human gene vs 1 based on random priming

## ----getClass------------------------------------------------------------
getClass("GeneFeatureSet")

## ----rawPeak-------------------------------------------------------------
exprs(rawData)[1:4,1:3] # values indicate that it's raw intensity measurements 16 bit scanner, 0 - 2^16

## ----maxExpr-------------------------------------------------------------
max(exprs(rawData)) # 65534 means there's basically a prob that maxes it out, so if we log2 it we'll get a value between 1 and 16

## ----pData---------------------------------------------------------------
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
pData(rawData)

## ----rawBox, plot=TRUE---------------------------------------------------
boxplot(rawData) # sumary of the dist of of the intensity of each array.  y axis is on the log scale, so a different of 1 or 2 is 
# a massive differences.  what about those little three?  maybe nothing hybridized to them...

## ----rma-----------------------------------------------------------------
normData <- rma(rawData)  # robust multichip average algo, all the probes that measure the same thing have been summarized
normData
featureNames(normData)[1:10]

## ----normBox, plot=TRUE--------------------------------------------------
boxplot(normData)  # same mean, same spread, ready for analysis
exprs(normData)[1:4,1:3]

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(limma) # linear models for microarrays
library(leukemiasEset)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("limma", "leukemiasEset"))

## ----load----------------------------------------------------------------
library(leukemiasEset)
data(leukemiasEset)
leukemiasEset
table(leukemiasEset$LeukemiaType)

## ----subset--------------------------------------------------------------
ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
ourData$LeukemiaType <- factor(ourData$LeukemiaType)

## ----limma---------------------------------------------------------------
design <- model.matrix(~ ourData$LeukemiaType)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
topTable(fit)

## ----level---------------------------------------------------------------
ourData$LeukemiaType

## ----FCbyHand------------------------------------------------------------
topTable(fit, n = 1)
genename <- rownames(topTable(fit, n=1))
typeMean <- tapply(exprs(ourData)[genename,], ourData$LeukemiaType, mean)
typeMean
typeMean["NoL"] - typeMean["ALL"]

## ----design2-------------------------------------------------------------
design <- model.matrix(~ ourData$LeukemiaType)

## ----headDesign----------------------------------------------------------
head(design)

## ----design3-------------------------------------------------------------
design2 <- model.matrix(~ ourData$LeukemiaType - 1)
head(design2)
colnames(design2) <- c("ALL", "NoL")

## ----design4-------------------------------------------------------------
fit2 <- lmFit(ourData, design2)
contrast.matrix <- makeContrasts("ALL-NoL", levels = design2)
contrast.matrix

## ----cont.fit------------------------------------------------------------
fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)
topTable(fit2C)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(minfi)
library(GEOquery)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("minfi", "GEOquery"))

## ----geoquery------------------------------------------------------------
library(GEOquery)
getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

## ----decompress----------------------------------------------------------
idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

## ----readExp-------------------------------------------------------------
rgSet <- read.450k.exp("GSE68777/idat")
rgSet
pData(rgSet)
head(sampleNames(rgSet))

## ----geoPheno------------------------------------------------------------
geoMat <- getGEO("GSE68777")
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)

## ----merge---------------------------------------------------------------
sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- pD
rgSet

## ----preprocess----------------------------------------------------------
grSet <- preprocessQuantile(rgSet)
grSet

## ----granges-------------------------------------------------------------
granges(grSet)

## ----getBeta-------------------------------------------------------------
getBeta(grSet)[1:3,1:3]

## ----getIslandStatus-----------------------------------------------------
head(getIslandStatus(grSet))

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(DESeq2)
library(edgeR)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("DESeq2", "edgeR))

## ----data----------------------------------------------------------------
library(airway)
data(airway)
airway
assay(airway, "counts")[1:3, 1:3]
airway$dex

## ----relevel-------------------------------------------------------------
airway$dex <- relevel(airway$dex, "untrt")
airway$dex

## ----granges-------------------------------------------------------------
granges(airway)

## ----edgeRsetup----------------------------------------------------------
library(edgeR)
dge <- DGEList(counts = assay(airway, "counts"),
               group = airway$dex)
dge$samples <- merge(dge$samples,
                     as.data.frame(colData(airway)),
                     by = 0)
dge$genes <- data.frame(name = names(rowRanges(airway)),
                        stringsAsFactors = FALSE)

## ----calcNormFactors-----------------------------------------------------
dge <- calcNormFactors(dge)

## ----disp----------------------------------------------------------------
design <- model.matrix(~dge$samples$group)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

## ----edgeRdesign---------------------------------------------------------
fit <- glmFit(dge, design)

## ----glmLRT--------------------------------------------------------------
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

## ----DESeq2setup---------------------------------------------------------
library(DESeq2)
dds <- DESeqDataSet(airway, design = ~ dex)

## ----deseqfit------------------------------------------------------------
dds <- DESeq(dds)

## ----deseqResults--------------------------------------------------------
res <- results(dds)
res <- res[order(res$padj),]
res

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

