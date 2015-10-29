#(1)
library(ShortRead) 
library(yeastRNASeq)
fastqFilePath <-  system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
fqFile <- FastqFile(fastqFilePath)
reads <- readFastq(fqFile)
reads_set <- sread(reads)
sum(DNAStringSet(reads_set,5,5) == "A") / length(reads_set) # 0.3638

#(2)
qm <- as(quality(reads), "matrix")
mean(qm[,5:5]) # 28.93

#(3)
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
# In this interval, how many reads are duplicated by position?
bamFile <- BamFile(bamFilePath)
seqinfo(bamFile)
aln <- scanBam(bamFile)
aln <- aln[[1]]  
names(aln)
lapply(aln, function(xx) xx[1])
unique(aln$rname)

gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = 800000, end = 801000))
params <- ScanBamParam(which = gr, what = scanBamWhat())
aln <- scanBam(bamFile, param = params)
aln <- aln[[1]]  
aln$pos # 327 total 
duplicatedValues = unique(aln$pos[duplicated(aln$pos)]) # 49 positions and their names
sum(aln$pos %in% duplicatedValues)  # 129

## another way to do it..
table(table(aln$pos))  # shows 198 values with only 1 position, and 49 with > 1 position
length(aln$pos) - 198 # 129 positions overall with 1 or more 

#(4)
library(GenomicRanges)
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = 807762, end = 808068))
# Scchr13:807762-808068
## ----BamViews------------------------------------------------------------
bamView <- BamViews(bpaths)
bamRanges(bamView) <- gr
aln <- scanBam(bamView)
names(aln)
names(aln[[1]])

lens <- list()
for(i in 1:length(aln)) {
  lens[i] <- length(aln[[i]][[1]]$seq)
}

mean(unlist(lens)) # 90.25

#(5)
library(oligo) # affy and nimblegen arrays GE and snp arrays
library(GEOquery)
geoMat <- getGEO("GSE38792")
pD.all <- pData(geoMat[[1]])
getGEOSuppFiles("GSE38792") # remember, raw data is in supplementary
list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL") # samples from control to samples with sleep apnea
celfiles <- list.files("GSE38792/CEL", full = TRUE)

rawData <- read.celfiles(celfiles)
rawData  # pd.hugene.1.0.st.v1  human gene vs 1 based on random priming

filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),"OSA", "Control")
pData(rawData)

normData <- rma(rawData)
expr <- exprs(normData)
mean(expr["8149273",1:8])  # 7.0218

#(6)
library(limma)
design <- model.matrix(~ normData$group)
fit <- lmFit(normData, design)
fit <- eBayes(fit)
topTable(fit)
abs(topTable(fit, n=1)$logFC)

#(7)
topTable(fit, p.value = 0.05)  # I think that's adj  0

#(8)
library(minfi)
require(minfiData)
data(RGsetEx)
p <- preprocessFunnorm(RGsetEx)
b <- getBeta(p)
is <- getIslandStatus(p)
pData(p)$status
norm <- b[,c(1,2,5)]
can <- b[,c(3,4,6)]
norm_os <- norm[is == "OpenSea",]
can_os <- can[is == "OpenSea",]
mean(norm_os) - mean(can_os)

