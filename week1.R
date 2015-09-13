## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(IRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("IRanges"))

## ----iranges1------------------------------------------------------------
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir1
ir2 <- IRanges(start = c(1,3,5), width = 3)
all.equal(ir1, ir2)

## ----ir_width------------------------------------------------------------
start(ir1)
width(ir2) <- 1
ir2

## ----ir_names------------------------------------------------------------
names(ir1) <- paste("A", 1:3, sep = "")
ir1

## ----ir_dim--------------------------------------------------------------
dim(ir1)
length(ir1)

## ----ir_subset-----------------------------------------------------------
ir1[1]
ir1["A1"]

## ----concatenate---------------------------------------------------------
c(ir1, ir2)

## ----irNormal1, echo=FALSE-----------------------------------------------
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))

## ----irNormal2, echo=FALSE, fig.height=2, small.mar=TRUE-----------------
plotRanges(ir)

## ----irNormal3, echo=FALSE, fig.height=1.75, small.mar=TRUE--------------
plotRanges(reduce(ir))

## ----irNormal4-----------------------------------------------------------
ir
reduce(ir)

## ----irDisjoin1, eval=FALSE----------------------------------------------
## disjoin(ir1)

## ----irDisjoin2, echo=FALSE, fig.height=2, small.mar=TRUE----------------
plotRanges(ir)

## ----irDisjoin3, echo=FALSE, fig.height=1.75, small.mar=TRUE-------------
plotRanges(disjoin(ir))

## ----ir_resize-----------------------------------------------------------
resize(ir, width = 1, fix = "start")
resize(ir, width = 1, fix = "center")

## ----ir_sets-------------------------------------------------------------
ir1 <- IRanges(start = c(1, 3, 5), width = 1)
ir2 <- IRanges(start = c(4, 5, 6), width = 1)
union(ir1, ir2)
intersect(ir1, ir2)

## ----union2--------------------------------------------------------------
reduce(c(ir1, ir2))

## ----findOverlaps--------------------------------------------------------
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3)
ov <- findOverlaps(ir1, ir2)
ov

## ----findOverlaps_ill----------------------------------------------------
intersect(ir1[subjectHits(ov)[1]],
          ir2[queryHits(ov)[2]])

## ----subjectHits---------------------------------------------------------
queryHits(ov)
unique(queryHits(ov))

## ----argsFindOverlaps, tidy=TRUE-----------------------------------------
args(findOverlaps)

## ----countOverlaps-------------------------------------------------------
countOverlaps(ir1, ir2)

## ----nearest-------------------------------------------------------------
ir1
ir2
nearest(ir1, ir2)


######### GenomicRanges ###########
##
## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomeInfoDb)
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomeInfoDb", "GenomicRanges"))

## ----seqlevelsForce------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
seqlevels(gr, force=TRUE) <- "chr1"
gr

## ----dropSeqlevels-------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
dropSeqlevels(gr, "chr1")
keepSeqlevels(gr, "chr2")

## ----keepStandard--------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chrU345"),
              ranges = IRanges(start = 1:2, end = 4:5))
keepStandardChromosomes(gr)

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 2))

## ----seqStyle------------------------------------------------------------
newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
gr <- renameSeqlevels(gr, newStyle)

## ----DataFrame-----------------------------------------------------------
ir <- IRanges(start = 1:2, width = 3)
df1 <- DataFrame(iranges = ir)
df1
df1$iranges
df2 <- data.frame(iranges = ir)
df2

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
values(gr) <- DataFrame(score = c(0.1, 0.5, 0.3))
gr

## ----grdollar------------------------------------------------------------
gr$score
gr$score2 = gr$score * 0.2
gr

## ----findOverlaps_setup--------------------------------------------------
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*",
               ranges = IRanges(start = c(1, 3, 5), width = 3))
gr2
gr

## ----findOverlaps--------------------------------------------------------
findOverlaps(gr, gr2)

## ----subsetByOverlaps----------------------------------------------------
subsetByOverlaps(gr, gr2)

## ----makeGRangesFromDataFrame--------------------------------------------
df <- data.frame(chr = "chr1", start = 1:3, end = 4:6, score = 7:9)
makeGRangesFromDataFrame(df)
makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomeInfoDb)
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomeInfoDb", "GenomicRanges"))

## ----seqlevelsForce------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
seqlevels(gr, force=TRUE) <- "chr1"
gr

## ----dropSeqlevels-------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
dropSeqlevels(gr, "chr1")
keepSeqlevels(gr, "chr2")

## ----keepStandard--------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chrU345"),
              ranges = IRanges(start = 1:2, end = 4:5))
keepStandardChromosomes(gr)

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 2))

## ----seqStyle------------------------------------------------------------
newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
gr <- renameSeqlevels(gr, newStyle)

##
###### AnnotationHub GRanges USE CASE
###

library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomicRanges", "rtracklayer", "AnnotationHub"))

## ----ahub_species--------------------------------------------------------
ah <- AnnotationHub()
ah <- subset(ah, species == "Homo sapiens")

## ----ahub_histone--------------------------------------------------------
qhs <- query(ah, "H3K4me3")
qhs <- query(qhs, "Gm12878")

## ----ahub_look-----------------------------------------------------------
qhs

## ----ahub_closerlook-----------------------------------------------------
qhs$title
qhs$dataprovider

## ----ahub_twoGR----------------------------------------------------------
gr1 <- subset(qhs, title == "wgEncodeUwHistoneGm12878H3k4me3StdPkRep1.narrowPeak.gz")[[1]]
gr1
gr2 <- subset(qhs, title == "E116-H3K4me3.narrowPeak.gz")[[1]]
gr2

## ----ahub_summary--------------------------------------------------------
summary(width(gr1))
table(width(gr1))
summary(width(gr2))

## ----ahub_refseq---------------------------------------------------------
qhs <- query(ah, "RefSeq")
qhs

## ----ahub_refseq_genome--------------------------------------------------
qhs$genome

## ----ahub_histone_genome-------------------------------------------------
genome(gr1)

## ----ahub_get_refseq-----------------------------------------------------
refseq <- qhs[qhs$genome == "hg19" & qhs$title == "RefSeq Genes"]
refseq
refseq <- refseq[[1]] ## Downloads

## ----refseq--------------------------------------------------------------
refseq

## ----ahub_refseq_name----------------------------------------------------
table(table(refseq$name))

## ----promoters-----------------------------------------------------------
promoters <- promoters(refseq)
table(width(promoters))
args(promoters)

## ----findOverlaps--------------------------------------------------------
ov <- findOverlaps(promoters, gr1)
ov

## ----queryHits-----------------------------------------------------------
length(unique(queryHits(ov))) / length(gr1)

## ----subjectHits---------------------------------------------------------
length(unique(subjectHits(ov))) / length(promoters)

## ----widthPercentage-----------------------------------------------------
sum(width(reduce(gr1))) / 10^6
sum(width(reduce(promoters))) / 10^6

## ----size----------------------------------------------------------------
sum(width(intersect(gr1, promoters))) / 10^6

## ----size2---------------------------------------------------------------
sum(width(intersect(gr1, promoters, ignore.strand = TRUE))) / 10^6

## ----widthPercentage2----------------------------------------------------
sum(width(reduce(promoters))) / 10^6
sum(width(reduce(promoters, ignore.strand = TRUE))) / 10^6

## ----promInOut-----------------------------------------------------------
prom <- reduce(promoters, ignore.strand = TRUE)
peaks <- reduce(gr1)
both <- intersect(prom, peaks)
only.prom <- setdiff(prom, both)
only.peaks <- setdiff(peaks, both)
overlapMat <- matrix(0,ncol = 2, nrow = 2)
colnames(overlapMat) <- c("in.peaks", "out.peaks")
rownames(overlapMat) <- c("in.promoters", "out.promoter")
overlapMat[1,1] <- sum(width(both))
overlapMat[1,2] <- sum(width(only.prom))
overlapMat[2,1] <- sum(width(only.peaks))
overlapMat[2,2] <- 3*10^9 - sum(overlapMat)
round(overlapMat / 10^6, 2)

## ----oddsratio-----------------------------------------------------------
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio

## ----oddsRatio2----------------------------------------------------------
overlapMat[2,2] <- 1.5*10^9
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio


