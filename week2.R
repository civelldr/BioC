## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(Biostrings)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Biostrings"))

## ----DNAString, error=TRUE-----------------------------------------------
dna1 <- DNAString("ACGT-N")
dna1
DNAStringSet("ADE")
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTccA"))
dna2

## ----IUPAC---------------------------------------------------------------
IUPAC_CODE_MAP

## ----DNAStringSubset-----------------------------------------------------
dna1[2:4]
dna2[2:3]

## ----DNAStringSubset2----------------------------------------------------
dna2[[2]] # how to get a string from a set

## ----DNAStringSetNames---------------------------------------------------
names(dna2) <- paste0("seq", 1:3)
dna2

## ----basicFunc-----------------------------------------------------------
width(dna2)
sort(dna2)
rev(dna2)
rev(dna1)

## ----bioFunc-------------------------------------------------------------
translate(dna2)
reverseComplement(dna1)

## ----counting------------------------------------------------------------
alphabetFrequency(dna1)
alphabetFrequency(dna2)
letterFrequency(dna2, "GC")
consensusMatrix(dna2, as.prob = TRUE)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

##
##### BSGenome 
##
## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("BSgenome", "BSgenome.Scerevisiae.UCSC.sacCer2"))

## ----BSgenome------------------------------------------------------------
available.genomes()
library("BSgenome.Scerevisiae.UCSC.sacCer2")
Scerevisiae

## ----BSgenomeLength------------------------------------------------------
seqlengths(Scerevisiae)
seqnames(Scerevisiae)

## ----BSgenomeLoad--------------------------------------------------------
Scerevisiae$chrI

## ----gcChrI--------------------------------------------------------------
letterFrequency(Scerevisiae$chrI, "CG", as.prob = TRUE)

## ----gcGenome------------------------------------------------------------
param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency)
head(bsapply(param, letters = "GC"))

## ----gcGenome2-----------------------------------------------------------
param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency, simplify = TRUE)
bsapply(param, letters = "GC")

## ----gcGenome3-----------------------------------------------------------
sum(bsapply(param, letters = "GC")) / sum(seqlengths(Scerevisiae))

library(Biostrings)
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)


## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Biostrings", "BSgenome",
##            "BSgenome.Scerevisiae.UCSC.sacCer2", "AnnotationHub"))

## ----mmatchPattern-------------------------------------------------------
dnaseq <- DNAString("ACGTACGT")
matchPattern(dnaseq, Scerevisiae$chrI)
countPattern(dnaseq, Scerevisiae$chrI)
vmatchPattern(dnaseq, Scerevisiae)
head(vcountPattern(dnaseq, Scerevisiae))

## ----revCompCheck--------------------------------------------------------
dnaseq == reverseComplement(dnaseq)


## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
library(AnnotationHub)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("BSgenome",
##            "BSgenome.Scerevisiae.UCSC.sacCer2", "AnnotationHub"))

## ----views1--------------------------------------------------------------
library("BSgenome.Scerevisiae.UCSC.sacCer2")
dnaseq <- DNAString("ACGTACGT")
vi <- matchPattern(dnaseq, Scerevisiae$chrI)
vi

## ----views2--------------------------------------------------------------
ranges(vi)

## ----views3--------------------------------------------------------------
vi
Scerevisiae$chrI[ start(vi):end(vi) ]

## ----views4--------------------------------------------------------------
alphabetFrequency(vi)

## ----views5--------------------------------------------------------------
shift(vi, 10)

## ----viewsVMatchPattern--------------------------------------------------
gr <- vmatchPattern(dnaseq, Scerevisiae)
vi2 <- Views(Scerevisiae, gr)

## ----annotationHub-------------------------------------------------------
ahub <- AnnotationHub()
qh <- query(ahub, c("sacCer2", "genes"))
qh
genes <- qh[[which(qh$title == "SGD Genes")]]
genes

## ----promoterGCcontent---------------------------------------------------
prom <- promoters(genes)
head(prom, n = 3)

## ----promoterGCcontent2--------------------------------------------------
prom <- trim(prom)
promViews <- Views(Scerevisiae, prom)
gcProm <- letterFrequency(promViews, "GC", as.prob = TRUE)
head(gcProm)

## ----genomeGC------------------------------------------------------------
params <- new("BSParams", X = Scerevisiae, FUN = letterFrequency, simplify = TRUE)
gccontent <- bsapply(params, letters = "GC")
gcPercentage <- sum(gccontent) / sum(seqlengths(Scerevisiae))
gcPercentage

## ----plotGC, fig=TRUE----------------------------------------------------
plot(density(gcProm))
abline(v = gcPercentage, col = "red")


## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomicRanges"))

## ----RleEx1--------------------------------------------------------------
rl <- Rle(c(1,1,1,1,2,2,3,3,2,2))
rl
runLength(rl)
runValue(rl)
as.numeric(rl)

## ----aggregate-----------------------------------------------------------
ir <- IRanges(start = c(2,6), width = 2)
aggregate(rl, ir, FUN = mean)

## ----coverage------------------------------------------------------------
ir <- IRanges(start = 1:10, width = 3)
rl <- coverage(ir)
rl

## ----slice---------------------------------------------------------------
slice(rl, 2)

## ----views---------------------------------------------------------------
vi <- Views(rl, start = c(3,7), width = 3)
vi

## ----viewsMean-----------------------------------------------------------
mean(vi)

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
rl <- coverage(gr)
rl

## ----GRangesViews, error=TRUE--------------------------------------------
grView <- GRanges("chr1", ranges = IRanges(start = 2, end = 7))
vi <- Views(rl, grView)

## ----GRangesViews2-------------------------------------------------------
vi <- Views(rl, as(grView, "RangesList"))
vi
vi[[1]]

## ----usecase, eval=FALSE-------------------------------------------------
## bases <- reduce(exons)
## nBases <- sum(width(bases))
## nCoverage <- sum(Views(coverage(reads), bases))
## nCoverage / nBases


## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomicRanges"))

## ----CreateGrangesList---------------------------------------------------
gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4, width = 3))
gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 1:4, width = 3))
gL <- GRangesList(gr1 = gr1, gr2 = gr2)
gL

## ----GRangesAccess-------------------------------------------------------
start(gL)
seqnames(gL)

## ----elementLengths------------------------------------------------------
elementLengths(gL)

## ----endoapply-----------------------------------------------------------
shift(gL, 10)

## ----findOverlaps--------------------------------------------------------
findOverlaps(gL, gr2)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
   


## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

## ----txdb----------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

## ----gr------------------------------------------------------------------
gr <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 999999999))
subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)

subsetByOverlaps(genes(txdb), gr)
subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

## ----transcripts---------------------------------------------------------


## ----exons---------------------------------------------------------------
subsetByOverlaps(exons(txdb), gr)

## ----exonsBy-------------------------------------------------------------
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)

## ----cds-----------------------------------------------------------------
subsetByOverlaps(cds(txdb), gr)
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)

## ----transcriptLengths---------------------------------------------------
subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))

## ----help, eval=FALSE----------------------------------------------------
## ?import
## ?BigWigFile

## ----ahub----------------------------------------------------------------
library(AnnotationHub)
ahub <- AnnotationHub()
table(ahub$rdataclass)

## ----granges-------------------------------------------------------------
ahub.gr <- subset(ah, rdataclass == "GRanges" & species == "Homo sapiens")
gr <- ahub.gr[[1]]
gr
seqinfo(gr)

## ----BigWig--------------------------------------------------------------
ahub.bw <- subset(ahub, rdataclass == "BigWigFile" & species == "Homo sapiens")
ahub.bw
bw <- ahub.bw[[1]]
bw

## ----importBigWig--------------------------------------------------------
gr1 <- gr[1:3]
out.gr <- import(bw, which = gr1)
out.gr

## ----importBigWig2-------------------------------------------------------
out.rle <- import(bw, which = gr1, as = "Rle")
out.rle

## ----importBigWig3-------------------------------------------------------
gr.chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = seqlengths(gr)["chr22"]))
out.chr22 <- import(fc.sig, which = gr.chr22, as = "Rle")
out.chr22[["chr22"]]

## ----liftOver------------------------------------------------------------
ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
query(ahub.chain, c("hg18", "hg19"))
chain <- ahub.chain[ahub.chain$title == "hg19ToHg18.over.chain.gz"]
chain <- chain[[1]]
gr.hg18 <- liftOver(gr, chain)
gr.hg18

## ----liftOver2-----------------------------------------------------------
table(elementLengths(gr.hg18))

## ----tabixIndex----------------------------------------------------------
library(Rsamtools)
from <- system.file("extdata", "ex1.sam", package="Rsamtools",
                    mustWork=TRUE)
from
to <- tempfile()
zipped <- bgzip(from, to)
idx <- indexTabix(zipped, "sam")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


