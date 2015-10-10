library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
###

ah <- AnnotationHub()
ah <- subset(ah, genome == "hg19")
ah <- subset(ah, dataprovider == "UCSC")
ah <- query(ah, "CpG Islands")
cpgIslands <- ah[[1]]
cpgAuto <- keepSeqlevels(cpgIslands, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                          "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                          "chr20", "chr21", "chr22"))
cpgChr4 <- keepSeqlevels(cpgIslands, "chr4")
length(cpgAuto)
length(cpgChr4)

###
ah <- AnnotationHub()
qah_h1 <- query(ah, c("E003", "H3K4me3"))
h1 <- qah_h1[[2]]
seqlevels(h1)
h1 <- dropSeqlevels(h1, c("chrX", "chrY"))
sum(width(reduce(h1)))

####
qah_h2 <- query(ah, c("E003", "H3K27me3"))
h2 <- qah_h2[[2]]
h2 <- dropSeqlevels(h2, c("chrX", "chrY"))
mean(h2$signalValue)

##
bivReg <- intersect(h1, h2)
sum(width(reduce(bivReg)))

####
ov <- findOverlaps(bivReg, cpgAuto)
# number of bivReg that have cpg islands in them
length(unique(queryHits(ov))) / length(bivReg) 
# the number of cpg islands that have bivReg in them 
length(unique(subjectHits(ov))) / length(cpgAuto) 
# perc of bivalent regions that overlap a cpgIsland
length(subsetByOverlaps(bivReg, cpgAuto, ignore.strand = TRUE)) / length(bivReg)

# bases which part of cpgAuto that are also part of the bivalent regions 
sum(width(intersect(bivReg, cpgAuto, ignore.strand = TRUE))) / sum(width(reduce(cpgAuto)))

####### How many bases are bivalently marked within 10kb of CpG Islands?
cpgAuto_flank<- resize(cpgAuto,width=20000+width(cpgAuto),fix='center')
ov <- intersect(cpgAuto_flank, bivReg, ignore.strand = TRUE)
sum(width(reduce(ov)))

######
g <- ah[["AH5018"]]  # Assembly
g <- keepSeqlevels(g, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                       "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                       "chr20", "chr21", "chr22"))
gSize <- sum(as.numeric(seqlengths(g))) # hm....
sum(width(reduce(cpgAuto, ignore.strand = T))) / gSize

## odds ratio for the overlap of bivalent marks with CpG islands.
inOut <- matrix(0, ncol = 2, nrow = 2)
rownames(inOut) <- c("in", "out")
colnames(inOut) <- c("in", "out")
# number of items that are in common between both cpg and bivalent regions
inOut[1,1] <- sum(width(intersect(cpgAuto, bivReg, ignore.strand=TRUE)))
inOut[1,2] <- sum(width(setdiff(bivReg, cpgAuto, ignore.strand=TRUE)))
inOut[2,1] <- sum(width(setdiff(cpgAuto, bivReg, ignore.strand=TRUE)))
inOut[2,2] <- gSize - sum(inOut)
or <- inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2]) 
or
# i.e. 169x more enrichment of the bivalent regions contained within CpG Islands
