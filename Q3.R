## ExpressionSet => eSet (ex: 1 matrix is methylated channel, another matrix is the unmethylated channel)
#
#
# (1)
sample_5 <- ALL[,5] # get just sample #5
mean(exprs(sample_5))


# (6)
library(GenomicRanges)
library(airway)
data(airway)
mean(airway$avgLength)

# (7)
sample_3 <- airway[,3]
counts <- assay(sample_3, "counts")
sum(counts >= 1)

#(8)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb_auto <- keepSeqlevels(txdb, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                       "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                       "chr20", "chr21", "chr22"))
exons <- exons(txdb_auto)
subset <- subsetByOverlaps(airway, exons)  # problem is that exons use "chr22" and airway is only "22" so need to remember how to convert


