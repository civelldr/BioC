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
exons <- exons(txdb)
auto_ucsc <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC", group="auto")
exons <- keepSeqlevels(exons, auto_ucsc)
ncbiStyleLevels <- mapSeqlevels(seqlevels(exons),"NCBI")
exons <- renameSeqlevels(exons, newStyle)
subset <- subsetByOverlaps(airway, exons) 
dim(subset)

#(9) 
colnames(airway)  # SRR1039508 is sample 1
sample_1 <- airway[,1]
subset_sample_1 <- subsetByOverlaps(sample_1, exons) # from above (8)
counts <- assay(sample_1, "counts")
subset_sample_1_counts <- assay(subset_sample_1, "counts")
sum(subset_sample_1_counts) / sum(counts)

#(10) dividing subset_sample_1 into expressed and non-expressed.  expr should be marked by H3K4me3 at promoter
# Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package
library(AnnotationHub)
ah <- AnnotationHub()
qah_h1 <- query(ah, c("E096", "H3K4me3"))
h1 <- qah_h1[["AH30596"]] # AH30596 | E096-H3K4me3.narrowPeak.gz 
h1 <- keepSeqlevels(h1, auto_ucsc)
h1 <- promoters(h1)
h1 <- renameSeqlevels(h1, newStyle)

transcripts <- transcripts(txdb)
auto_ucsc <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC", group="auto")
transcripts <- keepSeqlevels(transcripts, auto_ucsc)
ncbiStyleLevels <- mapSeqlevels(seqlevels(transcripts),"NCBI")
transcripts <- renameSeqlevels(transcripts, newStyle)
p <- promoters(transcripts)

subset_sample_1 <- subsetByOverlaps(sample_1, transcripts) 
subset_sample_1 <- subsetByOverlaps(subset_sample_1, h1)

subset_sample_1_counts <- assay(subset_sample_1, "counts")
median(subset_sample_1_counts)
median(subset_sample_1_counts[subset_sample_1_counts > 0])
