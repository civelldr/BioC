## ExpressionSet => eSet (ex: 1 matrix is methylated channel, another matrix is the unmethylated channel)
#
#
# (1)
library(ALL)
library(hgu95av2.db)
data(ALL)
sample_5 <- ALL[,5] # get just sample #5
mean(exprs(sample_5)) # 5.629627

#(2)
library(biomaRt)
library(dplyr)
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
fnames <- featureNames(ALL)
ALL  # says annotation hgu95av2
attributes <- listAttributes(ensembl, page = "feature_page") # there's an annotation called affy_hg_u95av2
query.result <- getBM(attributes = c("affy_hg_u95av2", "ensembl_gene_id", "chromosome_name"),
                      filters = "affy_hg_u95av2", values = fnames, mart = ensembl)

qr <- query.result %>% group_by(affy_hg_u95av2) %>% summarize(gene_count = n()) 
sum(qr$gene_count > 1) # 1045

#(3)
## using query.result from above because I stuck in the chr name

auto.df <- data.frame(chromosome_name = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))
query.result.auto <- merge(query.result, auto.df)  # 12543
dim(unique(query.result.auto))

# should be 11016


#(4)
library(minfiData)
data(MsetEx)
pData(MsetEx)  # 5723646052_R04C01
sample_2 <- MsetEx[,2]  #returning a MethylSet for sample #2
mean(getMeth(sample_2)) # 7228.277

## (5)
library(GEOquery)
eList <- getGEO("GSE788")
class(eList)
length(eList)
names(eList)
eData <- eList[[1]]
pData(eData)
sample_2 <- eData[,2]
mean(exprs(sample_2)) # 756.432

# (6)
library(GenomicRanges)
library(airway)
data(airway)
mean(airway$avgLength) # 113.75

# (7)
sample_3 <- airway[,3]
counts <- assay(sample_3, "counts")
sum(counts >= 1) # 25699

#(8)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exons(txdb)
auto_ucsc <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC", group="auto")
exons <- keepSeqlevels(exons, auto_ucsc)
ncbiStyleLevels <- mapSeqlevels(seqlevels(exons),"NCBI")
exons <- renameSeqlevels(exons, ncbiStyleLevels)
subset <- subsetByOverlaps(airway, exons) 
dim(subset) # 26276

#(9) 
colnames(airway)  # SRR1039508 is sample 1
sample_1 <- airway[,1]
subset_sample_1 <- subsetByOverlaps(sample_1, exons) # from above (8)
counts <- assay(sample_1, "counts")
subset_sample_1_counts <- assay(subset_sample_1, "counts")
sum(subset_sample_1_counts) / sum(counts) # 0.9004193

#(10.3)
library(AnnotationHub)
ah <- AnnotationHub()
# qah_h1 <- query(ah, c("E096", "H3K4me3"))
h1 <- qah_h1[["AH30596"]] # AH30596 | E096-H3K4me3.narrowPeak.gz 
h1 <- keepSeqlevels(h1, auto_ucsc)
h1 <- renameSeqlevels(h1, ncbiStyleLevels)

t <- range(rowRanges(subset_sample_1))
auto_ncbi <- extractSeqlevelsByGroup(species="Homo sapiens", style="NCBI", group="auto")
t <- keepSeqlevels(t, auto_ncbi)
p <- promoters(t)

ov <- subsetByOverlaps(p, h1)
t2 <- subsetByOverlaps(sample_1, ov)
counts <- assay(t2, "counts")
median(counts)

## should be 232 

#(10) dividing subset_sample_1 into expressed and non-expressed.  expr should be marked by H3K4me3 at promoter
# Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package
library(AnnotationHub)
ah <- AnnotationHub()
# qah_h1 <- query(ah, c("E096", "H3K4me3"))
h1 <- qah_h1[["AH30596"]] # AH30596 | E096-H3K4me3.narrowPeak.gz 
h1 <- keepSeqlevels(h1, auto_ucsc)
h1 <- renameSeqlevels(h1, ncbiStyleLevels)

transcripts <- transcripts(txdb)
auto_ucsc <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC", group="auto")
transcripts <- keepSeqlevels(transcripts, auto_ucsc)
ncbiStyleLevels <- mapSeqlevels(seqlevels(transcripts),"NCBI")
transcripts <- renameSeqlevels(transcripts, ncbiStyleLevels)
p <- promoters(transcripts)

subset_sample_1 <- subsetByOverlaps(sample_1, transcripts) 
subset_sample_1 <- subsetByOverlaps(subset_sample_1, h1)

subset_sample_1_counts <- assay(subset_sample_1, "counts")
median(subset_sample_1_counts)
median(subset_sample_1_counts[subset_sample_1_counts > 0])

# (10.2)
library(AnnotationHub)
ah <- AnnotationHub()
qah_h1 <- query(ah, c("E096", "H3K4me3"))
h1 <- qah_h1[["AH30596"]] # AH30596 | E096-H3K4me3.narrowPeak.gz 
auto_ucsc <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC", group="auto")
h1 <- keepSeqlevels(h1, auto_ucsc)
h1 <- promoters(h1)
h1 <- renameSeqlevels(h1, ncbiStyleLevels)

tx_by_gn <- transcriptsBy(txdb, by="gene")
unlisted <- unlist(tx_by_gn)
TSS <- ifelse(strand(unlisted) == "+", start(unlisted), end(unlisted))
TSS <- GRanges(seqnames(unlisted), IRanges(TSS, width=1), strand(unlisted))
TSS_by_gn <- relist(TSS, tx_by_gn)
mcols(TSS) <- mcols(unlisted)
TSS_by_gn <- relist(TSS, tx_by_gn) 
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