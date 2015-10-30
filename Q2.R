library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

validBaseCount <- letterFrequency(Hsapiens$chr22, "A") +
                  letterFrequency(Hsapiens$chr22, "C") +
                  letterFrequency(Hsapiens$chr22, "G") +
                  letterFrequency(Hsapiens$chr22, "T")
                   
letterFrequency(Hsapiens$chr22, "GC") / validBaseCount

# 0.4798807

##### (2)
library(AnnotationHub)
ah <- AnnotationHub()
qah_h1 <- query(ah, c("E003", "H3K27me3"))
h1 <- qah_h1[[2]]
seqlevels(h1)
h1 <- keepSeqlevels(h1, c("chr22"))

v1 <- Views(Hsapiens, h1)
mean(letterFrequency(v1, "GC", as.prob = T))

# 0.528866

####### (3) 

sv <- mcols(v1)$signalValue
gcFreq <- letterFrequency(v1, "GC", as.prob = T)
cor(sv, gcFreq)
# 0.004467924

##### (4)
fc.sig <- qah_h1[[4]] # this is the bigwig file
gr.chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 51304566))
rle.allChr <- import(fc.sig, which = gr.chr22, as = "Rle")
rle.chr22 <- rle.allChr$chr22
fc.sig.view <- Views(rle.chr22, start=start(h1), end=end(h1))
fc.sig.mean <- mean(fc.sig.view)
cor(fc.sig.mean, sv)

# 0.9149614

### (5)

sum(rle.chr22 >= 1)

# 10914671

####### (6)

qah_h2 <- query(ah, c("E055", "H3K27me3"))
fc.sig.E055 <- qah_h2[[4]]
gr.chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 51304566))
rle.E055.allChr <- import(fc.sig.E055, which = gr.chr22, as = "Rle")
rle.E003.allChr <- rle.allChr 

 #  Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.
E003_eval <- slice(rle.E003.allChr[["chr22"]],upper=0.5)
E055_eval <- slice(rle.E055.allChr[["chr22"]],lower=2)

E003_ir <- as(E003_eval, "IRanges")
E055_ir <- as(E055_eval, "IRanges")

comb_ir <- intersect(E003_ir, E055_ir)

sum(width(comb_ir))

# 1869937

####### (7) (new q10)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb, force=TRUE) <- c("chr22")
gr <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 52330658))
gr.trans.chr22 <- subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
length(gr.trans.chr22) 
gr.prom <- promoters(gr.trans.chr22, upstream = 900, downstream = 100)
tl.chr22 <- transcriptLengths(txdb, with.cds_len = TRUE) #rtn df
tl.chr22  <- tl.chr22[tl.chr22$cds_len > 0,]
trans.eval <- gr.prom[mcols(gr.prom)$tx_id %in% tl.chr22$tx_id]
sum(coverage(trans.eval) > 1)

# 306920

library(AnnotationHub)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb, force=TRUE) <- c("chr22")
ahub <- AnnotationHub()
ahub.gr <- subset(ahub, rdataclass == "GRanges" & species == "Homo sapiens")
gr <- ahub.gr[[1]]

gr <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = seqlengths(gr)["chr22"]))
gr.trans.chr22 <- subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
length(gr.trans.chr22)
gr.prom <- promoters(gr.trans.chr22, upstream = 900, downstream = 100)
tl.chr22 <- transcriptLengths(txdb, with.cds_len = TRUE) #rtn df
tl.chr22  <- tl.chr22[tl.chr22$cds_len > 0,]
trans.eval <- gr.prom[mcols(gr.prom)$tx_id %in% tl.chr22$tx_id]
sum(coverage(trans.eval) > 1)


### new (7)

ah <- AnnotationHub()
ah <- subset(ah, genome == "hg19")
ah <- subset(ah, dataprovider == "UCSC")
# write.csv(mcols(ah), "ah.csv")
# g <- ah[["AH5018"]] # assembly
cpg <- ah[["AH5086"]] # CpG islands
cpg_chr22 <- keepSeqlevels(cpg, "chr22")
v_cpg_chr22 <- Views(Hsapiens, cpg_chr22)
# gchr22 <- keepSeqlevels(g, "chr22")
# Hschr22 <- Hsapiens$chr22 # actually chr22 sequence

# for some reason the as.prob isn't working...
oe <- dinucleotideFrequency(v_cpg_chr22, as.prob = T)[,7] / (letterFrequency(v_cpg_chr22, "C", as.prob = T) * letterFrequency(v_cpg_chr22, "G", as.prob = T))
# mean(oe)

oe2 <- (dinucleotideFrequency(v_cpg_chr22)[,7] / width(v_cpg_chr22)) /  
  (letterFrequency(v_cpg_chr22, "C") / width(v_cpg_chr22)  * letterFrequency(v_cpg_chr22, "G") / width(v_cpg_chr22))
mean(oe2)

### (8)
dnaseq <- DNAString("TATAAA")
countPattern(dnaseq, Hsapiens$chr22) + countPattern(dnaseq, reverseComplement(Hsapiens$chr22))

### (9)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb, force=TRUE) <- c("chr22")
# trans.eval  <- from earlier
v_prom <- Views(Hsapiens$chr22, start=start(trans.eval), end=end(trans.eval))
countPattern(dnaseq, v_prom)
# beats me but answer within 2

##(10)

