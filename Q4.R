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

