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

