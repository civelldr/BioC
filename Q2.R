###
library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

validBaseCount <- letterFrequency(Hsapiens$chr22, "A") +
                  letterFrequency(Hsapiens$chr22, "C") +
                  letterFrequency(Hsapiens$chr22, "G") +
                  letterFrequency(Hsapiens$chr22, "T")
                   
letterFrequency(Hsapiens$chr22, "GC") / validBaseCount
