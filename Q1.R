library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
###

# In this assessment we will study features of the human genome, version "hg19".
# We will only consider data which has been mapped to the autosomes (chr 1 to 22).
# We will use data from the H1 cell line as assayed and quantified by the Roadmap Epigenomics project. 
# The Roadmap Epigenomics project code for the H1 cell line is E003. 
# For histone modification data, the Roadmap project makes several types of quantification available. 
# We will use the co-called "narrowPeak" quantification.
# # 
# Note: This entire quiz is one long analysis, so later questions refer to results generated in earlier questions. 
# Some biology Bivalent chromatin is marked by a combination of active and repressive histone marks. 
# A number of slightly different definitions exists; we will say a region is bivalent if it is enriched in both H3K4me3 and H3K27me3. 
# Note that histone modification marks does not have a strand.
# # 
# Bivalent chromatin has especially been considered in embryonic stem cells. 
# An example of such a cell is the ENCODE Tier 1 cell line called H1. The Roadmap Epigenomics id for this cell line is "E003".
# # 
# We will examine the relationship between bivalent chromatin and CpG Islands. 
# CpG Islands are clusters of many CpG (this is juts CG dinucleotides). 
# Several definitions exists of what is an “CpG Island”; we will use the UCSC definition. 
# Because the CG dinucleotide is its own reverse complement, a CpG cluster exists on the forward strand if and only if it exists on the reverse strand. 
# In other words, CpG Islands does not have a strand.

ah <- AnnotationHub()
qah_h1 <- query(ah, c("E003", "H3K4me3"))
qah_h2 <- query(ah, c("E003", "H3K27me3"))

h1 <- qah_h1[[2]]
h2 <- qah_h2[[2]]
h1 <- dropSeqlevels(h1, c("chrX", "chrY"))
h2 <- dropSeqlevels(h2, c("chrX", "chrY"))

dropSeqlevels(h2, "chrY")
