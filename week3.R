## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(Biobase)
library(ALL)
library(hgu95av2.db)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Biobase", "ALL", "hgu95av2.db"))

## ----ALL-----------------------------------------------------------------
library(ALL)
data(ALL)
ALL # expressionset 

## ----help, eval=FALSE----------------------------------------------------
?ALL

## ----experimentData------------------------------------------------------
experimentData(ALL)

## ----exprs---------------------------------------------------------------
exprs(ALL)[1:4, 1:4]  # just first 4 rows and cols

## ----names---------------------------------------------------------------
head(sampleNames(ALL))
head(featureNames(ALL))

## ----pData---------------------------------------------------------------
head(pData(ALL))  #covar about the samples

## ----dollar--------------------------------------------------------------
head(pData(ALL)$sex)
head(ALL$sex)  # $ opertor can be used directly on the ExpressionSet

## ----subset--------------------------------------------------------------
ALL[,1:5]  # all 12k features, first 5 samples  $ FIRST DIM FEATURES, SECOND DIM SAMPLES
ALL[1:10,] # first 10 features, all samples
ALL[1:10,1:5] # first 10 features, first 5 samples

## ----subset2-------------------------------------------------------------
ALL[, c(3,2,1)]
ALL$sex[c(1,2,3)]
ALL[, c(3,2,1)]$sex

## ----featureData---------------------------------------------------------
featureData(ALL)  # seems to be null most of the time

## ----annotation----------------------------------------------------------
ids <- featureNames(ALL)[1:5]  # these need to be mapped into gene symbols
ids

## ----annotation2---------------------------------------------------------
library(hgu95av2.db)
as.list(hgu95av2ENTREZID[ids])  # this gets the entrez id's for the above ids for the features


## ----varLabels-----------------------------------------------------------
pD <- phenoData(ALL)  ## so basically just use pData!  i.e. pData(ALL[,5])  (pData for sample5)
varLabels(pD)

## ----varLabels2----------------------------------------------------------
varLabels(pD)[2] <- "Age at diagnosis"
pD
colnames(pD)[1:3]
varLabels(pD)[1:3]

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

###
##### SummarizedExperiment (newer version of expressionset)
###

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicRanges)
library(airway)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomicRanges", "airway"))

## ----airway--------------------------------------------------------------
library(airway)
data(airway)
airway

## ----colData-------------------------------------------------------------
colData(airway)  # this returns a DataFrame and was pData in ExpressionSet

## ----getColumn-----------------------------------------------------------
airway$cell

## ----exptData- Experiment data----------------------------------------------
exptData(airway)
as.list(exptData(airway))

## ----names---------------------------------------------------------------
colnames(airway)
head(rownames(airway))

## ----assay---------------------------------------------------------------
airway
assayNames(airway)  #there is 1 and it's called "counts"
assays(airway)
head(assay(airway, "counts"))[1:4, 1:4]  # first 4 features and samples

## ----rowRanges-----------------------------------------------------------
length(rowRanges(airway))  # one granges for each feature b/c each row is a gene
rowRanges(airway) # returns a list, and each element in the list is a list of exons for a gene

dim(airway)
rowRanges(airway)

## ----numberOfExons-------------------------------------------------------
length(rowRanges(airway))
sum(elementLengths(rowRanges(airway)))  # number of exons b/c each grange are the exons for the features

## ----start---------------------------------------------------------------
start(rowRanges(airway))
start(airway)

## ----subsetByOverlaps----------------------------------------------------
gr <- GRanges(seqnames = "1", ranges = IRanges(start = 1, end = 10^7))  # say a 10 MB region
subsetByOverlaps(airway, gr) # genes inside airways that overlap the above gr

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GEOquery)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GEOquery"))

## ----getData-------------------------------------------------------------
eList <- getGEO("GSE11675")
class(eList)
length(eList)
names(eList)
eData <- eList[[1]]
eData

## ----pData---------------------------------------------------------------
names(pData(eData))

## ----getGEOsupp----------------------------------------------------------
eList2 <- getGEOSuppFiles("GSE11675")
eListeList2
tarArchive <- rownames(eList2)[1]
tarArchive

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(biomaRt)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("biomaRt"))

## ----listMarts-----------------------------------------------------------
head(listMarts())
mart <- useMart("ensembl")
mart
head(listDatasets(mart))
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
ensembl

## ----getBMex-------------------------------------------------------------
values <- c("202763_at","209310_s_at","207500_at")

## ----listAttributes------------------------------------------------------
attributes <- listAttributes(ensembl)
head(attributes)
nrow(attributes)
filters <- listFilters(ensembl)
head(filters)
nrow(filters)

# i.e. just read browseVignettes("biomaRt")

## ----listPages-----------------------------------------------------------
attributePages(ensembl)
attributes <- listAttributes(ensembl, page = "feature_page")
head(attributes)
nrow(attributes)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(ALL)
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("ALL" "GenomicRanges"))

## ----lm------------------------------------------------------------------
df <- data.frame(y = rnorm(10), x = rnorm(10))
lm.object <- lm(y ~ x, data = df)
lm.object
names(lm.object)
class(lm.object)

## ----lm2-----------------------------------------------------------------
xx <- list(a = letters[1:3], b = rnorm(3))
xx
class(xx) <- "lm"
xx

## ----ALL-----------------------------------------------------------------
library(ALL)
data(ALL)
ALL
class(ALL)
isS4(ALL)

## ----help, eval=FALSE----------------------------------------------------
## ?"ExpressionSet-class"
## class?ExpressionSet

## ----list----------------------------------------------------------------
xx <- list(a = 1:3)

## ----ExpressionSet-------------------------------------------------------
ExpressionSet()

## ----help2,eval=FALSE----------------------------------------------------
## ?ExpressionSet

## ----newExpressionSet----------------------------------------------------
new("ExpressionSet")

## ----getClass------------------------------------------------------------
getClass("ExpressionSet")

## ----slots---------------------------------------------------------------
ALL@annotation
slot(ALL, "annotation")

## ----accessor------------------------------------------------------------
annotation(ALL)

## ----updateObject, eval=FALSE--------------------------------------------
## new_object <- updateObject(old_object)

## ----updateObject2, eval=FALSE-------------------------------------------
## object <- updateObject(object)

## ----validity------------------------------------------------------------
validObject(ALL)

## ----mimicMethod---------------------------------------------------------
mimicMethod <- function(x) {
  if (is(x, "matrix"))
    method1(x)
  if (is(x, "data.frame"))
    method2(x)
  if (is(x, "IRanges"))
    method3(x)
}

## ----as.data.frame-------------------------------------------------------
as.data.frame

## ----showMethods---------------------------------------------------------
showMethods("as.data.frame")

## ----getMethod-----------------------------------------------------------
getMethod("as.data.frame", "DataFrame")

## ----base_as.data.frame--------------------------------------------------
base::as.data.frame

## ----helpMethod,eval=FALSE-----------------------------------------------
## method?as.data.frame,DataFrame
## ?"as.data.frame-method,DataFrame"

## ----findOverlaps--------------------------------------------------------
showMethods("findOverlaps")

## ----ignore.strand-------------------------------------------------------
getMethod("findOverlaps", signature(query = "Ranges", subject = "Ranges"))
getMethod("findOverlaps", signature(query = "GenomicRanges", subject = "GenomicRanges"))

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(ALL)
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("ALL" "GenomicRanges"))

## ----lm------------------------------------------------------------------
df <- data.frame(y = rnorm(10), x = rnorm(10))
lm.object <- lm(y ~ x, data = df)
lm.object
names(lm.object)
class(lm.object)

## ----lm2-----------------------------------------------------------------
xx <- list(a = letters[1:3], b = rnorm(3))
xx
class(xx) <- "lm"
xx

## ----ALL-----------------------------------------------------------------
library(ALL)
data(ALL)
ALL
class(ALL)
isS4(ALL)

## ----help, eval=FALSE----------------------------------------------------
## ?"ExpressionSet-class"
## class?ExpressionSet

## ----list----------------------------------------------------------------
xx <- list(a = 1:3)

## ----ExpressionSet-------------------------------------------------------
ExpressionSet()

## ----help2,eval=FALSE----------------------------------------------------
## ?ExpressionSet

## ----newExpressionSet----------------------------------------------------
new("ExpressionSet")

## ----getClass------------------------------------------------------------
getClass("ExpressionSet")

## ----slots---------------------------------------------------------------
ALL@annotation
slot(ALL, "annotation")

## ----accessor------------------------------------------------------------
annotation(ALL)

## ----updateObject, eval=FALSE--------------------------------------------
## new_object <- updateObject(old_object)

## ----updateObject2, eval=FALSE-------------------------------------------
## object <- updateObject(object)

## ----validity------------------------------------------------------------
validObject(ALL)

## ----mimicMethod---------------------------------------------------------
mimicMethod <- function(x) {
  if (is(x, "matrix"))
    method1(x)
  if (is(x, "data.frame"))
    method2(x)
  if (is(x, "IRanges"))
    method3(x)
}

## ----as.data.frame-------------------------------------------------------
as.data.frame

## ----showMethods---------------------------------------------------------
showMethods("as.data.frame")

## ----getMethod-----------------------------------------------------------
getMethod("as.data.frame", "DataFrame")

## ----base_as.data.frame--------------------------------------------------
base::as.data.frame

## ----helpMethod,eval=FALSE-----------------------------------------------
## method?as.data.frame,DataFrame
## ?"as.data.frame-method,DataFrame"

## ----findOverlaps--------------------------------------------------------
showMethods("findOverlaps")

## ----ignore.strand-------------------------------------------------------
getMethod("findOverlaps", signature(query = "Ranges", subject = "Ranges"))
getMethod("findOverlaps", signature(query = "GenomicRanges", subject = "GenomicRanges"))

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()



