library(gplots)
library(ggplot2)
library(dplyr)
library(assertthat)
library(sva)

source('analysis.r')
source('metadata.r')

### Get bulk data
bulkTPMData <- read.table('bulk_tpm.txt', header = TRUE)
celltypesToExclude <- getCelltypesToExclude()
bulkExperimentMetadata <- getExperimentMetadata()
bulkExperimentMetadata <- filter(bulkExperimentMetadata,
                                 !Celltype %in% celltypesToExclude,
                                 SampleID %in% names(bulkTPMData))
row.names(bulkExperimentMetadata) <- bulkExperimentMetadata$SampleID
bulkTPMData <- dplyr::select(bulkTPMData,
                             match(bulkExperimentMetadata$SampleID, names(bulkTPMData)))
logBulkTPM <- log2(bulkTPMData + 1)

assert_that(all(bulkExperimentMetadata$SampleID == colnames(bulkTPMData)))

# Remove all genes that don't have logTPM of at least 4 at some point
maxTPM <- apply(logBulkTPM, 1, max)
minTPM <- apply(logBulkTPM, 1, min)
logBulkTPM <- logBulkTPM[(maxTPM >= 2) & ((maxTPM - minTPM) >= 2), ]

# Batch correction
TFExprSet <- ExpressionSet(
  assayData = as.matrix(logBulkTPM),
  phenoData = new(
    "AnnotatedDataFrame", 
    data = bulkExperimentMetadata))

phenoData <- pData(TFExprSet)
batch <- phenoData$Batch
modcombat <- model.matrix(~as.factor(Celltype), data = phenoData)
combatExpr <- ComBat(
  dat = exprs(TFExprSet), 
  batch = batch, 
  mod = modcombat, 
  par.prior = FALSE, 
  prior.plots = FALSE)
combatExpr[combatExpr < 0] <- 0

logBulkTPM <- combatExpr

### Get scRNA data
# Parameters
minTPM <- 10
minCellsWithMinTPM <- 20

# Read in data
scTPMData <- read.table('sc_tpm.txt',
                        header = TRUE)
scExperimentMetadata <- getScExperimentMetadata()

# Make sure the genes line up
assert_that(all(rownames(scTPMData) == rownames(bulkTPMData)))

# Filter out all genes that don't have at least minTPM TPM in at least minCellsWithMinTPM cells
idxToKeep <- which(rowSums(scTPMData >= minTPM) >= minCellsWithMinTPM)
scTPMData <- scTPMData[idxToKeep, ]
logScTPM <- log2(scTPMData + 1)

### Combine bulk and scRNA data
# Filter out all genes that do not appear in both logScTPM and logBulkTPM
idxToKeep <- which(rownames(logScTPM) %in% rownames(logBulkTPM))
logScTPM <- logScTPM[idxToKeep, ]
idxToKeep <- which(rownames(logBulkTPM) %in% rownames(logScTPM))
logBulkTPM <- logBulkTPM[idxToKeep, ]

# Plot average expression in single cells vs. average expression in bulk
for (celltype in c('D0 H7 hESC', 'D1 APS', 'D1 MPS', 'D2 DLL1+ PXM', 'D3 Somite', 'D5 Dermomyotome', 'D6 Sclerotome')){

  scIdxOfCelltype = (scExperimentMetadata$Celltype == celltype) 
  logScTPMCelltype = logScTPM[, scIdxOfCelltype]
  bulkIdxOfCelltype = (bulkExperimentMetadata$Celltype == celltype) 
  logBulkTPMCelltype = logBulkTPM[, bulkIdxOfCelltype]
  
  statsToPlot = data.frame(
    bulkMean = rowMeans(logBulkTPMCelltype),
    scMean = rowMeans(logScTPMCelltype))
  
  ggplot(statsToPlot, aes_string(x = "bulkMean", y = "scMean")) +
    geom_point(size = 1, alpha = 0.40) +
    theme(text = element_text(size=20)) +   
    xlab("Average expression in bulk") +
    ylab("Average expression in single cells") 
  
  print(paste(celltype, cor(statsToPlot)[1,2]))
}

