library(ggplot2)
library(dplyr)
library(assertthat)

getGeneIdWithPrefix <- function(geneData, targetGeneIdPrefix){
  return(as.character(with(geneData, geneid[which(startsWith(geneid, targetGeneIdPrefix))])))
}

getGeneId <- function(geneData, targetGeneSymbol){
  return(as.character(with(geneData, geneid[which(geneSymbol == targetGeneSymbol)])))
}

getGeneSymbol <- function(geneData, targetGeneId){
  return(with(geneData, geneSymbol[which(geneid == targetGeneId)]))
}


# PCA
# expr is a matrix of dimension # samples (rows) by # genes (columns)
PCA <- function(expr, experimentMetadata, PCs = c(1,2), intgroup) {
  
  assert_that(length(PCs) == 2)
  assert_that(length(intgroup) <= 2)
  pca <- prcomp(expr)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  intgroup.df <- experimentMetadata[intgroup]
  
  d <- data.frame(PCy = pca$x[, PCs[1]], PCx = pca$x[, PCs[2]], 
                  intgroup.df, name = rownames(experimentMetadata))
  percentVar <- round(100 * c(percentVar[PCs[1]], percentVar[PCs[2]]))
  
  if (length(intgroup) == 1) {
    p <- ggplot(d, aes_string("PCx", "PCy", color = intgroup[1])) +
      geom_point(size = 6, alpha = 0.60) +
      theme(text = element_text(size=20)) + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +      
      xlab(sprintf("PC%s: %s%% variance ", PCs[2], percentVar[2])) +
      ylab(sprintf("PC%s: %s%% variance ", PCs[1], percentVar[1]))
  } else if (length(intgroup) == 2) {
    p <- ggplot(d, aes_string("PCx", "PCy", color = intgroup[1], shape = intgroup[2])) +
      geom_point(size = 6, alpha = 0.60) +
      theme(text = element_text(size=20)) + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +      
      xlab(sprintf("PC%s: %s%% variance ", PCs[2], percentVar[2])) +
      ylab(sprintf("PC%s: %s%% variance ", PCs[1], percentVar[1]))
  }
  print(p)
  return(pca)
}

