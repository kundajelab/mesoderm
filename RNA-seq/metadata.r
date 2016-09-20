library(dplyr)
library(RCurl)
library(stringr)

getCSVFromURL <- function(url){
  read.csv(text = getURL(url),
           colClasses = c('character'))
}

getScExperimentMetadata <- function(){
  experimentNames <- read.csv('sc_tpm.csv', header = TRUE, nrows = 1)
  experimentMetadata <- data.frame(
    SampleID = colnames(experimentNames),
    stringsAsFactors = FALSE)
  experimentMetadata<- mutate(
    experimentMetadata,
    Celltype = str_split_fixed(SampleID, "[.]", 2)[, 1])
   
  experimentMetadata <- experimentMetadata[3:nrow(experimentMetadata), ] # Remove geneID and geneSymbol
  
  experimentMetadata[experimentMetadata$Celltype == 'H7hESC', 'Celltype'] <- 'D0 H7 hESC'
  experimentMetadata[experimentMetadata$Celltype == 'APS', 'Celltype'] <- 'D1 APS'
  experimentMetadata[experimentMetadata$Celltype == 'MPS3', 'Celltype'] <- 'D1 MPS'
  experimentMetadata[experimentMetadata$Celltype == 'DLL1PXM', 'Celltype'] <- 'D2 DLL1+ PXM'
  experimentMetadata[experimentMetadata$Celltype == 'Earlysomite', 'Celltype'] <- 'D3 Somite'
  experimentMetadata[experimentMetadata$Celltype == 'Sclerotome', 'Celltype'] <- 'D6 Sclerotome'  
  experimentMetadata[experimentMetadata$Celltype == 'cDM', 'Celltype'] <- 'D5 Dermomyotome'
  experimentMetadata[experimentMetadata$Celltype == 'LatM', 'Celltype'] <- 'D2 LatM'
  
  return(experimentMetadata)
}

getExperimentMetadata <- function(){
  experimentMetadata <- getCSVFromURL('https://docs.google.com/spreadsheets/d/1k1gDoEEeT3v2Ot2GEli1aKFx6WdcffuAK8SuHGFqN3Q/pub?gid=442056520&single=true&output=csv')

  experimentMetadata <- mutate(
    experimentMetadata,
    label = paste(SampleID, Celltype, Batch))

  return(experimentMetadata)
}
  
getCelltypesToExclude <- function(){
  # In our original experiments we also profiled a different line of hESCs,
  # which we ended up not using in our final paper for consistency reasons.
  # Here, we exclude them from the subsequent analysis.
  return(c('NG hESC', 'NG Cardiac'))  
}

getCelltypesToCompare <- function(){
  celltypesToCompare <- data.frame(
    matrix(
      c('D1 APS', 'D0 H7 hESC',
        'D1 MPS', 'D0 H7 hESC',
        'D1 APS', 'D1 MPS',
        'D2 DLL1pos PXM', 'D2 DLL1neg PXM',
        'D2 DLL1pos PXM', 'D0 H7 hESC',
        'D2 DLL1pos PXM', 'D1 APS',
        'D3 Somite', 'D2 DLL1pos PXM', 
        'D6 Sclerotome', 'D3 Somite',
        'D5 Dermomyotome', 'D3 Somite',
        'D6 Sclerotome', 'D5 Dermomyotome',
        'D2 LatM', 'D0 H7 hESC',
        'D2 LatM', 'D1 MPS',
        'D3 Cardiac', 'D2 LatM'
        ),
      ncol = 2,
      byrow = TRUE),
    stringsAsFactors = FALSE)
  colnames(celltypesToCompare) <- c('later', 'earlier')
  celltypesToCompare <- mutate(celltypesToCompare, 
                               both = paste(later, earlier, sep = "-"))
  
  return(celltypesToCompare)
}