library(tidyverse)
library(cBioPortalData)
library(AnVIL)

################
# input vars  #
################
selStudyList = c("gbm_tcga_pan_can_atlas_2018", "chol_nccs_2013") # a vector of studyIds to investigate (pull e.g. from manually imported data)
selGenes = c("MYC", "BASP1", "PHB") # a vector of hugoGeneSymbols

cbio <- cBioPortal()

# translates hugoGeneSymbols to entrezGeneIds
# this function takes LONG, therefore cache when possible
translateGeneSymbols <- function(hugoGeneSymbols) {
  geneTablePage <- 1
  translatorTable <- tibble()
  while (dim(translatorTable)[[1]] != length(selGenes)) {
    geneTable <- geneTable(cbio, pageNumber = geneTablePage)
    filteredGeneTable <-
      geneTable[geneTable$hugoGeneSymbol %in% selGenes, ]
    translatorTable <- rbind(translatorTable, filteredGeneTable)
    geneTablePage <- geneTablePage + 1
    if (geneTablePage > 170) {
      stop("Invalid gene symbols")
    }
  }
  #cache for faster reloading, TODO: check for updates in input data
  write.csv(translatorTable, file = "translatorTable.csv", row.names = FALSE)
  return(translatorTable)
}

getMolecularData <- function(study) {
  # make sampleId list
  samples <- allSamples(cbio, studyId = study)
  
  # make molecularProfileId (mRNA)
  molecularProfile <- paste0(study, "_rna_seq_v2_mrna")
  
  # fetch data
  mols <- molecularData(
    api = cbio,
    molecularProfileIds = molecularProfile,
    entrezGeneIds = translatorTable$entrezGeneId,
    sampleIds = samples$sampleId
  )
  return(mols)
}

#################
# main function #
#################

#TODO try-catch
translatorTable <- read.csv("translatorTable.csv")
if(!exists("translatorTable")) {
  translatorTable <- translateGeneSymbols(selGenes)
}

# get data from all selected studies
all.data <- tibble()
for (study in selStudyList) {
  new.data <- getMolecularData(study)
  all.data <- rbind(all.data, new.data[[1]])
}

# restructure output table

