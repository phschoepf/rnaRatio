library(tibble)
library(tidyr)
library(cBioPortalData)
library(AnVIL)

################
# input vars  #
################
selStudyList = c("gbm_tcga_pan_can_atlas_2018", "chol_nccs_2013") # a vector of studyIds to investigate (pull e.g. from manually imported data)
selGenes = c("MYC", "BASP1", "PHB", "MTOR") # a vector of hugoGeneSymbols

cbio <- cBioPortal()

# Generates a translation table for Entrez to Hugo gene symbols.
# Making this table from scratch takes a long time, therefore the function caches its output in the PWD.
generateTranslatorTable <- function(inputGeneList) {
  # Check cache existence, and if yes, consistency
  if(file.exists("translatorTable.csv")) {
    translatorTable <- read.csv("translatorTable.csv")
    if(all(inputGeneList %in% translatorTable$hugoGeneSymbol)) {
      return(translatorTable)
    }
  }
  
  # Generate from scratch if cache miss
  geneTablePage <- 1
  translatorTable <- tibble()
  while (dim(translatorTable)[[1]] != length(inputGeneList)) {
    geneTable <- geneTable(cbio, pageNumber = geneTablePage)
    filteredGeneTable <-
      geneTable[geneTable$hugoGeneSymbol %in% inputGeneList, ]
    translatorTable <- rbind(translatorTable, filteredGeneTable)
    geneTablePage <- geneTablePage + 1
    if (geneTablePage > 170) {
      stop("Invalid gene symbols")
    }
  }
  
  #cache result in PWD
  write.csv(translatorTable, file = "translatorTable.csv", row.names = FALSE)
  return(translatorTable)
}


getMolecularData <- function(study) {
  # make sampleId list
  samples <- allSamples(cbio, studyId = study)
  
  # get mRNA expression molecular profiles, or empty tibble if there are no data
  rnaProfiles <-
    molecularProfiles(cbio, studyId = study, projection = "ID") %>%
    filter(grepl("_rna_seq_v2_mrna$", molecularProfileId))
  
  #return if no RNA expression profile exists
  if (!isEmpty(rnaProfiles$molecularProfileId)) {
    return(
      molecularData(
        api = cbio,
        molecularProfileIds = rnaProfiles$molecularProfileId,
        entrezGeneIds = translatorTable$entrezGeneId,
        sampleIds = samples$sampleId
      )
    )
  }
  else {
    return(NULL)
  }
}

translate <- function(input) {
  if(class(input) == "integer") { # translate from Entrez to Hugo
    return(translatorTable[translatorTable$entrezGeneId == input, ]$hugoGeneSymbol)
  }
  else if(class(input) == "character") { # translate from Hugo to Entrez
    return(translatorTable[translatorTable$hugoGeneSymbol == input, ]$entrezGeneId)
  }
  else {
    return(NULL)
  }
}

#################
# main function #
#################

translatorTable <- generateTranslatorTable(selGenes)

# get data from all selected studies
all.data <- tibble()
for (study in selStudyList) {
  new.data <- getMolecularData(study)
  all.data <- rbind(all.data, new.data[[1]])
}
# restructure output table
wideOutput <- all.data %>% 
  mutate(all.data, hugoGeneSymbol = sapply(entrezGeneId, translate)) %>%
  select(sampleId, hugoGeneSymbol, value) %>% 
  pivot_wider(names_from = hugoGeneSymbol, values_from = value)