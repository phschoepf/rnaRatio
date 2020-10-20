# Method script for downloading data from the cBioPortal API.
# @author Philemon Schöpf <philemon.schoepf@student.ubik.ac.at>

library(tibble)
library(tidyr)
library(cBioPortalData)
library(AnVIL)

cbio <- cBioPortal()

getStudiesFromUrl <- function(url) {
  split <- url %>%
    gsub(pattern = "%25", replacement = "%") %>%
    gsub(pattern = "%2C", replacement = ",") %>%
    gsub(pattern = "%20", replacement = " ") %>%
    strsplit("&")
  studyList <- split[[1]][4] %>%
    gsub(pattern = "cancer_study_list=", replacement = "") %>%
    strsplit(",")
  return(studyList[[1]])
  
}

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

getMolecularData <- function(study, molecularProfile, translatorTable) {
  # make sampleId list
  samples <- allSamples(cbio, studyId = study)
  
  # get mRNA expression molecular profiles, or empty tibble if there are no data
  rnaProfiles <-
    molecularProfiles(cbio, studyId = study, projection = "ID") %>%
    filter(grepl(paste0("_", molecularProfile, "$"), molecularProfileId))
  
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

translate <- function(input, table) {
  if(class(input) == "integer") { # translate from Entrez to Hugo
    return(table[table$entrezGeneId == input, ]$hugoGeneSymbol)
  }
  else if(class(input) == "character") { # translate from Hugo to Entrez
    return(table[table$hugoGeneSymbol == input, ]$entrezGeneId)
  }
  else {
    return(NULL)
  }
}

#################
# main function #
#################
fetchPat <- function(url, selGenes) {
  translatorTable <- generateTranslatorTable(selGenes)
  selStudyList <- getStudiesFromUrl(url)
  
  # get data from all selected studies
  all.data <- tibble()
  for (study in selStudyList) {
    new.data <- getMolecularData(study, molecularProfile = "rna_seq_v2_mrna", translatorTable =  translatorTable)
    all.data <- rbind(all.data, new.data[[1]])
  }
  # restructure output table
  wideOutput <- all.data %>%
    mutate(all.data, hugoGeneSymbol = sapply(X = entrezGeneId, FUN = translate, translatorTable)) %>%
    select(studyId, sampleId, hugoGeneSymbol, value) %>%
    pivot_wider(names_from = hugoGeneSymbol, values_from = value)
  return(wideOutput)
}

fetchCells <- function(url, selGenes) {
  translatorTable <- generateTranslatorTable(selGenes)
  selStudyList <- getStudiesFromUrl(url)
  
  # get data from all selected studies
  all.data <- tibble()
  for (study in selStudyList) {
    new.data <- getMolecularData(study, molecularProfile = "rna_seq_mrna", translatorTable = translatorTable)
    all.data <- rbind(all.data, new.data[[1]])
  }
  # restructure output table
  wideOutput <- all.data %>%
    mutate(all.data, hugoGeneSymbol = sapply(X = entrezGeneId, FUN = translate, translatorTable)) %>%
    select(studyId, sampleId, hugoGeneSymbol, value) %>%
    pivot_wider(names_from = hugoGeneSymbol, values_from = value)
  return(wideOutput)
}