# Method script for downloading data from the cBioPortal API.
# @author Philemon Schöpf <philemon.schoepf@student.ubik.ac.at>

library(tibble)
library(cBioPortalData)
library(AnVIL)
library(org.Hs.eg.db)
library(dplyr)
library(httr)
library(stringr)

#################################
#    input variables            #
#################################
DO_FILTER = TRUE # filter cell lines for selected ones?
INPUT_DATA_DIR = "input_data/"
patientUrls <-  read.csv(paste0(INPUT_DATA_DIR, "source_pat.csv")) # list of bookmark URLs from cBioPortal
selGenes <- c("MYC", "BASP1", "PHB") # genes of interest

#################################
#################################

# Study IDs of cell line studies in cBioPortal.
# Should not change often, therefore hardcoded here.
CELL_LINE_STUDIES = c("cellline_nci60",
                      "cellline_ccle_broad",
                      "ccle_broad_2019"
                      )

# Create a CBioPortal API client
cbio <- cBioPortal()


#' getStudiesFromUrl
#'
#' @description Extracts the studyIds from a cBioPortal URL. This enables use of the web GUI for study selection by exporting the URLs.
#' This is intended to be used with patient data only, for the cell lines the study IDs are provided as constants.
#'
#' @param url a URL generated by a cBioPortal query
#'
#' @return a vector of studyIds
#' @export
#'
#' @examples
getStudiesFromUrl <- function(url) {

  cancerStudyString <- httr::parse_url(url)$query$cancer_study_list
  cancerStudiesAsVector <- cancerStudyString %>% strsplit(",")
  return(cancerStudiesAsVector)

}

#' toEntrez
#'
#' @description Converts a HUGO gene symbol to Entrez using org.Hs.eg.db.
#' Returns identity if it is already a Entrez symbol.
#'
#' @param gene Either a HUGO gene symbol string or a Entrez ID
#'
#' @return A Entrez Id
#'
#' @examples
#' translate("MYC")
#' translate(4609)
#' translate("4609")
toEntrez <- function(gene) {

  # AnnotationDbi only handles strings
  gene <- toString(gene)

  if(!str_starts(gene, "[0-9]")) {
    gene <- AnnotationDbi::select(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "SYMBOL")$ENTREZID
  }

  # convert back to int
  return(strtoi(gene))
}


#' toHugo
#'
#' @description Converts a Entrez Id to HUGO using org.Hs.eg.db.
#' Returns identity if it is already a HUGO symbol.
#'
#' @param gene Either a HUGO gene symbol string or a Entrez ID
#'
#' @return A HUGO symbol.
#'
#' @examples
#' translate("MYC")
#' translate(4609)
#' translate("4609")
toHugo <- function(gene) {

  # AnnotationDbi only handles strings
  gene <- toString(gene)

  if(str_starts(gene, "[0-9]")) {
    gene <- AnnotationDbi::select(org.Hs.eg.db, keys = gene, column = "SYMBOL", keytype = "ENTREZID")$SYMBOL
  }

  return(gene)
}


#' getMolecularData
#'
#' @description gets RNA expression profiles from cBioPortal for the given study and genes
#'
#' @param study a studyId string, e.g. "ccle_broad_2019"
#' @param molecularProfile a molecularProfile stub to search for. Usually "rna_seq_v2_mrna" or similar
#'
#' @return
#' @export
#'
#' @examples
getMolecularData <- function(study, molecularProfile) {
  # list all sampleIds from the study
  samples <- cBioPortalData::allSamples(cbio, studyId = study)

  # get mRNA expression molecular profiles, or empty tibble if there are no data
  rnaProfiles <-
    molecularProfiles(cbio, studyId = study, projection = "ID") %>%
    filter(str_detect(molecularProfileId, paste0(molecularProfile, "$")))

  #return if no RNA expression profile exists
  if (!isEmpty(rnaProfiles$molecularProfileId)) {
    return(
      molecularData(
        api = cbio,
        molecularProfileIds = rnaProfiles$molecularProfileId,
        entrezGeneIds = unlist(lapply(selGenes, toEntrez)),
        sampleIds = samples$sampleId
      )
    )
  }
  else {
    return(NULL)
  }
}

#' Make a table of RNA expression values.
#'
#' @param studyList List of studies
#' @param selGenes List of genes
#' @param molecularProfile The type of molecular profile to query. Default is "rna_seq_v2_mrna", for cell lines "rna_seq_mrna" has to be used.
#'
#' @return Table of RNA expression values for selGenes, one row per sample.
#' @export
#'
#' @examples
#' TODO
makeExpressionTable <- function(studyList, selGenes, molecularProfile = "rna_seq_v2_mrna") {

  # get data from all selected studies
  all.data <- accumulate(map(studyList, getMolecularData, molecularProfile = molecularProfile), rbind)
  #   tibble()
  # for (study in studyList) {
  #   new.data <- getMolecularData(study, molecularProfile = molecularProfile)
  #   all.data <- rbind(all.data, new.data[[1]])
  # }
  # # remove NA and restructure output
  # wideOutput <- all.data %>%
  #   mutate(all.data, hugoGeneSymbol = sapply(X = entrezGeneId, FUN = translate, translatorTable)) %>%
  #   select(studyId, sampleId, hugoGeneSymbol, value) %>%
  #   pivot_wider(names_from = hugoGeneSymbol, values_from = value)
  # return(wideOutput)
}



allPatientData <- data.frame()

#annotate patient data w/ origin and combine
lapply(patientUrls["url"], fetchPat, selGenes = selGenes)
for (row in 1:nrow(patientUrls)) {
  dataframe <- fetchPat(patientUrls[row, "url"], selGenes)
  dataframe["Origin"] = patientUrls[row, "Type"] #add row with patient cancer origin
  allPatientData <-
    rbind(allPatientData, dataframe)
}
for (origin in unique(allPatientData$Origin)) {
  numOfOccurances = sum(allPatientData$Origin == origin)
  allPatientData$Origin <-
    gsub(
      paste0("\\b", origin, "\\b"),
      paste0(origin, " (n = ", numOfOccurances, ")"),
      allPatientData$Origin
    )
}
allPatientData <- allPatientData %>% filter(BASP1 != 0)

#annotate cell data w/ origin and merge to patient dataframe
cellLines <- fetchCells(cellLineUrl[1, "url"], selGenes)
cellLines[which(cellLines$MYC != "NA" &
                  cellLines$MYC != 0 &
                  cellLines$BASP1 != "NA" & cellLines$BASP1 != 0),]
Origin <- gsub("^.*?_", "", cellLines$sampleId)
cellLines <- cbind(cellLines, Origin)
for (origin in unique(cellLines$Origin)) {
  numOfOccurances = sum(cellLines$Origin == origin)
  cellLines$Origin <-
    gsub(
      paste0("\\b", origin, "\\b"),
      paste0(origin, " (n = ", numOfOccurances, ")"),
      cellLines$Origin
    )
}

#subset to interesting cancer types, and rename them
allPatientData <- allPatientData %>%
  filter(grepl("Bowel|Myeloma|Lymphoma|Cervix|Breast|Glioma|Melanoma", Origin, ignore.case = T))
cellLines <- cellLines %>%
  filter(grepl("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|LARGE_INTESTINE|CERVIX|BREAST|SKIN|CENTRAL_NERVOUS_SYSTEM", Origin, ignore.case = T)) %>%
  mutate(Origin = gsub("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HAEMATOPOIETIC/LYMPHOID_TISSUE", Origin)) %>%
  mutate(Origin = gsub("_", " ", str_to_sentence(Origin)))

#get subset of selected cell lines
if(DO_FILTER) {
  selectedCellLines <- cellLines %>%
    filter(grepl("K562|MOLT4|SW480|MCF7|LN18|HEK293T|HL60|HeLa|HT29|CACO2|A375|MDAMB231|SUM159|HCC159", sampleId, ignore.case = T))%>%
    mutate(Name = gsub("_.*", "", sampleId, ignore.case = T))

}

