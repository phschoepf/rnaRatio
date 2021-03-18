# Parses data downloaded by cBioPortalData into plottable table.
# @author Philemon Schöpf <philemon.schoepf@student.ubik.ac.at>

library(tidyverse)
source("fetchCbioData.r")

#################################
#    input variables            #
#################################
DO_FILTER = TRUE # filter cell lines for selected ones?
INPUT_DATA_DIR = "input_data/"
patientUrls <-  read.csv(paste0(INPUT_DATA_DIR, "source_pat.csv")) # list of bookmark URLs from cBioPortal
selGenes <- c("MYC", "BASP1", "PHB") # genes of interest

#################################
#################################


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

