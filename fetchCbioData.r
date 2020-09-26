library(ggplot2)
library(gridExtra)
library(tibble)
library(cBioPortalData)
library(AnVIL)

#pull data from API and get gene expression information
cbio <- cBioPortal()

all.studies <- getStudies(cbio)
filtered.studies <- all.studies %>% filter(grepl(".*myelo.*|.*lympho.*", name, ignore.case = TRUE))
filtered.studies <- filtered.studies[- grepl("cell line", filtered.studies$name, ignore.case = TRUE), ]

sampleList <- list()
expList <- list()
i <- 1
for(studyId in filtered.studies$studyId) {
  sampleList[[i]] <- allSamples(cbio, studyId)
  expList[[i]] <- molecularData(api = cbio,
                molecularProfileIds = c("acc_tcga_rna_seq_v2_mrna"),
                entrezGeneIds = c(4609, 10409),
                sampleIds = sampleList[[i]]$sampleId)
  i <- i + 1
}

concatSampleList <- do.call(rbind, sampleList)
sampleIdList <- paste0(samples$patientId, '-01')
testList <- c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")

expressionData <- molecularData(api = cbio,
                      molecularProfileIds = c("acc_tcga_rna_seq_v2_mrna"),
                      entrezGeneIds = c(4609, 10409),
                      sampleIds = sampleIdList
)