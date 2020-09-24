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

# Data source
# https://www.cbioportal.org/results/oncoprint?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=ccle_broad_2019%2Ccellline_ccle_broad%2Ccellline_nci60&case_set_id=all&data_priority=0&gene_list=MYC%250ABASP1&geneset_list=%20&plots_coloring_selection=%7B%22selectedOption%22%3A%224609_undefined%22%7D&plots_horz_selection=%7B%22dataType%22%3A%22MRNA_EXPRESSION%22%2C%22selectedDataSourceOption%22%3A%22rna_seq_mrna%22%7D&plots_vert_selection=%7B%22selectedGeneOption%22%3A10409%2C%22selectedDataSourceOption%22%3A%22rna_seq_mrna%22%7D&profileFilter=0&tab_index=tab_visualize
# https://www.cbioportal.org/results/cancerTypesSummary?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=all_stjude_2015%2Call_stjude_2016%2Clcll_broad_2013%2Ccll_broad_2015%2Ccll_iuopa_2015%2Ccllsll_icgc_2011%2Cctcl_columbia_2015%2Cdlbcl_dfci_2018%2Cdlbc_broad_2012%2Cdlbcl_duke_2017%2Cdlbc_tcga_pan_can_atlas_2018%2Cnhl_bcgsc_2013%2Cmcl_idibips_2013%2Cmbn_mdacc_2013%2Cmm_broad%2Cnhl_bcgsc_2011%2Call_phase2_target_2018_pub%2Cpcnsl_mayo_2015%2Caml_ohsu_2018%2Claml_tcga_pan_can_atlas_2018%2Cmnm_washu_2016%2Chistiocytosis_cobi_msk_2019%2Call_stjude_2013%2Cmds_tokyo_2011%2Cmds_mskcc_2020%2Cmpn_cimr_2013%2Caml_target_2018_pub&case_set_id=all&data_priority=0&gene_list=MYC%2520BASP1&geneset_list=%20&profileFilter=0&tab_index=tab_visualize
cellLines<- read.delim("input_data/cell lines - MYC_BASP1 - mRNA expression (RNA-Seq RPKM).txt")
lymphoidPatients <- read.delim("input_data/lymphoid_patients - MYC_BASP1 -  (RNA-Seq RPKM).txt")
myeloidPatients <- read.delim("input_data/myeloid_patients - MYC_BASP1 -  (RNA-Seq RPKM).txt")

#filter NaN and zero values
attach(cellLines)
cellLines <- cellLines[which(MYC != "NA" & MYC != 0 & BASP1 != "NA" & BASP1 != 0), ]
detach(cellLines)

attach(myeloidPatients)
myeloidPatients <- myeloidPatients[which(MYC != "NA" & MYC != 0 & BASP1 != "NA" & BASP1 != 0), ]
myeloidPatients["Origin"] <- "myeloid"
detach(myeloidPatients)


attach(lymphoidPatients)
lymphoidPatients <- lymphoidPatients[which(MYC != "NA" & MYC != 0 & BASP1 != "NA" & BASP1 != 0), ]
lymphoidPatients["Origin"] <- "lymphoid"
detach(lymphoidPatients)

combinedTablePatients <- add_row(myeloidPatients, lymphoidPatients)

#extract sites of origin and combine
siteOfOriginCells <- gsub("^.*?_", "", cellLines$SAMPLE_ID)
combinedTableCells <- add_column(cellLines, siteOfOriginCells, .after = "SAMPLE_ID")

#expression plot (deprecated, doesn't show data nicely)
#expressions <- ggplot(combinedTableCells, aes(x = MYC, y = BASP1)) +
#  geom_point(aes(col = siteOfOriginCells)) +
#  labs(title = "BASP1 vs. MYC expression levels by primary cancer type", x = "MYC [RNAseq PRKM]", y = "BASP1 [RNAsqe RPKM]")
#plot(expressions)

#plot cell lines by cancer origin
ratiosCells <- ggplot(combinedTableCells, aes(y = siteOfOriginCells, x = log10(MYC/BASP1))) +
  #geom_point() +
  #xlim(0, 7000) +
  labs(title = "MYC to BASP1 in cell lines", y = "", x= "MYC/BASP1 [log10 RNAseq RPKM ratio]") +
  geom_boxplot(outlier.color = "red", outlier.shape = 8) +
  geom_jitter(shape=16, position=position_jitter(w=0, h=0.15))

#plot patient data by cancer origin
ratiosPatients <- ggplot(combinedTablePatients, aes(y = Origin, x = log10(MYC/BASP1))) +
  #geom_point() +
  #xlim(0, 7000) +
  labs(title = "MYC to BASP1 in patient samples", y = "", x = "MYC/BASP1 [log10 RNAseq RPKM ratio]") +
  geom_boxplot(outlier.color = "red", outlier.shape = 8) +
  geom_jitter(shape=16, position=position_jitter(w=0, h=0.2))

grid.arrange(plot(ratiosCells), plot(ratiosPatients))