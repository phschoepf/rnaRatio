library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)

## Data sources
# Cell lines: 
# https://www.cbioportal.org/results/download?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=ccle_broad_2019%2Ccellline_ccle_broad%2Ccellline_nci60&case_set_id=all&data_priority=0&gene_list=MYC%2520BASP1&geneset_list=%20&profileFilter=0&tab_index=tab_visualize

# Lymphoid patients:
# https://www.cbioportal.org/results/download?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=all_stjude_2015%2Call_stjude_2016%2Clcll_broad_2013%2Ccll_broad_2015%2Ccll_iuopa_2015%2Ccllsll_icgc_2011%2Cctcl_columbia_2015%2Cdlbcl_dfci_2018%2Cdlbc_broad_2012%2Cdlbcl_duke_2017%2Cdlbc_tcga_pan_can_atlas_2018%2Cnhl_bcgsc_2013%2Cmcl_idibips_2013%2Cmbn_mdacc_2013%2Cmm_broad%2Cnhl_bcgsc_2011%2Call_phase2_target_2018_pub%2Cpcnsl_mayo_2015&case_set_id=all&data_priority=0&gene_list=MYC%2520BASP1&geneset_list=%20&profileFilter=0&tab_index=tab_visualize

# Myeloid patients:
# https://www.cbioportal.org/results/download?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=aml_ohsu_2018%2Claml_tcga_pan_can_atlas_2018%2Cmnm_washu_2016%2Chistiocytosis_cobi_msk_2019%2Call_stjude_2013%2Cmds_tokyo_2011%2Cmds_mskcc_2020%2Cmpn_cimr_2013%2Caml_target_2018_pub&case_set_id=all&data_priority=0&gene_list=MYC%2520BASP1&geneset_list=%20&profileFilter=0&tab_index=tab_visualize

# Bowel patients:
# https://www.cbioportal.org/results/download?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=coad_caseccc_2015%2Ccoad_cptac_2019%2Ccoadread_dfci_2016%2Ccoadread_genentech%2Ccoadread_tcga_pan_can_atlas_2018%2Ccoadread_mskcc%2Ccrc_msk_2017%2Crectal_msk_2019&case_set_id=all&data_priority=0&gene_list=MYC%2520BASP1&geneset_list=%20&profileFilter=0&tab_index=tab_visualize

INPUT_DATA_DIR = "input_data/"
inputCellLines <-
  "cell lines - MYC_BASP1 - mRNA expression (RNA-Seq RPKM).txt"
inputPatients <- c(
  "lymphoid_patients - MYC_BASP1 - mRNA expression (RNA-Seq RPKM).txt",
  "myeloid_patients - MYC_BASP1 - mRNA expression (RNA-Seq RPKM).txt",
  "bowel_patients - MYC_BASP1 - mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt"
)

# parse and filter patient file
parseCsv <- function(file) {
  csv <- read.delim(paste0(INPUT_DATA_DIR, file))
  csv[which(csv$MYC != "NA" &
              csv$MYC != 0 &
              csv$BASP1 != "NA" & csv$BASP1 != 0),] #filter invalid values
}

allPatientData <- data.frame()
#annotate patient data w/ origin and combine
for (filename in inputPatients) {
  dataframe <- parseCsv(filename)
  dataframe["Origin"] = strsplit(filename, " - .*")[[1]] #add row with patient cancer origin
  allPatientData <-
    rbind(allPatientData, dataframe) #add data to global list
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

#annotate cell data w/ origin and merge to patient dataframe
cellLines <-  parseCsv(inputCellLines)
Origin <- gsub("^.*?_", "", cellLines$SAMPLE_ID)
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

# plotting cells
cellplot <-
  ggplot(cellLines, aes(
    y = reorder(Origin, MYC / BASP1, FUN = median),
    x = log10(MYC / BASP1),
    na.rm = TRUE
  )) +
  xlim(-4, 4) +
  labs(title = "Cell lines", y = "", x = "") +
  stat_boxplot(geom = "errorbar", width = 0.6, lwd = 1) +
  geom_boxplot(outlier.shape = NA,
               lwd = 1,
               fatten = 1.2) +
  geom_jitter(
    shape = 16,
    position = position_jitter(w = 0, h = 0.15),
    size = 0.9,
    color = "#f39200"
  )

# plotting patients
patplot <-
  ggplot(allPatientData, aes(
    y = reorder(Origin, MYC / BASP1, FUN = median),
    x = log10(MYC / BASP1),
    na.rm = TRUE
  )) +
  xlim(-4, 4) +
  labs(title = "Patient samples",
       y = "" ,
       x = expression(log[10](MYC / BASP1) ^ {
         "[1]"
       })) +
  stat_boxplot(geom = "errorbar", width = 0.6, lwd = 1) +
  geom_boxplot(outlier.shape = NA,
               lwd = 1,
               fatten = 1.2) +
  geom_jitter(
    shape = 16,
    position = position_jitter(w = 0, h = 0.15),
    size = 0.8,
    color = "#f39200"
  )

#display plot
aligned <- align_plots(cellplot, patplot, align = "v", axis = "lr")
grid.arrange(aligned[[1]], aligned[[2]], layout_matrix = rbind(1, 1, 1, 1, 2))
