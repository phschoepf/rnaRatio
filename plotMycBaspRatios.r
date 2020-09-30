library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)

INPUT_DATA_DIR = "input_data/"
inputCellLines <-
  "cell lines - MYC_BASP1_PHB - mRNA expression (RNA-Seq RPKM).txt"

inputPatients <- dir(path = INPUT_DATA_DIR, pattern = "*RSEM*")

#inputPatients <- c(
  #"lymphoid_patients - MYC_BASP1_PHB - mRNA expression (RNA-Seq RPKM).txt",
  #"myeloid_patients - MYC_BASP1_PHB - mRNA expression (RNA-Seq RPKM).txt",
  #"lymphoid_patients - MYC_BASP1_PHB - mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt",
  #"myeloid_patients - MYC_BASP1_PHB - mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt", 
  #"bowel_patients - MYC_BASP1_PHB - mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt"
#)

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

#export to png
aligned <- align_plots(cellplot, patplot, align = "v", axis = "lr")
png("images/myc_basp1_ratios_combined.png", width = 2300, height = 2000)
grid.arrange(aligned[[1]], 
             aligned[[2]],
             textGrob(expression(paste({}^{
         "[1]"
       }, "Cell lines: RNAseq RKPM, patient data: RNAseq RSEM")), x = 1, hjust = 1),
             layout_matrix = rbind(1, 1, 1, 1, 2, 2, 2, 2, 3))
dev.off()