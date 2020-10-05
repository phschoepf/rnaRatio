library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)

INPUT_DATA_DIR = "input_data/"


inputPatients <- c("lymphoid_patients - MYC_BASP1_PHB - mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt",
                   "lymphoid_patients - MYC_BASP1_PHB - mRNA expression (RNA-Seq RPKM).txt",
                   "myeloid_patients - MYC_BASP1_PHB - mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt",
                   "myeloid_patients - MYC_BASP1_PHB - mRNA expression (RNA-Seq RPKM).txt")

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
  dataframe["Origin"] = filename #add row with patient cancer origin
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

# plotting patients
patplot <-
  ggplot(allPatientData, aes(
    y = Origin,
    x = log10(MYC / BASP1),
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(-2,4))) +
  labs(title = "Patient samples",
       y = "" ,
       x = expression(paste(log[10](MYC / BASP1), ", RNAseq v2 RSEM"))) +
  theme(axis.text.y = element_text(size = 20)) +
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
png("images/rpkm_rsem_demo.png", width = 2400, height = 1100)
plot(patplot)
dev.off()