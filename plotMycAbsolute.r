library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)

INPUT_DATA_DIR <- "input_data/"

inputPatients <- dir(path = INPUT_DATA_DIR, pattern = "*RSEM*")

# parse and filter patient file
parseCsv <- function(file) {
  csv <- read.delim(paste0(INPUT_DATA_DIR, file))
  csv[which(csv$BASP1 != "NA" &
              csv$BASP1 != 0 &
              csv$MYC != "NA" & csv$MYC != 0),] #filter invalid values
}

allPatientData <- data.frame()
#annotate patient data w/ origin and combine
for (filename in inputPatients) {
  dataframe <- parseCsv(filename)
  dataframe["Origin"] = strsplit(filename, "_.*")[[1]] #add row with patient cancer origin
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
    x = reorder(Origin, -MYC, FUN = "median"),
    y = MYC,
    na.rm = TRUE
  )) +
  ylim(0,15000) +
  labs(title = "Patient samples",
       x = "" ,
       y = "MYC RNA expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"
  ) +
  stat_boxplot(geom = "errorbar", width = 0.6, lwd = 1) +
  geom_boxplot(outlier.shape = NA,
               lwd = 1,
               fatten = 1.2) +
  geom_jitter(
    shape = 16,
    position = position_jitter(w = 0.25, h = 0),
    size = 0.8,
    color = "#f39200"
  )

#export to png
png("images/MYC_pat_absolute.png", width = 2300, height = 1000)
grid.arrange(patplot)
dev.off()
