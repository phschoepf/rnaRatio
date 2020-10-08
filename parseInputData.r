# Common input data parser for all plotting scripts.
# Author: Philemon Schöpf philemon.schoepf@student.uibk.ac.at

INPUT_DATA_DIR = "input_data/"

inputCellLines <-
  "cell lines - MYC_BASP1_PHB - mRNA expression (RNA-Seq RPKM).txt"

inputPatients <- dir(path = INPUT_DATA_DIR, pattern = "*RSEM*")

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

#get subset of selected cell lines
selectedCellLines <- cellLines %>% 
  filter(grepl("K562|MOLT4|SW480|MCF7|LN18|HEK293T|HL60|HeLa|HT29|CACO2|A375", SAMPLE_ID, ignore.case = T))%>% 
  mutate(Name = gsub("_.*", "", SAMPLE_ID, ignore.case = T))

