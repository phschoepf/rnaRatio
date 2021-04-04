# Sets up input data for the MYC/BASP1 ratio plot.
# Author: Philemon Schoepf <philemon.schoepf@student.ubik.ac.at>
# Date: 2021-03-22

library(rnaRatio)

OUTPUT_PATH = "../images/mtor-demo/"
OUTPUT_FORMAT = "png"
DO_WRITE = T # if T, write output to file, otherwise show it on GUI

# Input Variables ---------------------------------------------------------

patientUrls <-  read.csv("../input_data/source_pat.csv") # list of bookmark URLs from cBioPortal
patientFilter <- c("Ovary","Breast", "Cervix", "Lung", "Melanoma", "Endometrial cancer")
cellFilter <- c("OVARY", "BREAST", "CERVIX", "LUNG", "SKIN", "ENDOMETRIUM")
selectedCellFilter <- c("Hela", "MCF7", "EFM19", "LC1F")


# Run ---------------------------------------------------------------------

rnaRatio("MTOR", "IRS1", patientUrls, patientFilter, cellFilter, selectedCellFilter)
