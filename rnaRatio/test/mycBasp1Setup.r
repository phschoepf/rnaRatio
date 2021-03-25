# Sets up input data for the MYC/BASP1 ratio plot.
# Author: Philemon Schoepf <philemon.schoepf@student.ubik.ac.at>
# Date: 2021-03-22

library(rnaRatio)

OUTPUT_PATH = "../images/"
OUTPUT_FORMAT = "png"
DO_WRITE = F # if T, write output to file, otherwise show it on GUI

# Input Variables ---------------------------------------------------------

patientUrls <-  read.csv("../input_data/source_pat.csv") # list of bookmark URLs from cBioPortal
patientFilter <- c("Bowel","Myeloma", "Lymphoma", "Cervix", "Breast", "Glioma", "Melanoma")
cellFilter <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LARGE_INTESTINE",
                "CERVIX", "BREAST", "SKIN", "CENTRAL_NERVOUS_SYSTEM")
selectedCellFilter <- c("K562", "MOLT4", "SW480", "MCF7", "LN18", "HEK293T", "HL60", "HeLa",
                        "HT29", "CACO2", "A375", "MDAMB231", "SUM159", "HCC159")


# Run ---------------------------------------------------------------------

rnaRatio("MYC", "BASP1", patientUrls, patientFilter, cellFilter, selectedCellFilter)
