# Main script to call all other functions.
# Author: Philemon Schoepf <philemon.schoepf@student.ubik.ac.at>
# Date: 2021-03-22

# Global Variables -----------------------------------------------------------

INPUT_PATH = "../input_data/"
OUTPUT_PATH = "../images/"
OUTPUT_FORMAT = "png"

#image dims in mm (300 dpi)
WIDTH = 100
HEIGHT = 56

DO_WRITE <- F # if T, write output to file, otherwise show it on GUI

# Input Variables (for testing) ---------------------------------------------------------

patientUrls <-  read.csv(paste0(INPUT_PATH, "source_pat.csv")) # list of bookmark URLs from cBioPortal
patientFilter <- c("Bowel","Myeloma", "Lymphoma", "Cervix", "Breast", "Glioma", "Melanoma")
cellFilter <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LARGE_INTESTINE",
                "CERVIX", "BREAST", "SKIN", "CENTRAL_NERVOUS_SYSTEM")
selectedCellFilter <- c("K562", "MOLT4", "SW480", "MCF7", "LN18", "HEK293T", "HL60", "HeLa",
                        "HT29", "CACO2", "A375", "MDAMB231", "SUM159", "HCC159")


# Run ---------------------------------------------------------------------

#' RNARatioR
#'
#' @description The main function of rnaRatioR. You should not need to call
#' any other functions for normal use.
#'
#' @param gene1 The "numerator" gene.
#' @param gene2 The "denominator" gene,
#' @param patientInput A 2-column table with colnames =c("origin", "studies")
#' @param patientFilter A character vector of cancer types to include for patients.
#' @param cellFilter A character vector of cancer types to include for cells.
#' @param selectedCellFilter A character vector of cell lines to specially annotate..
#' @param forceReload Force reload of the cBioPortalData.
#'
#' @return Either outputs images or plot to the GUI.
#'
#' @examples rnaRatioR("MYC", "BASP1", patientUrls, patientFilter, cellFilter, selectedCellFilter)
#' @export
#'
rnaRatioR <- function(gene1, gene2, patientInput, patientFilter = "", cellFilter = "", selectedCellFilter = "", forceReload = TRUE){

  source("R/fetchCbioData.r")
  source("R/makePlots.r")

  if(!exists("filteredPatientData") | !exists("filteredCellData") | forceReload == T) {
    fetchCbioData(gene1, gene2)
  }
  plotRatio(gene1, gene2)
}

