# Main script to call all other functions.
# Author: Philemon Schoepf <philemon.schoepf@student.ubik.ac.at>
# Date: 2021-03-22

# Imports -----------------------------------------------------------------
#' @import cBioPortalData
NULL

# Global Variables -----------------------------------------------------------

#image dims in mm (300 dpi)
WIDTH = 100
HEIGHT = 56

# Main function ---------------------------------------------------------------------

#' rnaRatio
#'
#' @description The main function of rnaRatio. You should not need to call
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
#' @examples rnaRatio("MYC", "BASP1", patientUrls, patientFilter, cellFilter, selectedCellFilter)
#'
#' @export
rnaRatio <- function(gene1, gene2, patientInput, patientFilter = "", cellFilter = "", selectedCellFilter = "", forceReload = FALSE){

  if(!exists("filteredPatientData") | !exists("filteredCellData") | forceReload == T) {
    fetchCbioData(gene1, gene2, ...)
  }
  plotRatio(gene1, gene2)
}
