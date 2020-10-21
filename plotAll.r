# Wrapper script that sources/calls methods from all other scripts.
# Set up what to plot here.
# @author Philemon Schöpf <philemon.schoepf@student.ubik.ac.at>

library(cowplot)
library(tidyverse)
source("plotMethods.r")

#generate input data anew if it doesn't exist
if(!exists("allPatientData") | !exists("cellLines")) {
  source("parseInputData.r")
}

mycbasp_cell <- plotCells(cellLines, "MYC", "BASP1", c(-2,4))
mycbasp_pat <- plotPatients(allPatientData, "MYC", "BASP1", c(-2,4))
#exportToPng(baseName = "myc_basp1_ratios", mycbasp_cell, mycbasp_pat)

mycphb_cell <- plotCells(cellLines, "MYC", "PHB", c(-3,2))
mycphb_pat <- plotPatients(allPatientData, "MYC", "PHB", c(-3,2))
#exportToPng(baseName = "myc_phb_ratios", mycphb_cell, mycphb_pat)

myc_cell <- plotCells(cellLines, "MYC", lims = c(NA,400))
myc_pat <- plotPatients(allPatientData, "MYC", lims =  c(NA,15000))
#exportToPng(baseName = "myc_absolute", myc_cell, myc_pat)

basp1_cell <- plotCells(cellLines, "BASP1", lims = c(NA,700))
basp1_pat <- plotPatients(allPatientData, "BASP1", lims =  c(NA,15000))
#exportToPng(baseName = "basp1_absolute", basp1_cell, basp1_pat)


#source("RPKMvsRSEMdemo.r") does not work currently w/ new APIrequest input system

#combined large plot
aligned_plots <- align_plots( mycbasp_cell, mycbasp_pat,
                              myc_cell, myc_pat,
                              basp1_cell, basp1_pat,
                              align = "h", axis = "tb"
)

ggarrange(plotlist = aligned_plots, nrow = 3, ncol = 2)
ggsave("allPlots.tiff", width = 170, height = 190, units = "mm")
plotCellsTable(selectedCellLines, "MYC", "BASP1")
ggsave("cellValueTable.tiff", width = 130, height = 40, units = "mm")

# export study list
library(cBioPortalData)
studies <- getStudies(cBioPortal())
studyList <- studies[which(studies$studyId %in% unique(allPatientData$studyId)),] %>%
  select(name, description, citation)
write.csv(studyList, "studyList.csv")