library(rlang)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)
library(ggpmisc)

#source("parseInputData.r")
source("plotMethods.r")

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

#combined large plot
tiff("allPlots.tiff", width = 170, height = 200, units = "mm", res = 300)
grid.arrange(myc_cell, myc_pat,
             basp1_cell, basp1_pat,
             mycbasp_cell, mycbasp_pat,
             layout_matrix = cbind(c(1,3,5), c(2,4,6))
)
dev.off()
#source("RPKMvsRSEMdemo.r") does not work currently w/ new APIrequest input system