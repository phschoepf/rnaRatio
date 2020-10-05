library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)

source("parseInputData.r")

# plotting cells
cellplot <-
  ggplot(cellLines, aes(
    y = reorder(Origin, MYC / PHB, FUN = median),
    x = log10(MYC / PHB),
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(-3,2))) +
  labs(title = "Cell lines", y = "", x = expression(paste(log[10](MYC / PHB), ", RNAseq RPKM"))) +
  stat_boxplot(geom = "errorbar", width = 0.6, lwd = 1) +
  geom_boxplot(outlier.shape = NA,
               lwd = 1,
               fatten = 1.2) +
  geom_jitter(
    shape = 16,
    position = position_jitter(w = 0, h = 0.15),
    size = 0.9,
    color = "#f39200") +
  geom_point(data = selectedCellLines, color = "red", size = 4) +
  geom_text(data = selectedCellLines, aes(label = Name), nudge_x = 0.02, nudge_y = 0.65)

# plotting patients
patplot <-
  ggplot(allPatientData, aes(
    y = reorder(Origin, MYC / PHB, FUN = median),
    x = log10(MYC / PHB),
    na.rm = TRUE
  )) +
  coord_cartesian(xlim = (c(-3,2))) + 
  theme_classic(base_size = 20) +
  labs(title = "Patient samples",
       y = "" ,
       x = expression(paste(log[10](MYC / PHB), ", RNAseq v2 RSEM"
       ))) +
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
png("images/myc_phb_ratios_cells_annotation.png", width = 2400, height = 1100)
plot(aligned[[1]])
dev.off()
png("images/myc_phb_ratios_patients.png", width = 2400, height = 1100)
plot(aligned[[2]])
dev.off()
png("images/myc_phb_ratios_combined.png", width = 2400, height = 2400)
grid.arrange(aligned[[1]],
             aligned[[2]],
             textGrob("Cell lines: RNAseq RPKM, patient data: RNAseq v2 RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)", x = 1, hjust = 1),
             layout_matrix = rbind(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3))
dev.off()