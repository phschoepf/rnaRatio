library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)
library(ggpmisc)

source("parseInputData.r")

# plotting cells
cellplot <-
  ggplot(cellLines, aes(
    y = reorder(Origin, -BASP1, FUN = median),
    x = BASP1,
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(0, 700)), expand = F) +
  labs(title = "Cell lines",
       y = "" ,
       x = "BASP1 RNA expression, RNAseq RPKM") +
  stat_boxplot(geom = "errorbar", width = 0.6, lwd = 1) +
  geom_boxplot(outlier.shape = NA,
               lwd = 1,
               fatten = 1.2) +
  geom_jitter(
    shape = 16,
    position = position_jitter(w = 0, h = 0.15),
    size = 0.9,
    color = "#f39200"
  ) +
  geom_point(data = selectedCellLines, color = "red", size = 4) +
  geom_text(
    data = selectedCellLines,
    aes(label = Name),
    nudge_x = 8,
    nudge_y = 0.25
  ) +
  annotate(
    geom = "table",
    x = Inf,
    y = Inf,
    label = list(
      selectedCellLines %>% select(Name, BASP1) %>% mutate(BASP1 = round(BASP1, 2))
    ),
    size = 10
  )

# plotting patients
patplot <-
  ggplot(allPatientData, aes(
    y = reorder(Origin, -BASP1, FUN = "median"),
    x = BASP1,
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(-0, 15000)), expand = F) +
  labs(title = "Patient samples",
       y = "" ,
       x = "BASP1 RNA expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)") +
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
aligned <- align_plots(cellplot, patplot, align = "v", axis = "lr")
png(
  "images/basp1_absolute_cells_annotation.png",
  width = 2400,
  height = 1100
)
plot(aligned[[1]])
dev.off()
png("images/basp1_absolute_patients.png",
    width = 2400,
    height = 1100)
plot(aligned[[2]])
dev.off()
png("images/basp1_absolute_combined.png",
    width = 2400,
    height = 2400)
grid.arrange(aligned[[1]],
             aligned[[2]],
             layout_matrix = rbind(1, 2))
dev.off()