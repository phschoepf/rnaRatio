# Method script containing common code shared between all plotting scripts.
# @author Philemon Schöpf <philemon.schoepf@student.ubik.ac.at>

library(tidyverse)
library(ggpubr)
library(ggrepel)

# Arrange the plots and export to PNG.
# Generates 3 files: cell lines, patients, and a combined image.
exportToPng <- function(baseName, cellplot, patplot) {
  
  OUTPUT_PATH = "images/"
  OUTPUT_FORMAT = "png"
  
  #image dims in mm (300 dpi)
  WIDTH = 100
  HEIGHT = 40
  
  aligned <- align_plots(cellplot, patplot, align = "v", axis = "lr")
  ggsave(paste0(OUTPUT_PATH, baseName, "_cells_annotation.", OUTPUT_FORMAT), plot = aligned[[1]], width = WIDTH, height = HEIGHT, units = "mm")
  ggsave(paste0(OUTPUT_PATH, baseName, "_patients.", OUTPUT_FORMAT), plot = aligned[[2]], width = WIDTH, height = HEIGHT, units = "mm")
  ggarrange(plotlist = aligned, nrow = 1, ncol = 2)
  ggsave(paste0(OUTPUT_PATH, baseName, "_combined.", OUTPUT_FORMAT), width = WIDTH, height = 2* HEIGHT, units = "mm")
}


plotPatients <- function(in_data, gene1, gene2 = NULL, lims) {
  if(is.null(gene2)) {
    # only 1 gene, plot absolutes
    plotTitle <- paste0(gene1, " RNA expression")
    axisTitle <- "RNA expression, RNAseq v2 RSEM"
    plotBy <- expr(!!sym(gene1))
    
  } 
  else {
    # 2 genes, plot gene1 to gene2 ratio on log scale
    plotTitle <- paste0(gene1, " / ", gene2, " RNA expression ratio")
    axisTitle <- expression(paste(log[10], " RNA expression ratio, RNAseq v2 RSEM")) #M (Batch normalized from Illumina HiSeq_RNASeqV2)
    plotBy <- expr(log10(!!sym(gene1) / !!sym(gene2)))
  }
  
  patplot <-
    ggplot(in_data, aes(
      x = reorder(Origin, -(!!plotBy), FUN = median),
      y = !!plotBy,
      na.rm = TRUE
    )) +
    coord_cartesian(ylim = lims) +
    theme_grey(base_size = 5) +
    rotate_x_text(30, hjust = 1) +
    labs(
      title = plotTitle,
      x = "",
      y = axisTitle
    ) +
    geom_jitter(
      shape = 16,
      position = position_jitter(w = 0.25, h = 0),
      size = 0.6,
      color = "#f39200"
    ) +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_boxplot(alpha = 0.4,
                 outlier.shape = NA,
                 fatten = 1.2
    ) 
    
}


plotCells <- function(in_data, gene1, gene2 = NULL, lims) {
  if(is.null(gene2)) {
    # only 1 gene, plot absolutes
    plotTitle <- paste0(gene1, " RNA expression")
    axisTitle <- "RNA expression, RNAseq RPKM"
    plotBy <- expr(!!sym(gene1))
    
  } 
  else {
    # 2 genes, plot gene1 to gene2 ratio on log scale
    plotTitle <- paste0(gene1, " / ", gene2, " RNA expression ratio")
    axisTitle <- expression(paste(log[10], " RNA expression ratio, RNAseq RPKM"))
    plotBy <- expr(log10(!!sym(gene1) / !!sym(gene2)))
  }
             
  cellplot <-
    ggplot(in_data, aes(
      x = reorder(Origin, -(!!plotBy), FUN = median),
      y = !!plotBy,
      na.rm = TRUE
    )) +
    coord_cartesian(ylim = lims) +
    theme_grey(base_size = 5) +
    theme(plot.margin = margin(l = 35))+
    rotate_x_text(30, hjust = 1) +
    labs(
      title = plotTitle,
      x = "",
      y = axisTitle
    ) +
    geom_jitter(
      shape = 16,
      position = position_jitter(w = 0.25, h = 0),
      size = 0.6,
      color = "#f39200"
    ) +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_boxplot(alpha = 0.5,
                 outlier.shape = NA,
                 fatten = 1.2) +
    geom_point(data = selectedCellLines, color = "red", size = 1) +
    geom_label_repel(
      data = selectedCellLines,
      aes(label = Name),
      size = 1.2,
      label.padding = unit(0.12, "lines"),
      nudge_x = 0.3,
      direction = "y"
    ) 
}

# Plot a table of expression values for certain cells. Input selected cell lines (as table from parseInputData.r),
# and 2 genes. The table will contain the 2 absolute expression values and the log-fold ratio.
plotCellsTable <- function(in_data, gene1, gene2) {
  
    firstgene <- expr(!!sym(gene1))
    secondgene <- expr(!!sym(gene2))
    ratio <- expr(log10(!!sym(gene1) / !!sym(gene2)))
    
  selectedCellsTable <- in_data %>% 
    transmute(Name, !!firstgene, !!secondgene, !!ratio) %>%
    mutate_if(is.numeric, round, digits = 2) %>% 
    arrange(Name)
  
  ggtexttable(t(selectedCellsTable), theme = ttheme("classic", base_size = 5))
}
