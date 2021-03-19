# Plotting settings and methods for RNA expression ratio plots.
# Author: Philemon Schoepf <philemon.schoepf@student.ubik.ac.at>
# Date: 2021-03-19

# Imports -----------------------------------------------------------------

library(cowplot)
library(tidyverse)
library(ggpubr)
library(ggrepel)

# Helper functions -------------------------------------------------------------

plotPatients <- function(.data, gene1, gene2 = NULL, ...) {
  if(is.null(gene2)) {
    # only 1 gene, plot absolutes
    plotTitle <- paste0(gene1, " RNA expression")
    axisTitle <- "RNA expression, RNAseq v2 RSEM"
    plotBy <- expr(!!sym(gene1))

    #outlier removal
    valuesOnly<- .data[[gene1]]
    boxplotStats <- boxplot.stats(valuesOnly)$stats
    ylim <- c(0, boxplotStats[4]+(boxplotStats[4]-boxplotStats[2])*2)

  }
  else {
    # 2 genes, plot gene1 to gene2 ratio on log scale
    plotTitle <- paste0(gene1, " / ", gene2, " RNA expression ratio")
    axisTitle <- expression(paste(log[10], " RNA expression ratio, RNAseq v2 RSEM"))
    plotBy <- expr(log10(!!sym(gene1) / !!sym(gene2)))
    ylim <- NULL # let R do the y-limit for ratio plots
  }

  patplot <-
    ggplot(.data, aes(
      x = reorder(origin, -(!!plotBy), FUN = median),
      y = !!plotBy,
      na.rm = TRUE
    )) +
    coord_cartesian(ylim = ylim) +
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

plotCells <- function(.data, gene1, gene2 = NULL, ...) {
  if(is.null(gene2)) {
    # only 1 gene, plot absolutes
    plotTitle <- paste0(gene1, " RNA expression")
    axisTitle <- "RNA expression, RNAseq RPKM"
    plotBy <- expr(!!sym(gene1))

    #outlier removal
    valuesOnly<- .data[[gene1]]
    boxplotStats <- boxplot.stats(valuesOnly)$stats
    ylim <- c(0, boxplotStats[4]+(boxplotStats[4]-boxplotStats[2])*2)

  }
  else {
    # 2 genes, plot gene1 to gene2 ratio on log scale
    plotTitle <- paste0(gene1, " / ", gene2, " RNA expression ratio")
    axisTitle <- expression(paste(log[10], " RNA expression ratio, RNAseq RPKM"))
    plotBy <- expr(log10(!!sym(gene1) / !!sym(gene2)))
    ylim <- NULL # let R do the y-limit for ratio plots
  }

  cellplot <-
    ggplot(.data, aes(
      x = reorder(origin, -(!!plotBy), FUN = median),
      y = !!plotBy,
      na.rm = TRUE
    )) +
    coord_cartesian(ylim = ylim) +
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
    geom_point(data = selectedCellData, color = "red", size = 1) +
    geom_label_repel(
      data = selectedCellData,
      aes(label = Name),
      size = 1.2,
      label.padding = unit(0.12, "lines"),
      nudge_x = 0.3,
      direction = "y"
    )
}

# Plot a table of expression values for certain cells. Input selected cell lines (as table from parseInputData.r),
# and 2 genes. The table will contain the 2 absolute expression values and the log-fold ratio.
plotCellsTable <- function(.data, gene1, gene2) {

  ratio <- expr(!!sym(gene1) / !!sym(gene2))

  selectedCellsTable <- .data %>%
    dplyr::select(Name, {{gene1}}, {{gene2}}) %>%
    dplyr::mutate(ratio = !!ratio) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    arrange(-ratio)

  ggtexttable(t(selectedCellsTable), theme = ttheme("classic"))
}

# Main function -----------------------------------------------------------------

plotRatio <- function(gene1, gene2) {

  ratioCell <- plotCells(filteredCellData, gene1, gene2, ylim =  c(-2,4))
  ratioPat <- plotPatients(filteredPatientData, gene1, gene2, ylim = c(-2,4))

  abs1Cell <- plotCells(filteredCellData, gene1)
  abs1Pat <- plotPatients(filteredPatientData, gene1)

  abs2Cell <- plotCells(filteredCellData, gene2)
  abs2Pat <- plotPatients(filteredPatientData, gene2)

  # align plots to each other
  alignedPlots <- align_plots(ratioCell, ratioPat,
                              abs1Cell, abs1Pat,
                              abs2Cell, abs2Pat,
                              align = "h", axis = "tb")
  combinedPlot <- ggarrange(plotlist = alignedPlots, nrow = 3, ncol = 2)

  cellsTable <- plotCellsTable(selectedCellData, gene1, gene2)

  #write to disk
  if(DO_WRITE == T) {
    ggsave(paste0(OUTPUT_PATH, gene1, "_", gene2, "_cells.", OUTPUT_FORMAT), plot = alignedPlots[[1]], width = WIDTH, height = HEIGHT, units = "mm")
    ggsave(paste0(OUTPUT_PATH, gene1, "_", gene2, "_patients.", OUTPUT_FORMAT), plot = alignedPlots[[2]], width = WIDTH, height = HEIGHT, units = "mm")

    ggsave(paste0(OUTPUT_PATH, gene1, "_cells.", OUTPUT_FORMAT), plot = alignedPlots[[3]], width = WIDTH, height = HEIGHT, units = "mm")
    ggsave(paste0(OUTPUT_PATH, gene1,"_patients.", OUTPUT_FORMAT), plot = alignedPlots[[4]], width = WIDTH, height = HEIGHT, units = "mm")

    ggsave(paste0(OUTPUT_PATH, gene2, "_cells.", OUTPUT_FORMAT), plot = alignedPlots[[5]], width = WIDTH, height = HEIGHT, units = "mm")
    ggsave(paste0(OUTPUT_PATH, gene2, "_patients.", OUTPUT_FORMAT), plot = alignedPlots[[6]], width = WIDTH, height = HEIGHT, units = "mm")

    ggsave(paste0(OUTPUT_PATH, gene1, "_", gene2, "_combined.", OUTPUT_FORMAT), plot = combinedPlot, width = 170, height = 190, units = "mm")

    ggsave(paste0(OUTPUT_PATH, gene2, "_cellValueTable.", OUTPUT_FORMAT), plot = cellsTable, width = 130, height = 40, units = "mm")
  }
  # or output directly to a viewport
  else {
    for (plot in alignedPlots) {
      plot(plot)
    }
    plot(combinedPlot)
    plot(cellsTable)
  }
}

