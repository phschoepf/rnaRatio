# Method script containing common code shared between all plotting scripts.
# @author Philemon Schöpf <philemon.schoepf@student.ubik.ac.at>

OUTPUT_PATH = "images/"
WIDTH = 2400
HEIGHT = 1100

# Arrange the plots and export to PNG.
# Generates 3 files: cell lines, patients, and a combined image.
exportToPng <- function(baseName, cellplot, patplot) {
  aligned <- align_plots(cellplot, patplot, align = "v", axis = "lr")
  png(paste0(OUTPUT_PATH, baseName, "_cells_annotation.png"), width = WIDTH, height = HEIGHT)
  plot(aligned[[1]])
  dev.off()
  png(paste0(OUTPUT_PATH, baseName, "_patients.png"), width = WIDTH, height = HEIGHT)
  plot(aligned[[2]])
  dev.off()
  png(paste0(OUTPUT_PATH, baseName, "_combined.png"), width = WIDTH, height = 2 * HEIGHT)
  grid.arrange(aligned[[1]],
               aligned[[2]],
               layout_matrix = rbind(1, 2))
  dev.off()
}


plotPatients <- function(in_data, gene1, gene2 = NULL, lims) {
  if(is.null(gene2)) {
    # only 1 gene, plot absolutes
    plotTitle <- paste0(gene1, " RNA expression")
    axisTitle <- "RNA expression, RNAseq v2 RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"
    plotBy <- expr(!!sym(gene1))
    
  } 
  else {
    # 2 genes, plot gene1 to gene2 ratio on log scale
    plotTitle <- paste0(gene1, " / ", gene2, " RNA expression ratio")
    axisTitle <- expression(paste(log[10], " RNA expression ratio, RNAseq v2 RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"))
    plotBy <- expr(log10(!!sym(gene1) / !!sym(gene2)))
  }
  
  plotdata <- in_data %>% 
    group_by(Origin) %>% 
    arrange(desc(!!plotBy, by_group = T)) %>%
    mutate(toPlot = !!plotBy)
  
  patplot <-
    ggplot(plotdata, aes(
      x = Origin,
      y = !!plotBy,
      na.rm = TRUE
    )) +
    coord_cartesian(ylim = lims) +
    theme_grey(base_size = 5) +
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
    geom_boxplot(alpha = 0.8,
                 outlier.shape = NA,
                 fatten = 1.2
    ) 
    
}


plotCells <- function(in_data, gene1, gene2 = NULL, lims) {
  if(is.null(gene2)) {
    # only 1 gene, plot absolutes
    plotTitle <- paste0(gene1, " RNA expression")
    axisTitle <- "RNA expression, RNAseq v2 RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"
    plotBy <- expr(!!sym(gene1))
    
  } 
  else {
    # 2 genes, plot gene1 to gene2 ratio on log scale
    plotTitle <- paste0(gene1, " / ", gene2, " RNA expression ratio")
    axisTitle <- expression(paste(log[10], " RNA expression ratio, RNAseq v2 RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"))
    plotBy <- expr(log10(!!sym(gene1) / !!sym(gene2)))
  }
  
  plotdata <- in_data %>% 
    group_by(Origin) %>% 
    arrange(desc(!!plotBy, by_group = T)) %>%
    mutate(toPlot = !!plotBy)
  
  cellplot <-
    ggplot(plotdata, aes(
      x = Origin,
      y = !!plotBy,
      na.rm = TRUE
    )) +
    coord_cartesian(ylim = lims) +
    theme_grey(base_size = 5) +
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
    geom_boxplot(alpha = 0.8,
                 outlier.shape = NA,
                 fatten = 1.2) +
    geom_point(data = selectedCellLines, color = "red", size = 2) +
    geom_text(
      data = selectedCellLines,
      aes(label = Name),
      nudge_x = 0.5,
      nudge_y = 0
    ) + 
    annotate(
      geom = "table",
      x = Inf,
      y = Inf,
      label = list(
        selectedCellLines %>% mutate(rounded = round(!!plotBy, 2)) %>% select(Name, rounded) 
      ),
      size = 0.7
    )
}
