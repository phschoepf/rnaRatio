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

plotPatients <- function(data, genes) {
  geneExp <- substitute(genes)
  patplot <-
    ggplot(allPatientData, aes(
      y = reorder(Origin, MYC/PHB, FUN = median),
      x = log10(MYC / PHB),
      na.rm = TRUE
    )) +
    coord_cartesian(xlim = (c(-3,2)), expand = F) + 
    theme_classic(base_size = 20) +
    labs(title = "Patient samples",
         y = "" ,
         x = expression(paste(log[10](MYC / PHB), " RNA expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"
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
}
