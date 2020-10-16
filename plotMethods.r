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

plotPatients <- function(data, is_ratio, xlab) {
  patplot <-
    ggplot(data, aes(
      y = reorder(Origin, ~data[ncol(data)], FUN = median), # plot the last col (possibly TODO)
      x = ~log10(data[ncol(data)]),
      na.rm = TRUE
    )) +
    coord_cartesian(xlim = (c(-3,2)), expand = F) + 
    theme_grey(base_size = 20) +
    labs(title = "Patient samples", y = "" , x = xlab) +
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

plotCells <- function(data, is_ratio, xlab) {
  cellplot <-
    ggplot(cellLines, aes(
      y = reorder(Origin, data[ncol(data)], FUN = median), # plot the last col (possibly TODO)
      x = log10(!!data[ncol(data)]),
      na.rm = TRUE
    )) +
    coord_cartesian(xlim = (c(-3,2)), expand = F) +
    theme_grey(base_size = 20) +
    labs(title = "Cell lines", y = "", x = xlab) +
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
}
