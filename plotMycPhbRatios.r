source("plotMethods.r")

# plotting cells
cellplot <-
  ggplot(cellLines, aes(
    y = reorder(Origin, MYC / PHB, FUN = median),
    x = log10(MYC / PHB),
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(-3,2)), expand = F) +
  labs(title = "Cell lines", y = "", x = expression(paste(log[10](MYC / PHB), " RNA expression, RNAseq RPKM"))) +
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

patplot <- plotPatients(allPatientData, "MYC/PHB")

exportToPng(baseName = "myc_phb_ratios", cellplot, patplot)