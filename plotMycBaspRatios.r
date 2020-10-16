source("plotMethods.r")

# plotting cells
cellplot <-
  ggplot(cellLines, aes(
    y = reorder(Origin, MYC / BASP1, FUN = median),
    x = log10(MYC / BASP1),
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(-2, 4)), expand = F) +
  labs(title = "Cell lines",
       y = "",
       x = expression(paste(
         log[10](MYC / BASP1), " RNA expression, RNAseq RPKM"
       ))) +
  theme(axis.text.y = element_text(size = 20)) +
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
    nudge_x = 0.02,
    nudge_y = 0.65
  )

# plotting patients
patplot <-
  ggplot(allPatientData, aes(
    y = reorder(Origin, MYC / BASP1, FUN = median),
    x = log10(MYC / BASP1),
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(-2, 4)), expand = F) +
  labs(title = "Patient samples",
       y = "" ,
       x = expression(
         paste(
           log[10](MYC / BASP1),
           " RNA expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"
         )
       )) +
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
exportToPng(baseName = "myc_basp1_ratios", cellplot, patplot)
