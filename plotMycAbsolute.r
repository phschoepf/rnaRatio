source("plotMethods.r")

# plotting cells
cellplot <-
  ggplot(cellLines, aes(
    y = reorder(Origin, MYC, FUN = median),
    x = MYC,
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(0, 400)), expand = F) +
  labs(title = "Cell lines",
       y = "" ,
       x = "MYC RNA expression, RNAseq RPKM") +
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
    nudge_x = 4,
    nudge_y = 0.25
  ) +
  annotate(
    geom = "table",
    x = Inf,
    y = Inf,
    label = list(
      selectedCellLines %>% select(Name, MYC) %>% mutate(MYC = round(MYC, 2))
    ),
    size = 10
  )



# plotting patients
patplot <-
  ggplot(allPatientData, aes(
    y = reorder(Origin, MYC, FUN = "median"),
    x = MYC,
    na.rm = TRUE
  )) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = (c(0, 15000)), expand = F) +
  labs(title = "Patient samples",
       y = "" ,
       x = "MYC RNA expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)") +
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
exportToPng(baseName = "myc_absolute", cellplot, patplot)