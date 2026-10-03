# =========================================================
# Ranked lollipop plot of differentially expressed circRNAs
# Highlight hsa_circRNA_000554
# =========================================================

library(ggplot2)
library(dplyr)

outdir <- "18_circRNA_Ranked_Lollipop_Plot"
dir.create(outdir, showWarnings = FALSE)

circ_df <- data.frame(
  logFC = c(
    -2.01199, -1.84330, -1.71891, -1.44681, -1.43287,
    -1.41095, -1.36077, -1.32382, -1.15750, -1.11105,
    -1.09949
  ),
  P.Value = c(
    0.044298, 0.030574, 0.031712, 0.044010, 0.022428,
    0.011127, 0.004637, 0.001512, 0.045276, 0.042219,
    0.024173
  ),
  circRNA = c(
    "ASCRP003643", "ASCRP002964", "ASCRP000212", "ASCRP001740",
    "ASCRP004906", "ASCRP005049", "ASCRP000064", "ASCRP001375",
    "ASCRP001410", "ASCRP003932", "ASCRP003189"
  ),
  Rank = 1:11,
  DisplayName = c(
    "ASCRP003643", "ASCRP002964", "ASCRP000212", "ASCRP001740",
    "ASCRP004906", "ASCRP005049", "hsa_circRNA_000554",
    "ASCRP001375", "ASCRP001410", "ASCRP003932", "ASCRP003189"
  ),
  stringsAsFactors = FALSE
)

circ_df <- circ_df %>%
  mutate(
    Candidate = ifelse(DisplayName == "hsa_circRNA_000554",
                       "hsa_circRNA_000554", "Other circRNAs"),
    negLog10P = -log10(P.Value),
    Label = paste0("Rank ", Rank)
  ) %>%
  arrange(logFC)

circ_df$DisplayName <- factor(circ_df$DisplayName, levels = circ_df$DisplayName)

p_lollipop <- ggplot(
  circ_df,
  aes(
    x = logFC,
    y = DisplayName
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = logFC,
      y = DisplayName,
      yend = DisplayName
    ),
    color = "grey65",
    linewidth = 0.9
  ) +
  geom_point(
    aes(
      color = Candidate,
      size = negLog10P
    ),
    alpha = 0.95
  ) +
  geom_text(
    aes(
      label = Label
    ),
    hjust = 1.15,
    size = 3.2,
    fontface = "bold",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "hsa_circRNA_000554" = "#E64B35",
      "Other circRNAs" = "#4DBBD5"
    )
  ) +
  scale_size_continuous(
    range = c(3.5, 7),
    name = "-log10(P value)"
  ) +
  scale_x_continuous(
    limits = c(-2.25, 0.15),
    breaks = seq(-2.0, 0, 0.5)
  ) +
  labs(
    title = "Ranked differentially expressed circRNAs in TNBC",
    x = "log2 fold change",
    y = "",
    color = ""
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.title.x = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_lollipop)

ggsave(
  file.path(outdir, "Figure_circRNA_Ranked_Lollipop_FINAL.png"),
  p_lollipop,
  width = 9,
  height = 6.5,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure_circRNA_Ranked_Lollipop_FINAL.tiff"),
  p_lollipop,
  width = 9,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure_circRNA_Ranked_Lollipop_FINAL.pdf"),
  p_lollipop,
  width = 9,
  height = 6.5
)

write.csv(
  circ_df,
  file.path(outdir, "Table_Differentially_Expressed_circRNAs_FINAL.csv"),
  row.names = FALSE
)

list.files(outdir)
# =========================================================
# Ranked lollipop plot of differentially expressed circRNAs
# Final publication-ready version
# Highlight hsa_circRNA_000554
# =========================================================

library(ggplot2)
library(dplyr)

outdir <- "18_circRNA_Ranked_Lollipop_Plot"
dir.create(outdir, showWarnings = FALSE)

# -----------------------------
# Input data
# -----------------------------

circ_df <- data.frame(
  logFC = c(
    -2.01199, -1.84330, -1.71891, -1.44681, -1.43287,
    -1.41095, -1.36077, -1.32382, -1.15750, -1.11105,
    -1.09949
  ),
  P.Value = c(
    0.044298, 0.030574, 0.031712, 0.044010, 0.022428,
    0.011127, 0.004637, 0.001512, 0.045276, 0.042219,
    0.024173
  ),
  circRNA = c(
    "ASCRP003643", "ASCRP002964", "ASCRP000212", "ASCRP001740",
    "ASCRP004906", "ASCRP005049", "ASCRP000064", "ASCRP001375",
    "ASCRP001410", "ASCRP003932", "ASCRP003189"
  ),
  Rank = 1:11,
  DisplayName = c(
    "ASCRP003643", "ASCRP002964", "ASCRP000212", "ASCRP001740",
    "ASCRP004906", "ASCRP005049", "hsa_circRNA_000554",
    "ASCRP001375", "ASCRP001410", "ASCRP003932", "ASCRP003189"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Prepare data
# -----------------------------

circ_df <- circ_df %>%
  mutate(
    Candidate = ifelse(
      DisplayName == "hsa_circRNA_000554",
      "hsa_circRNA_000554",
      "Other circRNAs"
    ),
    negLog10P = -log10(P.Value)
  ) %>%
  arrange(logFC)

circ_df$DisplayName <- factor(
  circ_df$DisplayName,
  levels = circ_df$DisplayName
)

candidate_df <- circ_df %>%
  filter(DisplayName == "hsa_circRNA_000554")

# -----------------------------
# Plot
# -----------------------------

p_lollipop_final <- ggplot(
  circ_df,
  aes(
    x = logFC,
    y = DisplayName
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = logFC,
      y = DisplayName,
      yend = DisplayName
    ),
    color = "grey65",
    linewidth = 0.9
  ) +
  geom_point(
    data = circ_df %>% filter(Candidate == "Other circRNAs"),
    aes(
      size = negLog10P,
      color = Candidate
    ),
    alpha = 0.9
  ) +
  geom_point(
    data = candidate_df,
    aes(
      size = negLog10P,
      color = Candidate
    ),
    alpha = 1
  ) +
  geom_point(
    data = candidate_df,
    aes(
      x = logFC,
      y = DisplayName
    ),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 1.3,
    size = 8,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = candidate_df,
    aes(
      x = logFC - 0.42,
      y = DisplayName,
      label = "Selected candidate"
    ),
    size = 3.6,
    fontface = "bold",
    color = "black",
    fill = "white",
    label.size = 0.3,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = candidate_df,
    aes(
      x = logFC - 0.18,
      xend = logFC - 0.02,
      y = DisplayName,
      yend = DisplayName
    ),
    arrow = arrow(length = unit(0.18, "cm")),
    color = "black",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "hsa_circRNA_000554" = "#E64B35",
      "Other circRNAs" = "#4DBBD5"
    )
  ) +
  scale_size_continuous(
    range = c(3.5, 7),
    name = "-log10(P value)"
  ) +
  scale_x_continuous(
    limits = c(-2.35, 0.15),
    breaks = seq(-2.0, 0, 0.5)
  ) +
  labs(
    title = "Ranked differentially expressed circRNAs in TNBC",
    x = "log2 fold change",
    y = "",
    color = ""
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.title.x = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_lollipop_final)

# -----------------------------
# Save files
# -----------------------------

ggsave(
  file.path(outdir, "Figure_circRNA_Ranked_Lollipop_PUBLICATION_FINAL.png"),
  p_lollipop_final,
  width = 9,
  height = 6.5,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure_circRNA_Ranked_Lollipop_PUBLICATION_FINAL.tiff"),
  p_lollipop_final,
  width = 9,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure_circRNA_Ranked_Lollipop_PUBLICATION_FINAL.pdf"),
  p_lollipop_final,
  width = 9,
  height = 6.5
)

write.csv(
  circ_df,
  file.path(outdir, "Table_Differentially_Expressed_circRNAs_PUBLICATION_FINAL.csv"),
  row.names = FALSE
)

list.files(outdir)