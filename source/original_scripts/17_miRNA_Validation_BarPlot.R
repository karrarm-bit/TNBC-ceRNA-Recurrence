install.packages("tidyr")
# =========================================================
# miRNA validation bar plot
# Tumor vs Normal mean expression
# =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)

outdir <- "17_miRNA_Validation_BarPlot"
dir.create(outdir, showWarnings = FALSE)

# -----------------------------
# Input validation results
# -----------------------------

val_df <- data.frame(
  miRNA = c(
    "hsa-miR-590-5p",
    "hsa-miR-182-5p",
    "hsa-miR-452-5p"
  ),
  Validation_ID = c(
    "hsa-miR-590-5p",
    "hsa-miR-182",
    "hsa-miR-452"
  ),
  Tumor_mean = c(
    70.60078,
    38.65105,
    15.84012
  ),
  Normal_mean = c(
    28.44018,
    6.065257,
    13.47831
  ),
  Direction = c(
    "Up",
    "Up",
    "Up"
  ),
  P_value = c(
    0.00105,
    1.34E-06,
    0.876
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Prepare data for plotting
# -----------------------------

plot_df <- val_df %>%
  select(miRNA, Tumor_mean, Normal_mean, P_value) %>%
  pivot_longer(
    cols = c(Normal_mean, Tumor_mean),
    names_to = "Group",
    values_to = "Mean_expression"
  ) %>%
  mutate(
    Group = ifelse(Group == "Normal_mean", "Normal", "Tumor"),
    Group = factor(Group, levels = c("Normal", "Tumor")),
    miRNA = factor(
      miRNA,
      levels = c(
        "hsa-miR-590-5p",
        "hsa-miR-182-5p",
        "hsa-miR-452-5p"
      )
    )
  )

# Significance labels
sig_df <- val_df %>%
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = pmax(Tumor_mean, Normal_mean) * 1.15
  )

# -----------------------------
# Bar plot
# -----------------------------

p_val <- ggplot(
  plot_df,
  aes(
    x = miRNA,
    y = Mean_expression,
    fill = Group
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    data = sig_df,
    aes(
      x = miRNA,
      y = y_pos,
      label = paste0("P = ", signif(P_value, 3), "\n", Significance)
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#4DBBD5",
      "Tumor" = "#E64B35"
    )
  ) +
  labs(
    title = "Independent validation of candidate miRNAs",
    x = "",
    y = "Mean expression",
    fill = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      face = "bold"
    ),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_val)

# -----------------------------
# Save files
# -----------------------------

ggsave(
  file.path(outdir, "Figure_miRNA_Validation_BarPlot_FINAL.png"),
  p_val,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure_miRNA_Validation_BarPlot_FINAL.tiff"),
  p_val,
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure_miRNA_Validation_BarPlot_FINAL.pdf"),
  p_val,
  width = 8,
  height = 6
)

write.csv(
  val_df,
  file.path(outdir, "Table_miRNA_Validation_Results_FINAL.csv"),
  row.names = FALSE
)

list.files(outdir)