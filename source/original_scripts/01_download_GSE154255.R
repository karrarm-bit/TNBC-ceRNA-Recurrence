if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)

gse <- getGEO("GSE154255", GSEMatrix = TRUE)
gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)

dim(expr)
head(pheno[, c("title", "geo_accession")])

write.csv(expr, "data/GSE154255_expression_matrix.csv")
write.csv(pheno, "data/GSE154255_sample_info.csv")
# ============================================================
# GSE154255 - TNBC cancer vs adjacent normal
# Check hsa-miR-590-5p and hsa-miR-182-5p
# ============================================================

# 1) Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c("GEOquery", "limma", "dplyr", "ggplot2", "pheatmap")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)

# 2) Download GEO dataset
gse_list <- getGEO("GSE154255", GSEMatrix = TRUE)
gse <- gse_list[[1]]

expr <- exprs(gse)
pheno <- pData(gse)

# 3) Inspect sample titles
pheno[, c("geo_accession", "title")]
# 4) Select TNBC samples only
tnbc_idx <- grepl("TNBC", pheno$title, ignore.case = TRUE)

expr_tnbc <- expr[, tnbc_idx]
pheno_tnbc <- pheno[tnbc_idx, ]

pheno_tnbc[, c("geo_accession", "title")]

# 5) Define groups: Normal vs Cancer
group <- ifelse(
  grepl("cancer", pheno_tnbc$title, ignore.case = TRUE),
  "TNBC_Cancer",
  "Normal"
)

group <- factor(group, levels = c("Normal", "TNBC_Cancer"))
group
table(group)
# 6) limma differential expression
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_tnbc, design)

contrast_matrix <- makeContrasts(
  TNBC_vs_Normal = TNBC_Cancer - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  coef = "TNBC_vs_Normal",
  number = Inf,
  adjust.method = "BH"
)

head(deg)
deg$miRNA <- rownames(deg)

upregulated_pvalue <- deg %>%
  filter(
    logFC > 0,
    P.Value < 0.05
  ) %>%
  arrange(P.Value)

head(upregulated_pvalue, 20)

write.csv(
  upregulated_pvalue,
  "TNBC_upregulated_miRNAs_Pvalue_only.csv",
  row.names = FALSE
)
selected_miRNAs_pvalue <- deg %>%
  filter(
    miRNA %in% c("hsa-miR-590-5p", "hsa-miR-182-5p"),
    logFC > 0,
    P.Value < 0.05
  ) %>%
  select(miRNA, logFC, AveExpr, P.Value, adj.P.Val)

selected_miRNAs_pvalue
selected_miRNAs_pvalue <- deg %>%
  filter(
    miRNA %in% c("hsa-miR-590-5p", "hsa-miR-182-5p"),
    logFC > 0,
    P.Value < 0.05
  ) %>%
  select(miRNA, logFC, AveExpr, P.Value, adj.P.Val)

selected_miRNAs_pvalue
# Full DEG table
write.csv(
  deg,
  "results/GSE154255_TNBC_vs_Normal_all_miRNAs.csv",
  row.names = FALSE
)

# Upregulated based on P.Value only
upregulated_pvalue <- deg %>%
  filter(logFC > 0, P.Value < 0.05) %>%
  arrange(P.Value)

write.csv(
  upregulated_pvalue,
  "results/GSE154255_TNBC_upregulated_Pvalue_only.csv",
  row.names = FALSE
)

# Selected miRNAs only
selected_miRNAs <- deg %>%
  filter(miRNA %in% c("hsa-miR-590-5p", "hsa-miR-182-5p")) %>%
  select(miRNA, logFC, AveExpr, t, P.Value, adj.P.Val, B)

selected_miRNAs

write.csv(
  selected_miRNAs,
  "results/GSE154255_selected_miRNAs_590_182.csv",
  row.names = FALSE
)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
write.csv(
  deg,
  "results/GSE154255_TNBC_vs_Normal_all_miRNAs.csv",
  row.names = FALSE
)
list.files()
# ============================================================
# GSE154255 - TNBC vs Normal
# Plots only, no file saving
# ============================================================

library(dplyr)
library(ggplot2)
library(pheatmap)

# Make sure the miRNA column exists
deg$miRNA <- rownames(deg)

# ============================================================
# 1) Results for selected miRNAs
# ============================================================

selected_miRNAs <- deg %>%
  filter(miRNA %in% c("hsa-miR-590-5p", "hsa-miR-182-5p")) %>%
  select(miRNA, logFC, AveExpr, t, P.Value, adj.P.Val, B)

selected_miRNAs
# ============================================================
# 2) Boxplot for selected miRNAs
# ============================================================

plot_miRNA_boxplot <- function(mirna_name) {
  
  hit <- grep(mirna_name, rownames(expr_tnbc), ignore.case = TRUE, value = TRUE)
  
  if (length(hit) == 0) {
    message("Not found: ", mirna_name)
    return(NULL)
  }
  
  hit <- hit[1]
  
  df <- data.frame(
    Sample = colnames(expr_tnbc),
    Expression = as.numeric(expr_tnbc[hit, ]),
    Group = group,
    miRNA = hit
  )
  
  print(df)
  
  p <- ggplot(df, aes(x = Group, y = Expression)) +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 3) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0(hit, " expression"),
      subtitle = "GSE154255: TNBC cancer vs adjacent normal",
      x = "",
      y = "Expression level"
    )
  
  print(p)
}

plot_miRNA_boxplot("hsa-miR-590-5p")
plot_miRNA_boxplot("hsa-miR-182-5p")
# ============================================================
# Colors
# ============================================================

my_colors <- c(
  "Normal" = "blue",
  "TNBC_Cancer" = "red"
)
plot_miRNA_boxplot <- function(mirna_name) {
  
  hit <- grep(mirna_name, rownames(expr_tnbc), ignore.case = TRUE, value = TRUE)
  
  if (length(hit) == 0) {
    message("Not found: ", mirna_name)
    return(NULL)
  }
  
  hit <- hit[1]
  
  df <- data.frame(
    Sample = colnames(expr_tnbc),
    Expression = as.numeric(expr_tnbc[hit, ]),
    Group = group,
    miRNA = hit
  )
  
  print(df)
  
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Group), width = 0.12, size = 3) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0(hit, " expression"),
      subtitle = "GSE154255: TNBC cancer vs adjacent normal",
      x = "",
      y = "Expression level"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p)
}

plot_miRNA_boxplot("hsa-miR-590-5p")
plot_miRNA_boxplot("hsa-miR-182-5p")
pca <- prcomp(t(expr_tnbc), scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(expr_tnbc),
  Group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

percentVar <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA plot",
    subtitle = "GSE154255: TNBC cancer vs adjacent normal",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_pca)
# ============================================================
# PCA plot - fixed version
# Remove miRNAs with zero variance before PCA
# ============================================================

# Calculate variance for each miRNA across TNBC samples
row_vars <- apply(expr_tnbc, 1, var)

# Keep only miRNAs with variance greater than zero
expr_tnbc_pca <- expr_tnbc[row_vars > 0, ]

# Run PCA
pca <- prcomp(t(expr_tnbc_pca), scale. = TRUE)

# Prepare PCA dataframe
pca_df <- data.frame(
  Sample = colnames(expr_tnbc_pca),
  Group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

percentVar <- round(100 * summary(pca)$importance[2, 1:2], 1)

# PCA plot
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA plot",
    subtitle = "GSE154255: TNBC cancer vs adjacent normal",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_pca)

pca_df
# ============================================================
# PCA plot - stronger cleaning version
# ============================================================

expr_tnbc_pca <- expr_tnbc

# Remove rows with NA or infinite values
expr_tnbc_pca <- expr_tnbc_pca[
  apply(expr_tnbc_pca, 1, function(x) all(is.finite(x))),
]

# Remove rows with zero variance
row_vars <- apply(expr_tnbc_pca, 1, var)
expr_tnbc_pca <- expr_tnbc_pca[row_vars > 0, ]

# Run PCA
pca <- prcomp(t(expr_tnbc_pca), scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(expr_tnbc_pca),
  Group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

percentVar <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA plot",
    subtitle = "GSE154255: TNBC cancer vs adjacent normal",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_pca)
pca_df
# ============================================================
# Cleaner PCA plot - no sample labels
# ============================================================

p_pca_clean <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA plot",
    subtitle = "GSE154255: TNBC cancer vs adjacent normal",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_pca_clean)
# ============================================================
# Differential Expression Volcano Plot
# Based on unadjusted P.Value < 0.05
# ============================================================

library(dplyr)
library(ggplot2)

deg$miRNA <- rownames(deg)

volcano_df <- deg %>%
  mutate(
    Regulation = case_when(
      logFC > 0 & P.Value < 0.05 ~ "Upregulated",
      logFC < 0 & P.Value < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    minus_log10_p = -log10(P.Value),
    label = ifelse(
      miRNA %in% c("hsa-miR-590-5p", "hsa-miR-182-5p"),
      miRNA,
      ""
    )
  )

volcano_colors <- c(
  "Upregulated" = "red",
  "Downregulated" = "blue",
  "Not significant" = "gray70"
)

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = minus_log10_p, color = Regulation)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(
    aes(label = label),
    vjust = -0.8,
    size = 4,
    color = "black"
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Differential Expression Volcano Plot",
    subtitle = "GSE154255: TNBC cancer vs adjacent normal",
    x = "logFC",
    y = "-log10(P.Value)",
    color = "Regulation"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_volcano)
# ============================================================
# Clean Volcano Plot highlighting only selected miRNAs
# ============================================================

library(dplyr)
library(ggplot2)

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

library(ggrepel)

deg$miRNA <- rownames(deg)

target_miRNAs <- c("hsa-miR-590-5p", "hsa-miR-182-5p")

volcano_selected <- deg %>%
  mutate(
    minus_log10_p = -log10(P.Value),
    Highlight = ifelse(miRNA %in% target_miRNAs, "Selected miRNAs", "Other miRNAs"),
    Label = ifelse(miRNA %in% target_miRNAs, miRNA, "")
  )

p_volcano_selected <- ggplot(volcano_selected, aes(x = logFC, y = minus_log10_p)) +
  geom_point(
    data = subset(volcano_selected, Highlight == "Other miRNAs"),
    color = "gray75",
    size = 1.8,
    alpha = 0.6
  ) +
  geom_point(
    data = subset(volcano_selected, Highlight == "Selected miRNAs"),
    color = "red",
    size = 4.5,
    alpha = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    aes(label = Label),
    size = 4.5,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.4
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Differential Expression Volcano Plot",
    subtitle = "Highlighted candidate upregulated miRNAs in TNBC",
    x = "logFC",
    y = "-log10(P.Value)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )

print(p_volcano_selected)
# ============================================================
# Zoomed Volcano Plot focusing on upregulated selected miRNAs
# ============================================================

p_volcano_zoom <- p_volcano_selected +
  coord_cartesian(
    xlim = c(0, max(volcano_selected$logFC, na.rm = TRUE) + 5),
    ylim = c(0, max(volcano_selected$minus_log10_p, na.rm = TRUE) + 1)
  ) +
  labs(
    title = "Selected Upregulated miRNAs",
    subtitle = "hsa-miR-590-5p and hsa-miR-182-5p in TNBC"
  )

print(p_volcano_zoom)
# ============================================================
# Volcano plot: 2 selected upregulated + 2 top downregulated
# Others in gray
# ============================================================

library(dplyr)
library(ggplot2)

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# Make sure miRNA names exist
deg$miRNA <- rownames(deg)

# Your selected upregulated miRNAs
up_targets <- c("hsa-miR-590-5p", "hsa-miR-182-5p")

# Automatically select top 2 downregulated miRNAs
down_targets <- deg %>%
  filter(logFC < 0, P.Value < 0.05) %>%
  arrange(logFC) %>%
  slice(1:2) %>%
  pull(miRNA)

down_targets
# Combine selected miRNAs
selected_targets <- c(up_targets, down_targets)

volcano_selected4 <- deg %>%
  mutate(
    minus_log10_p = -log10(P.Value),
    Category = case_when(
      miRNA %in% up_targets ~ "Upregulated",
      miRNA %in% down_targets ~ "Downregulated",
      TRUE ~ "Other"
    ),
    Label = ifelse(miRNA %in% selected_targets, miRNA, "")
  )

p_volcano_4 <- ggplot(volcano_selected4, aes(x = logFC, y = minus_log10_p)) +
  
  # Other miRNAs in gray
  geom_point(
    data = subset(volcano_selected4, Category == "Other"),
    color = "gray70",
    size = 2,
    alpha = 0.6
  ) +
  
  # Downregulated in blue
  geom_point(
    data = subset(volcano_selected4, Category == "Downregulated"),
    color = "blue",
    size = 4,
    alpha = 1
  ) +
  
  # Upregulated in red
  geom_point(
    data = subset(volcano_selected4, Category == "Upregulated"),
    color = "red",
    size = 4,
    alpha = 1
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  geom_text_repel(
    aes(label = Label),
    size = 4,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.4
  ) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot of Differentially Expressed miRNAs",
    subtitle = "Highlighted: 2 upregulated and 2 downregulated miRNAs",
    x = "logFC",
    y = "-log10(P.Value)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )

print(p_volcano_4)
# Combine selected miRNAs
selected_targets <- c(up_targets, down_targets)

volcano_selected4 <- deg %>%
  mutate(
    minus_log10_p = -log10(P.Value),
    Category = case_when(
      miRNA %in% up_targets ~ "Upregulated",
      miRNA %in% down_targets ~ "Downregulated",
      TRUE ~ "Other"
    ),
    Label = ifelse(miRNA %in% selected_targets, miRNA, "")
  )

p_volcano_4 <- ggplot(volcano_selected4, aes(x = logFC, y = minus_log10_p)) +
  
  # Other miRNAs in gray
  geom_point(
    data = subset(volcano_selected4, Category == "Other"),
    color = "gray70",
    size = 2,
    alpha = 0.6
  ) +
  
  # Downregulated in blue
  geom_point(
    data = subset(volcano_selected4, Category == "Downregulated"),
    color = "blue",
    size = 4,
    alpha = 1
  ) +
  
  # Upregulated in red
  geom_point(
    data = subset(volcano_selected4, Category == "Upregulated"),
    color = "red",
    size = 4,
    alpha = 1
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  geom_text_repel(
    aes(label = Label),
    size = 4,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.4
  ) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot of Differentially Expressed miRNAs",
    subtitle = "Highlighted: 2 upregulated and 2 downregulated miRNAs",
    x = "logFC",
    y = "-log10(P.Value)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )

print(p_volcano_4)
# ============================================================
# Create folders for saving results
# ============================================================

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

getwd()
list.files()
# ============================================================
# Save important tables
# ============================================================

deg$miRNA <- rownames(deg)

# Full differential expression table
write.csv(
  deg,
  "results/01_all_differential_expression_results.csv",
  row.names = FALSE
)

# Upregulated miRNAs based on unadjusted P.Value
upregulated_pvalue <- deg %>%
  filter(logFC > 0, P.Value < 0.05) %>%
  arrange(P.Value)

write.csv(
  upregulated_pvalue,
  "results/02_upregulated_miRNAs_Pvalue_only.csv",
  row.names = FALSE
)

# Downregulated miRNAs based on unadjusted P.Value
downregulated_pvalue <- deg %>%
  filter(logFC < 0, P.Value < 0.05) %>%
  arrange(P.Value)

write.csv(
  downregulated_pvalue,
  "results/03_downregulated_miRNAs_Pvalue_only.csv",
  row.names = FALSE
)

# Your selected upregulated miRNAs
up_targets <- c("hsa-miR-590-5p", "hsa-miR-182-5p")

# Top 2 downregulated miRNAs
down_targets <- deg %>%
  filter(logFC < 0, P.Value < 0.05) %>%
  arrange(logFC) %>%
  slice(1:2) %>%
  pull(miRNA)

selected_targets <- c(up_targets, down_targets)

# Selected 4 miRNAs table
selected_table <- deg %>%
  filter(miRNA %in% selected_targets) %>%
  select(miRNA, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write.csv(
  selected_table,
  "results/04_selected_2_upregulated_2_downregulated_miRNAs.csv",
  row.names = FALSE
)

selected_table
# ============================================================
# Save volcano plot
# ============================================================

ggsave(
  filename = "plots/01_volcano_selected_2_up_2_down.png",
  plot = p_volcano_4,
  width = 8,
  height = 6,
  dpi = 300
)
# ============================================================
# Save boxplots for selected upregulated miRNAs
# ============================================================

my_colors <- c(
  "Normal" = "blue",
  "TNBC_Cancer" = "red"
)

save_miRNA_boxplot <- function(mirna_name) {
  
  hit <- grep(mirna_name, rownames(expr_tnbc), ignore.case = TRUE, value = TRUE)
  
  if (length(hit) == 0) {
    message("Not found: ", mirna_name)
    return(NULL)
  }
  
  hit <- hit[1]
  
  df <- data.frame(
    Sample = colnames(expr_tnbc),
    Expression = as.numeric(expr_tnbc[hit, ]),
    Group = group,
    miRNA = hit
  )
  
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Group), width = 0.12, size = 3) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0(hit, " expression"),
      subtitle = "GSE154255: TNBC cancer vs adjacent normal",
      x = "",
      y = "Expression level"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p)
  
  safe_name <- gsub("-", "_", hit)
  
  ggsave(
    filename = paste0("plots/03_boxplot_", safe_name, ".png"),
    plot = p,
    width = 5,
    height = 5,
    dpi = 300
  )
  
  write.csv(
    df,
    paste0("results/05_expression_values_", safe_name, ".csv"),
    row.names = FALSE
  )
  
  return(p)
}

p_box_590 <- save_miRNA_boxplot("hsa-miR-590-5p")
p_box_182 <- save_miRNA_boxplot("hsa-miR-182-5p")
# ============================================================
# Save heatmap
# ============================================================

top20_miRNAs <- deg %>%
  filter(P.Value < 0.05) %>%
  arrange(P.Value) %>%
  slice(1:20)

top20_names <- top20_miRNAs$miRNA
heat_expr <- expr_tnbc[top20_names, ]

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(heat_expr)

ann_colors <- list(
  Group = my_colors
)

heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

png(
  filename = "plots/04_heatmap_top20_miRNAs.png",
  width = 1800,
  height = 1600,
  res = 300
)

pheatmap(
  heat_expr,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heat_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  main = "Top 20 differentially expressed miRNAs: TNBC vs Normal"
)

dev.off()

write.csv(
  top20_miRNAs,
  "results/06_top20_miRNAs_for_heatmap.csv",
  row.names = FALSE
)
# ============================================================
# Check saved files
# ============================================================

list.files("results")
list.files("plots")
# ============================================================
# Save ggplot figures in multiple formats
# PNG, PDF, TIFF, SVG
# ============================================================

dir.create("plots", showWarnings = FALSE)

save_plot_all_formats <- function(plot_object, file_name, width = 7, height = 6) {
  
  ggsave(
    filename = paste0("plots/", file_name, ".png"),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
  
  ggsave(
    filename = paste0("plots/", file_name, ".pdf"),
    plot = plot_object,
    width = width,
    height = height
  )
  
  ggsave(
    filename = paste0("plots/", file_name, ".tiff"),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    compression = "lzw"
  )
  
  ggsave(
    filename = paste0("plots/", file_name, ".svg"),
    plot = plot_object,
    width = width,
    height = height
  )
}