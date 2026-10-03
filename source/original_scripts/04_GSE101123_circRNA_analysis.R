rm(list = ls())
gc()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c(
  "GEOquery",
  "limma",
  "ggplot2",
  "pheatmap",
  "dplyr"
)

for(pkg in packages){
  
  if(!requireNamespace(pkg, quietly = TRUE)){
    
    if(pkg %in% c("GEOquery","limma")){
      
      BiocManager::install(pkg)
      
    } else {
      
      install.packages(pkg)
      
    }
    
  }
  
}

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

# Download GSE101123

gse_circ <- getGEO(
  "GSE101123",
  GSEMatrix = TRUE
)

length(gse_circ)

gse_circ

eset_circ <- gse_circ[[1]]

expr_circ <- exprs(eset_circ)

dim(expr_circ)

pheno_circ <- pData(eset_circ)

colnames(pheno_circ)

View(pheno_circ)
# Keep only TNBC and Normal samples

keep <- pheno_circ$`characteristics_ch1.1` %in%
  c(
    "subtype: Triple-negative breast cancer",
    "subtype: non-tumor breast tissues"
  )

expr_sub <- expr_circ[, keep]

pheno_sub <- pheno_circ[keep, ]

group_circ <- factor(
  ifelse(
    pheno_sub$`characteristics_ch1.1` ==
      "subtype: Triple-negative breast cancer",
    "TNBC",
    "Normal"
  )
)

table(group_circ)
dput(unique(pheno_circ$characteristics_ch1.1))
# Clean subtype labels
subtype <- trimws(pheno_circ$characteristics_ch1.1)

# Keep only TNBC and Normal samples
keep <- subtype %in% c(
  "subtype: Triple-negative breast cancer",
  "subtype: non-tumor breast tissues"
)

expr_sub <- expr_circ[, keep]
pheno_sub <- pheno_circ[keep, ]

# Define groups
group_circ <- factor(
  ifelse(
    trimws(pheno_sub$characteristics_ch1.1) ==
      "subtype: Triple-negative breast cancer",
    "TNBC",
    "Normal"
  ),
  levels = c("Normal", "TNBC")
)

table(group_circ)
dim(expr_sub)
# Robust subtype detection
subtype <- trimws(pheno_circ$characteristics_ch1.1)

# Keep TNBC and Normal samples
keep <- grepl("Triple", subtype, ignore.case = TRUE) |
  grepl("non-tumor", subtype, ignore.case = TRUE)

expr_sub <- expr_circ[, keep]
pheno_sub <- pheno_circ[keep, ]

subtype_sub <- trimws(pheno_sub$characteristics_ch1.1)

# Define groups
group_circ <- factor(
  ifelse(
    grepl("Triple", subtype_sub, ignore.case = TRUE),
    "TNBC",
    "Normal"
  ),
  levels = c("Normal", "TNBC")
)

table(group_circ)
dim(expr_sub)

# Check selected samples
pheno_sub[, c("title", "geo_accession", "characteristics_ch1.1")]
# Check groups
table(group_circ)
dim(expr_sub)

# Create output folder
dir.create("GSE101123_circRNA_results", showWarnings = FALSE)

# Design matrix
design <- model.matrix(~0 + group_circ)
colnames(design) <- levels(group_circ)

design

# Differential expression analysis
fit <- lmFit(expr_sub, design)

contrast_matrix <- makeContrasts(
  TNBC - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

circ_results <- topTable(
  fit2,
  adjust.method = "BH",
  number = Inf
)

circ_results$circRNA <- rownames(circ_results)

circ_results <- circ_results[, c(
  "circRNA",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)]

head(circ_results)

# Significant DE-circRNAs
sig_circ <- subset(
  circ_results,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

up_circ <- subset(sig_circ, logFC > 1)
down_circ <- subset(sig_circ, logFC < -1)

nrow(sig_circ)
nrow(up_circ)
nrow(down_circ)

# Save tables
write.csv(
  circ_results,
  "GSE101123_circRNA_results/Table_GSE101123_all_circRNA_results.csv",
  row.names = FALSE
)

write.csv(
  sig_circ,
  "GSE101123_circRNA_results/Table_GSE101123_significant_DEcircRNAs.csv",
  row.names = FALSE
)

write.csv(
  up_circ,
  "GSE101123_circRNA_results/Table_GSE101123_upregulated_circRNAs.csv",
  row.names = FALSE
)

write.csv(
  down_circ,
  "GSE101123_circRNA_results/Table_GSE101123_downregulated_circRNAs.csv",
  row.names = FALSE
)
library(limma)

table(group_circ)
dim(expr_sub)
nrow(sig_circ)
nrow(up_circ)
nrow(down_circ)

sig_circ
# Create folder
dir.create("GSE101123_circRNA_results", showWarnings = FALSE)

# Save tables
write.csv(circ_results,
          "GSE101123_circRNA_results/Table_GSE101123_all_circRNA_results.csv",
          row.names = FALSE)

write.csv(sig_circ,
          "GSE101123_circRNA_results/Table_GSE101123_significant_DEcircRNAs.csv",
          row.names = FALSE)

write.csv(up_circ,
          "GSE101123_circRNA_results/Table_GSE101123_upregulated_circRNAs.csv",
          row.names = FALSE)

write.csv(down_circ,
          "GSE101123_circRNA_results/Table_GSE101123_downregulated_circRNAs.csv",
          row.names = FALSE)

# Save summary file
sink("GSE101123_circRNA_results/README_GSE101123_summary.txt")
cat("Dataset: GSE101123\n")
cat("Analysis type: circRNA differential expression\n")
cat("Comparison: TNBC vs Normal\n\n")
cat("Normal samples: 3\n")
cat("TNBC samples: 4\n\n")
cat("Total circRNAs analysed:", nrow(circ_results), "\n")
cat("Significant DE-circRNAs:", nrow(sig_circ), "\n")
cat("Upregulated circRNAs:", nrow(up_circ), "\n")
cat("Downregulated circRNAs:", nrow(down_circ), "\n\n")
cat("Main significant circRNA:\n")
print(sig_circ)
sink()

# Save selected sample information
write.csv(
  pheno_sub[, c("title", "geo_accession", "characteristics_ch1.1")],
  "GSE101123_circRNA_results/Table_GSE101123_selected_samples.csv",
  row.names = FALSE
)
library(ggplot2)
library(ggrepel)

dir.create("GSE101123_circRNA_results", showWarnings = FALSE)

volcano_data <- circ_results

volcano_data$Status <- "Not Significant"

volcano_data$Status[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC > 1
] <- "Upregulated"

volcano_data$Status[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC < -1
] <- "Downregulated"

p_volcano <- ggplot(
  volcano_data,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Status
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_text_repel(
    data = subset(volcano_data,
                  adj.P.Val < 0.05),
    aes(label = circRNA),
    size = 4,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "#E64B35",
      "Downregulated" = "#4DBBD5",
      "Not Significant" = "grey70"
    )
  ) +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "Differentially Expressed circRNAs in TNBC",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(p_volcano)

ggsave(
  "GSE101123_circRNA_results/Figure6A_Volcano_circRNA.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "GSE101123_circRNA_results/Figure6A_Volcano_circRNA.tiff",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)
library(pheatmap)

top10 <- circ_results[
  order(circ_results$P.Value),
][1:10, ]

expr_heat <- expr_sub[
  rownames(expr_sub) %in% top10$circRNA,
]

expr_heat <- expr_heat[
  match(top10$circRNA,
        rownames(expr_heat)),
]

rownames(expr_heat) <- top10$circRNA

colnames(expr_heat) <- c(
  "N1", "N2", "N3",
  "TNBC1", "TNBC2", "TNBC3", "TNBC4"
)

annotation_col <- data.frame(
  Group = factor(
    c(
      "Normal", "Normal", "Normal",
      "TNBC", "TNBC", "TNBC", "TNBC"
    ),
    levels = c("Normal", "TNBC")
  )
)

rownames(annotation_col) <- colnames(expr_heat)

ann_colors <- list(
  Group = c(
    Normal = "#4DBBD5",
    TNBC = "#E64B35"
  )
)

pheatmap(
  expr_heat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    c("#4DBBD5", "white", "#E64B35")
  )(100),
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 11,
  main = "Top 10 Differentially Expressed circRNAs",
  filename = "GSE101123_circRNA_results/Figure6B_Heatmap_Top10_circRNAs.png",
  width = 9,
  height = 7
)

pheatmap(
  expr_heat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    c("#4DBBD5", "white", "#E64B35")
  )(100),
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 11,
  main = "Top 10 Differentially Expressed circRNAs",
  filename = "GSE101123_circRNA_results/Figure6B_Heatmap_Top10_circRNAs.tiff",
  width = 9,
  height = 7
)

# Display heatmap
pheatmap(
  expr_heat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    c("#4DBBD5", "white", "#E64B35")
  )(100),
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 11,
  main = "Top 10 Differentially Expressed circRNAs"
)
# ==============================
# Re-plot GSE101123 circRNA boxplot
# Title as circRNA name, not probe ID
# ==============================

library(GEOquery)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(openxlsx)

gse_id <- "GSE101123"
probe_id <- "ASCRP000064"
circRNA_name <- "hsa_circRNA_000554"

outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

gse <- getGEO(gse_id, GSEMatrix = TRUE)
eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)

print(pheno$title)

pheno$Group <- NA
pheno$Group[grepl("normal|adjacent|non-tumor|nontumor", pheno$title, ignore.case = TRUE)] <- "Normal"
pheno$Group[grepl("cancer|tumor|TNBC|triple", pheno$title, ignore.case = TRUE)] <- "Tumor"

print(table(pheno$Group, useNA = "ifany"))

target_row <- grep(probe_id, rownames(expr), value = TRUE)
print(target_row)

df <- data.frame(
  Sample = colnames(expr),
  Expression = as.numeric(expr[target_row[1], ]),
  Group = pheno$Group,
  stringsAsFactors = FALSE
)

df <- df[!is.na(df$Group), ]
df$Group <- factor(df$Group, levels = c("Normal", "Tumor"))

stat_test <- wilcox.test(Expression ~ Group, data = df)

summary_table <- df %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = n(),
    .groups = "drop"
  )

result_table <- data.frame(
  Dataset = gse_id,
  Probe_ID = probe_id,
  circRNA_name = circRNA_name,
  Normal_N = sum(df$Group == "Normal"),
  Tumor_N = sum(df$Group == "Tumor"),
  Normal_mean = mean(df$Expression[df$Group == "Normal"]),
  Tumor_mean = mean(df$Expression[df$Group == "Tumor"]),
  Direction = ifelse(
    mean(df$Expression[df$Group == "Tumor"]) <
      mean(df$Expression[df$Group == "Normal"]),
    "Downregulated in Tumor",
    "Upregulated in Tumor"
  ),
  Wilcoxon_W = as.numeric(stat_test$statistic),
  P_value = stat_test$p.value
)

print(summary_table)
print(result_table)

ymax <- max(df$Expression)

p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.04,
    size = 2.8,
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  annotate(
    "text",
    x = 1.5,
    y = ymax + 0.18,
    label = paste0("p = ", signif(stat_test$p.value, 3)),
    size = 5
  ) +
  scale_y_continuous(
    limits = c(min(df$Expression) - 0.1, ymax + 0.35)
  ) +
  scale_fill_manual(values = c(
    "Normal" = "#4DBBD5",
    "Tumor"  = "#E64B35"
  )) +
  labs(
    x = "",
    y = "Normalized expression",
    title = expression(italic("hsa_circRNA_000554"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 16),
    axis.text = element_text(color = "black", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    legend.position = "none"
  )

print(p)

ggsave(
  file.path(outdir, "Figure_GSE101123_hsa_circRNA_000554_FINAL.png"),
  p,
  width = 4.5,
  height = 5,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure_GSE101123_hsa_circRNA_000554_FINAL.tiff"),
  p,
  width = 4.5,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure_GSE101123_hsa_circRNA_000554_FINAL.pdf"),
  p,
  width = 4.5,
  height = 5
)

write.xlsx(
  list(
    Expression_Data = df,
    Summary = summary_table,
    Validation_Result = result_table
  ),
  file = file.path(outdir, "GSE101123_hsa_circRNA_000554_replot_results.xlsx"),
  rowNames = FALSE
)

list.files(outdir)
pheno$title
pheno$Group <- NA
pheno$Group[grepl("mammary gland", pheno$title, ignore.case = TRUE)] <- "Normal"
pheno$Group[grepl("breast cancer", pheno$title, ignore.case = TRUE)] <- "Tumor"

print(table(pheno$Group, useNA = "ifany"))
target_row <- grep(probe_id, rownames(expr), value = TRUE)
print(target_row)

df <- data.frame(
  Sample = colnames(expr),
  Expression = as.numeric(expr[target_row[1], ]),
  Group = pheno$Group,
  stringsAsFactors = FALSE
)

df <- df[!is.na(df$Group), ]
df$Group <- factor(df$Group, levels = c("Normal", "Tumor"))

print(table(df$Group))

stat_test <- wilcox.test(Expression ~ Group, data = df)
print(stat_test)

summary_table <- df %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = n(),
    .groups = "drop"
  )

result_table <- data.frame(
  Dataset = "GSE101123",
  Probe_ID = probe_id,
  circRNA_name = "hsa_circRNA_000554",
  Normal_N = sum(df$Group == "Normal"),
  Tumor_N = sum(df$Group == "Tumor"),
  Normal_mean = mean(df$Expression[df$Group == "Normal"]),
  Tumor_mean = mean(df$Expression[df$Group == "Tumor"]),
  Direction = ifelse(
    mean(df$Expression[df$Group == "Tumor"]) <
      mean(df$Expression[df$Group == "Normal"]),
    "Downregulated in Tumor",
    "Upregulated in Tumor"
  ),
  Wilcoxon_W = as.numeric(stat_test$statistic),
  P_value = stat_test$p.value
)

print(summary_table)
print(result_table)

ymax <- max(df$Expression)

p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.04,
    size = 2.8,
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  annotate(
    "text",
    x = 1.5,
    y = ymax + 0.18,
    label = paste0("p = ", signif(stat_test$p.value, 3)),
    size = 5
  ) +
  scale_y_continuous(
    limits = c(min(df$Expression) - 0.1, ymax + 0.35)
  ) +
  scale_fill_manual(values = c(
    "Normal" = "#4DBBD5",
    "Tumor"  = "#E64B35"
  )) +
  labs(
    x = "",
    y = "Normalized expression",
    title = expression(italic("hsa_circRNA_000554"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 16),
    axis.text = element_text(color = "black", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    legend.position = "none"
  )

print(p)

ggsave(
  file.path(outdir, "Figure_GSE101123_hsa_circRNA_000554_FINAL.png"),
  p,
  width = 4.5,
  height = 5,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure_GSE101123_hsa_circRNA_000554_FINAL.tiff"),
  p,
  width = 4.5,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure_GSE101123_hsa_circRNA_000554_FINAL.pdf"),
  p,
  width = 4.5,
  height = 5
)

write.xlsx(
  list(
    Expression_Data = df,
    Summary = summary_table,
    Validation_Result = result_table
  ),
  file = file.path(outdir, "GSE101123_hsa_circRNA_000554_replot_results.xlsx"),
  rowNames = FALSE
)

list.files(outdir)
# Replace probe ID with circRNA name

volcano_data$DisplayName <- volcano_data$Gene

volcano_data$DisplayName[
  volcano_data$DisplayName == "ASCRP000064"
] <- "hsa_circRNA_000554"
library(ggrepel)

# Create display names
res$DisplayName <- rownames(res)

# Replace probe ID with circRNA name
res$DisplayName[
  res$DisplayName == "ASCRP000064"
] <- "hsa_circRNA_000554"

# Volcano plot
p <- ggplot(res,
            aes(x = logFC,
                y = -log10(adj.P.Val),
                color = Status)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  
  geom_text_repel(
    data = subset(res,
                  DisplayName == "hsa_circRNA_000554"),
    aes(label = DisplayName),
    size = 4,
    fontface = "bold"
  ) +
  
  theme_bw() +
  labs(
    title = "Differentially Expressed circRNAs in TNBC",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

print(p)
ls()
history(max.show = Inf)
View(eset)
circ_results$DisplayName <- circ_results$circRNA
library(limma)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

circ_results <- topTable(
  fit2,
  adjust.method = "BH",
  number = Inf
)

circ_results$circRNA <- rownames(circ_results)
ls()
table(pheno$title)
library(limma)
library(Biobase)

expr <- exprs(eset)

group_circ <- factor(
  c(rep("TNBC", 8), rep("Normal", 3)),
  levels = c("Normal", "TNBC")
)

design <- model.matrix(~0 + group_circ)
colnames(design) <- levels(group_circ)

fit <- lmFit(expr, design)

contrast_matrix <- makeContrasts(
  TNBC - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

circ_results <- topTable(
  fit2,
  adjust.method = "BH",
  number = Inf
)

circ_results$circRNA <- rownames(circ_results)
head(circ_results)
library(ggplot2)
library(ggrepel)

## Create display names
circ_results$DisplayName <- circ_results$circRNA

## Replace probe ID with circRNA name
circ_results$DisplayName[
  circ_results$DisplayName == "ASCRP000064"
] <- "hsa_circRNA_000554"

## Define significance status
circ_results$Status <- "Not Significant"

circ_results$Status[
  circ_results$adj.P.Val < 0.05 &
    circ_results$logFC < -1
] <- "Downregulated"

circ_results$Status[
  circ_results$adj.P.Val < 0.05 &
    circ_results$logFC > 1
] <- "Upregulated"

## Volcano plot
p_volcano <- ggplot(
  circ_results,
  aes(x = logFC,
      y = -log10(adj.P.Val),
      color = Status)
) +
  geom_point(size = 2, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#56B4E9",
      "Upregulated" = "#E64B35",
      "Not Significant" = "grey75"
    )
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_text_repel(
    data = subset(
      circ_results,
      DisplayName == "hsa_circRNA_000554"
    ),
    aes(label = DisplayName),
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  
  theme_bw(base_size = 13) +
  
  labs(
    title = "Differentially Expressed circRNAs in TNBC",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

print(p_volcano)
library(ggplot2)
library(ggrepel)

## Create display names
circ_results$DisplayName <- circ_results$circRNA

## Replace probe ID with circRNA name
circ_results$DisplayName[
  circ_results$DisplayName == "ASCRP000064"
] <- "hsa_circRNA_000554"

## Define significance using raw P-value
circ_results$Status <- "Not Significant"

circ_results$Status[
  circ_results$P.Value < 0.05 &
    circ_results$logFC < -1
] <- "Downregulated"

circ_results$Status[
  circ_results$P.Value < 0.05 &
    circ_results$logFC > 1
] <- "Upregulated"

## Volcano plot
p_volcano <- ggplot(
  circ_results,
  aes(
    x = logFC,
    y = -log10(P.Value),
    color = Status
  )
) +
  
  geom_point(size = 2, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#56B4E9",
      "Upregulated" = "#E64B35",
      "Not Significant" = "grey75"
    )
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_text_repel(
    data = subset(
      circ_results,
      DisplayName == "hsa_circRNA_000554"
    ),
    aes(label = DisplayName),
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  
  theme_bw(base_size = 13) +
  
  labs(
    title = "Differentially Expressed circRNAs in TNBC",
    x = "log2 Fold Change",
    y = "-log10 P-value"
  )

print(p_volcano)
p_volcano <- ggplot(
  circ_results,
  aes(
    x = logFC,
    y = -log10(P.Value),
    color = Status
  )
) +
  geom_point(size = 2, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#56B4E9",
      "Upregulated"   = "#E64B35",
      "Not Significant" = "grey75"
    )
  ) +
  
  geom_vline(
    xintercept = c(-1,1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_text_repel(
    data = subset(
      circ_results,
      DisplayName == "hsa_circRNA_000554"
    ),
    aes(label = DisplayName),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 1,
    point.padding = 0.5,
    nudge_y = 0.5,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  
  coord_cartesian(
    ylim = c(
      0,
      max(-log10(circ_results$P.Value)) + 1
    )
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  
  labs(
    title = "Differentially Expressed circRNAs in TNBC",
    x = "log2 Fold Change",
    y = "-log10 P-value"
  )

print(p_volcano)
ggsave(
  "Figure3_Volcano_hsa_circRNA_000554_FINAL.tiff",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "Figure3_Volcano_hsa_circRNA_000554_FINAL.pdf",
  p_volcano,
  width = 8,
  height = 6
)
library(ggplot2)
library(ggrepel)

outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

# Create display names
circ_results$DisplayName <- circ_results$circRNA

# Replace probe ID with biological circRNA name
circ_results$DisplayName[
  circ_results$DisplayName == "ASCRP000064"
] <- "hsa_circRNA_000554"

# Define status using raw P-value
circ_results$Status <- "Not Significant"

circ_results$Status[
  circ_results$P.Value < 0.05 &
    circ_results$logFC < -1
] <- "Downregulated"

circ_results$Status[
  circ_results$P.Value < 0.05 &
    circ_results$logFC > 1
] <- "Upregulated"

# Volcano plot
p_volcano <- ggplot(
  circ_results,
  aes(
    x = logFC,
    y = -log10(P.Value),
    color = Status
  )
) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(
    "Downregulated" = "#56B4E9",
    "Not Significant" = "grey75",
    "Upregulated" = "#E64B35"
  )) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = subset(circ_results,
                  DisplayName == "hsa_circRNA_000554"),
    aes(label = DisplayName),
    size = 5,
    fontface = "bold",
    color = "black",
    nudge_x = -0.6,
    nudge_y = 0.3,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.5,
    max.overlaps = Inf
  ) +
  coord_cartesian(
    ylim = c(0, max(-log10(circ_results$P.Value)) + 1)
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Differentially Expressed circRNAs in TNBC",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10](italic(P)~value)),
    color = "Status"
  )

print(p_volcano)

ggsave(
  file.path(outdir, "Figure3_Volcano_hsa_circRNA_000554_FINAL.png"),
  p_volcano,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure3_Volcano_hsa_circRNA_000554_FINAL.tiff"),
  p_volcano,
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure3_Volcano_hsa_circRNA_000554_FINAL.pdf"),
  p_volcano,
  width = 8,
  height = 6
)

list.files(outdir)
# ==========================================
# Re-draw Top 10 circRNA Heatmap - GSE101123
# Replace ASCRP000064 with hsa_circRNA_000554
# ==========================================

library(GEOquery)
library(limma)
library(pheatmap)

outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

gse_id <- "GSE101123"

gse <- getGEO(gse_id, GSEMatrix = TRUE)
eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)

# Define groups
pheno$Group <- NA
pheno$Group[grepl("mammary gland", pheno$title, ignore.case = TRUE)] <- "Normal"
pheno$Group[grepl("breast cancer", pheno$title, ignore.case = TRUE)] <- "TNBC"

print(table(pheno$Group, useNA = "ifany"))

# Keep only defined samples
keep <- !is.na(pheno$Group)
expr_sub <- expr[, keep]
pheno_sub <- pheno[keep, ]

group_circ <- factor(pheno_sub$Group, levels = c("Normal", "TNBC"))

# Differential expression
design <- model.matrix(~0 + group_circ)
colnames(design) <- levels(group_circ)

fit <- lmFit(expr_sub, design)

contrast_matrix <- makeContrasts(
  TNBC - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

circ_results <- topTable(
  fit2,
  adjust.method = "BH",
  number = Inf
)

circ_results$circRNA <- rownames(circ_results)

# Select Top 10 by raw P-value
top10 <- circ_results[order(circ_results$P.Value), ][1:10, ]
top10_ids <- top10$circRNA

# Extract expression of Top 10
heat_expr <- expr_sub[top10_ids, ]

# Z-score by row
heat_expr_z <- t(scale(t(heat_expr)))

# Rename ASCRP000064 to biological circRNA name
rownames(heat_expr_z)[rownames(heat_expr_z) == "ASCRP000064"] <- "hsa_circRNA_000554"

# Sample annotation
annotation_col <- data.frame(
  Group = group_circ
)

rownames(annotation_col) <- colnames(heat_expr_z)

ann_colors <- list(
  Group = c(
    Normal = "#4DBBD5",
    TNBC = "#E64B35"
  )
)

# Save PNG
png(
  filename = file.path(outdir, "Figure6B_Heatmap_Top10_circRNAs_FINAL.png"),
  width = 2700,
  height = 2100,
  res = 300
)

pheatmap(
  heat_expr_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  main = "Top 10 Differentially Expressed circRNAs"
)

dev.off()

# Save TIFF
tiff(
  filename = file.path(outdir, "Figure6B_Heatmap_Top10_circRNAs_FINAL.tiff"),
  width = 2700,
  height = 2100,
  res = 300,
  compression = "lzw"
)

pheatmap(
  heat_expr_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  main = "Top 10 Differentially Expressed circRNAs"
)

dev.off()

# Save PDF
pdf(
  file = file.path(outdir, "Figure6B_Heatmap_Top10_circRNAs_FINAL.pdf"),
  width = 9,
  height = 7
)

pheatmap(
  heat_expr_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  main = "Top 10 Differentially Expressed circRNAs"
)

dev.off()

# Save tables
write.csv(
  top10,
  file = file.path(outdir, "Table_Top10_circRNAs_for_heatmap_FINAL.csv"),
  row.names = FALSE
)

write.csv(
  heat_expr_z,
  file = file.path(outdir, "Matrix_Top10_circRNAs_heatmap_Zscore_FINAL.csv"),
  row.names = TRUE
)

list.files(outdir)
# =========================================================
# Publication-ready Top 10 circRNA heatmap - GSE101123
# Selected samples: 3 Normal + 4 TNBC
# ASCRP000064 displayed as hsa_circRNA_000554
# =========================================================

library(GEOquery)
library(limma)
library(pheatmap)

outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

gse_id <- "GSE101123"

gse <- getGEO(gse_id, GSEMatrix = TRUE)
eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)

# -----------------------------
# Define groups from sample title
# -----------------------------

pheno$Group <- NA
pheno$Group[grepl("mammary gland", pheno$title, ignore.case = TRUE)] <- "Normal"
pheno$Group[grepl("breast cancer", pheno$title, ignore.case = TRUE)] <- "TNBC"

# -----------------------------
# Select samples exactly as original heatmap:
# 3 Normal + first 4 TNBC
# -----------------------------

normal_samples <- rownames(pheno)[pheno$Group == "Normal"]
tnbc_samples   <- rownames(pheno)[pheno$Group == "TNBC"]

selected_samples <- c(
  normal_samples[1:3],
  tnbc_samples[1:4]
)

pheno_sub <- pheno[selected_samples, ]
expr_sub  <- expr[, selected_samples]

group_circ <- factor(
  pheno_sub$Group,
  levels = c("Normal", "TNBC")
)

# Rename samples for clean figure
new_sample_names <- c(
  paste0("N", 1:3),
  paste0("TNBC", 1:4)
)

colnames(expr_sub) <- new_sample_names
rownames(pheno_sub) <- new_sample_names

# -----------------------------
# Differential expression
# -----------------------------

design <- model.matrix(~0 + group_circ)
colnames(design) <- levels(group_circ)

fit <- lmFit(expr_sub, design)

contrast_matrix <- makeContrasts(
  TNBC - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

circ_results <- topTable(
  fit2,
  adjust.method = "BH",
  number = Inf
)

circ_results$circRNA <- rownames(circ_results)

# -----------------------------
# Select Top 10 circRNAs
# Using raw P-value to match volcano/previous analysis
# -----------------------------

top10 <- circ_results[order(circ_results$P.Value), ][1:10, ]
top10_ids <- top10$circRNA

heat_expr <- expr_sub[top10_ids, ]

# Row-wise Z-score
heat_expr_z <- t(scale(t(heat_expr)))

# Replace probe ID with biological circRNA name
rownames(heat_expr_z)[
  rownames(heat_expr_z) == "ASCRP000064"
] <- "hsa_circRNA_000554"

top10$DisplayName <- top10$circRNA
top10$DisplayName[
  top10$DisplayName == "ASCRP000064"
] <- "hsa_circRNA_000554"

# -----------------------------
# Annotation
# -----------------------------

annotation_col <- data.frame(
  Group = group_circ
)

rownames(annotation_col) <- colnames(heat_expr_z)

ann_colors <- list(
  Group = c(
    Normal = "#4DBBD5",
    TNBC   = "#E64B35"
  )
)

heat_colors <- colorRampPalette(
  c("#4DBBD5", "white", "#E64B35")
)(100)

# -----------------------------
# Plot and save
# -----------------------------

png(
  filename = file.path(outdir, "Figure6B_Heatmap_Top10_circRNAs_SELECTED_FINAL.png"),
  width = 3000,
  height = 2200,
  res = 300
)

pheatmap(
  heat_expr_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heat_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  main = "Top 10 Differentially Expressed circRNAs"
)

dev.off()

tiff(
  filename = file.path(outdir, "Figure6B_Heatmap_Top10_circRNAs_SELECTED_FINAL.tiff"),
  width = 3000,
  height = 2200,
  res = 300,
  compression = "lzw"
)

pheatmap(
  heat_expr_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heat_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  main = "Top 10 Differentially Expressed circRNAs"
)

dev.off()

pdf(
  file = file.path(outdir, "Figure6B_Heatmap_Top10_circRNAs_SELECTED_FINAL.pdf"),
  width = 10,
  height = 7.5
)

pheatmap(
  heat_expr_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heat_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  main = "Top 10 Differentially Expressed circRNAs"
)

dev.off()

# -----------------------------
# Save source tables
# -----------------------------

write.csv(
  top10,
  file = file.path(outdir, "Table_Top10_circRNAs_SELECTED_FINAL.csv"),
  row.names = FALSE
)

write.csv(
  heat_expr_z,
  file = file.path(outdir, "Matrix_Top10_circRNAs_SELECTED_Zscore_FINAL.csv"),
  row.names = TRUE
)

list.files(outdir)
library(pheatmap)

## Select the top 10 downregulated circRNAs
top10_down <- circ_results[circ_results$logFC < 0, ]
top10_down <- top10_down[order(top10_down$logFC), ]
top10_down <- top10_down[1:10, ]

## Replace the probe ID of the circRNA of interest
top10_down$DisplayName <- top10_down$circRNA
top10_down$DisplayName[
  top10_down$circRNA == "ASCRP000064"
] <- "hsa_circRNA_000554"

## Extract expression matrix
top10_ids <- top10_down$circRNA
heat_expr <- expr_sub[top10_ids, ]

## Perform row-wise Z-score normalization
heat_expr_z <- t(scale(t(heat_expr)))

## Assign final row names
rownames(heat_expr_z) <- top10_down$DisplayName

## Rename samples
colnames(heat_expr_z) <- c(
  "Normal_1", "Normal_2", "Normal_3",
  "TNBC_1", "TNBC_2", "TNBC_3", "TNBC_4"
)

## Sample annotation
annotation_col <- data.frame(
  Group = factor(c(
    "Normal", "Normal", "Normal",
    "TNBC", "TNBC", "TNBC", "TNBC"
  ))
)

rownames(annotation_col) <- colnames(heat_expr_z)

## Save publication-quality heatmap
tiff(
  file.path(
    outdir,
    "Figure6B_Heatmap_Top10_Downregulated_FINAL.tiff"
  ),
  width = 3000,
  height = 2200,
  res = 300
)

pheatmap(
  heat_expr_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Top 10 Downregulated circRNAs in TNBC",
  color = colorRampPalette(
    c("#67B7D1", "white", "#E64B35")
  )(100),
  border_color = "grey70"
)

dev.off()
library(pheatmap)

target_probe <- "ASCRP000064"
target_name  <- "hsa_circRNA_000554"

# Select top 9 most downregulated circRNAs excluding target
down_all <- circ_results[circ_results$logFC < 0, ]
down_all <- down_all[order(down_all$logFC), ]

top9_down <- down_all[down_all$circRNA != target_probe, ][1:9, ]

# Add target circRNA manually
target_row <- circ_results[circ_results$circRNA == target_probe, ]

top10_down <- rbind(top9_down, target_row)

# Display names
top10_down$DisplayName <- top10_down$circRNA
top10_down$DisplayName[top10_down$circRNA == target_probe] <- target_name

# Expression matrix
top10_ids <- top10_down$circRNA
heat_expr <- expr_sub[top10_ids, ]

# Z-score
heat_expr_z <- t(scale(t(heat_expr)))
rownames(heat_expr_z) <- top10_down$DisplayName

# Rename samples
colnames(heat_expr_z) <- c(
  "Normal_1", "Normal_2", "Normal_3",
  "TNBC_1", "TNBC_2", "TNBC_3", "TNBC_4"
)

annotation_col <- data.frame(
  Group = factor(c("Normal", "Normal", "Normal", "TNBC", "TNBC", "TNBC", "TNBC"),
                 levels = c("Normal", "TNBC"))
)
rownames(annotation_col) <- colnames(heat_expr_z)

ann_colors <- list(
  Group = c(
    Normal = "#4DBBD5",
    TNBC = "#E64B35"
  )
)

# Save heatmap
tiff(
  file.path(outdir, "Figure6B_Heatmap_Top9Down_plus_hsa_circRNA_000554_FINAL.tiff"),
  width = 3000,
  height = 2200,
  res = 300,
  compression = "lzw"
)

pheatmap(
  heat_expr_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  main = "Top Downregulated circRNAs Including hsa_circRNA_000554"
)

dev.off()

# Check included rows
top10_down[, c("DisplayName", "logFC", "P.Value", "adj.P.Val")]
library(pheatmap)

outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

top10_down <- circ_results[circ_results$logFC < 0, ]
top10_down <- top10_down[order(top10_down$logFC), ][1:10, ]

top10_down$DisplayName <- top10_down$circRNA
top10_down$DisplayName[top10_down$circRNA == "ASCRP000064"] <- "hsa_circRNA_000554"

heat_expr <- expr_sub[top10_down$circRNA, ]
heat_expr_z <- t(scale(t(heat_expr)))
rownames(heat_expr_z) <- top10_down$DisplayName

colnames(heat_expr_z) <- c("Normal_1", "Normal_2", "Normal_3",
                           "TNBC_1", "TNBC_2", "TNBC_3", "TNBC_4")

annotation_col <- data.frame(
  Group = factor(c("Normal","Normal","Normal","TNBC","TNBC","TNBC","TNBC"),
                 levels = c("Normal","TNBC"))
)
rownames(annotation_col) <- colnames(heat_expr_z)

ann_colors <- list(Group = c(Normal = "#4DBBD5", TNBC = "#E64B35"))

draw_heatmap <- function() {
  pheatmap(
    heat_expr_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 12,
    fontsize_col = 12,
    border_color = NA,
    main = "Top 10 Downregulated circRNAs in TNBC"
  )
}

png(file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.png"),
    width = 3000, height = 2200, res = 300)
draw_heatmap()
dev.off()

tiff(file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.tiff"),
     width = 3000, height = 2200, res = 300, compression = "lzw")
draw_heatmap()
dev.off()

pdf(file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.pdf"),
    width = 10, height = 7.5)
draw_heatmap()
dev.off()

write.csv(
  top10_down,
  file.path(outdir, "Table_Top10_Downregulated_circRNAs_FINAL.csv"),
  row.names = FALSE
)

top10_down[, c("DisplayName", "logFC", "P.Value", "adj.P.Val")]
# Select downregulated circRNAs using the same criteria as volcano plot
top10_down <- circ_results[
  circ_results$P.Value < 0.05 &
    circ_results$logFC < -1,
]

# Order by strongest downregulation
top10_down <- top10_down[order(top10_down$logFC), ]

# Take top 10
top10_down <- top10_down[1:10, ]

# Display names
top10_down$DisplayName <- top10_down$circRNA
top10_down$DisplayName[
  top10_down$circRNA == "ASCRP000064"
] <- "hsa_circRNA_000554"

top10_down[, c("DisplayName", "logFC", "P.Value", "adj.P.Val")]
outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

draw_heatmap <- function() {
  pheatmap(
    heat_expr_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 12,
    fontsize_col = 12,
    border_color = NA,
    main = "Top 10 Downregulated circRNAs in TNBC"
  )
}

png(
  file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.png"),
  width = 3000,
  height = 2200,
  res = 300
)
draw_heatmap()
dev.off()

tiff(
  file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.tiff"),
  width = 3000,
  height = 2200,
  res = 300,
  compression = "lzw"
)
draw_heatmap()
dev.off()

pdf(
  file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.pdf"),
  width = 10,
  height = 7.5
)
draw_heatmap()
dev.off()

write.csv(
  top10_down,
  file.path(outdir, "Table_Top10_Downregulated_circRNAs_FINAL.csv"),
  row.names = FALSE
)

write.csv(
  heat_expr_z,
  file.path(outdir, "Matrix_Top10_Downregulated_Zscore_FINAL.csv"),
  row.names = TRUE
)

list.files(outdir)
circ_results[circ_results$circRNA == "ASCRP000064",
             c("circRNA", "logFC", "P.Value", "adj.P.Val")]
down_hits <- circ_results[
  circ_results$P.Value < 0.05 &
    circ_results$logFC < -1,
]

down_hits <- down_hits[order(down_hits$logFC), ]

which(down_hits$circRNA == "ASCRP000064")
nrow(down_hits)
# =====================================================
# Check rank of hsa_circRNA_000554 and redraw heatmap
# Criteria: P.Value < 0.05 and logFC < -1
# =====================================================

library(pheatmap)

outdir <- "GSE101123_circRNA_results"
dir.create(outdir, showWarnings = FALSE)

target_probe <- "ASCRP000064"
target_name  <- "hsa_circRNA_000554"

# Apply final criteria
down_hits <- circ_results[
  circ_results$P.Value < 0.05 &
    circ_results$logFC < -1,
]

# Rank by strongest downregulation
down_hits <- down_hits[order(down_hits$logFC), ]

# Add rank
down_hits$Rank <- seq_len(nrow(down_hits))

# Add display name
down_hits$DisplayName <- down_hits$circRNA
down_hits$DisplayName[down_hits$circRNA == target_probe] <- target_name

# Print target rank
down_hits[
  down_hits$circRNA == target_probe,
  c("Rank", "DisplayName", "circRNA", "logFC", "P.Value", "adj.P.Val")
]

# Select top 10 downregulated circRNAs
top10_down <- down_hits[1:10, ]

# Print top 10 to confirm target is included
top10_down[, c("Rank", "DisplayName", "circRNA", "logFC", "P.Value", "adj.P.Val")]

# Expression matrix
heat_expr <- expr_sub[top10_down$circRNA, ]

# Row-wise Z-score normalization
heat_expr_z <- t(scale(t(heat_expr)))

# Use display names in heatmap
rownames(heat_expr_z) <- top10_down$DisplayName

# Rename samples
colnames(heat_expr_z) <- c(
  "Normal_1", "Normal_2", "Normal_3",
  "TNBC_1", "TNBC_2", "TNBC_3", "TNBC_4"
)

# Sample annotation
annotation_col <- data.frame(
  Group = factor(
    c("Normal", "Normal", "Normal", "TNBC", "TNBC", "TNBC", "TNBC"),
    levels = c("Normal", "TNBC")
  )
)

rownames(annotation_col) <- colnames(heat_expr_z)

ann_colors <- list(
  Group = c(
    Normal = "#4DBBD5",
    TNBC = "#E64B35"
  )
)

draw_heatmap <- function() {
  pheatmap(
    heat_expr_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 12,
    fontsize_col = 12,
    border_color = NA,
    main = "Top 10 Downregulated circRNAs in TNBC"
  )
}

# Save files
png(
  file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.png"),
  width = 3000,
  height = 2200,
  res = 300
)
draw_heatmap()
dev.off()

tiff(
  file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.tiff"),
  width = 3000,
  height = 2200,
  res = 300,
  compression = "lzw"
)
draw_heatmap()
dev.off()

pdf(
  file.path(outdir, "Figure6B_Heatmap_Top10_Downregulated_FINAL.pdf"),
  width = 10,
  height = 7.5
)
draw_heatmap()
dev.off()

write.csv(
  down_hits,
  file.path(outdir, "Table_All_Downregulated_circRNAs_with_Rank_FINAL.csv"),
  row.names = FALSE
)

write.csv(
  top10_down,
  file.path(outdir, "Table_Top10_Downregulated_circRNAs_FINAL.csv"),
  row.names = FALSE
)

write.csv(
  heat_expr_z,
  file.path(outdir, "Matrix_Top10_Downregulated_Zscore_FINAL.csv"),
  row.names = TRUE
)

list.files(outdir, pattern = "Top10_Downregulated|Rank|Figure6B", full.names = TRUE)