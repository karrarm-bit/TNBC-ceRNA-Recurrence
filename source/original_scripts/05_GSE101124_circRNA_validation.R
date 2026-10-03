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

# Download GSE101124

gse_circ_val <- getGEO(
  "GSE101124",
  GSEMatrix = TRUE
)

length(gse_circ_val)

gse_circ_val

eset_circ_val <- gse_circ_val[[1]]

expr_circ_val <- exprs(eset_circ_val)

dim(expr_circ_val)

pheno_circ_val <- pData(eset_circ_val)

colnames(pheno_circ_val)

View(pheno_circ_val)
table(pheno_circ_val$title)
dput(unique(pheno_circ_val$`subtype:ch1`))
subtype <- trimws(pheno_circ_val$`subtype:ch1`)

keep <- subtype %in% c(
  "Triple-negative breast cancer",
  "non-tumor breast tissues"
)

expr_sub <- expr_circ_val[, keep]

pheno_sub <- pheno_circ_val[keep, ]

group_circ <- factor(
  ifelse(
    trimws(pheno_sub$`subtype:ch1`) ==
      "Triple-negative breast cancer",
    "TNBC",
    "Normal"
  ),
  levels = c("Normal", "TNBC")
)

table(group_circ)

dim(expr_sub)
dput(unique(pheno_circ_val$`subtype:ch1`))
subtype <- trimws(pheno_circ_val$`subtype:ch1`)

keep <- subtype %in% c(
  "subtype: Triple-negative breast cancer",
  "subtype: non-tumor breast tissues"
)

expr_sub <- expr_circ_val[, keep]

pheno_sub <- pheno_circ_val[keep, ]

group_circ <- factor(
  ifelse(
    trimws(pheno_sub$`subtype:ch1`) ==
      "subtype: Triple-negative breast cancer",
    "TNBC",
    "Normal"
  ),
  levels = c("Normal", "TNBC")
)

table(group_circ)
dim(expr_sub)
# Check exact values
dput(unique(pheno_circ_val$characteristics_ch1.1))

# Robust filtering using grepl
subtype <- trimws(pheno_circ_val$characteristics_ch1.1)

keep <- grepl("Triple-negative", subtype, ignore.case = TRUE) |
  grepl("non-tumor", subtype, ignore.case = TRUE)

expr_sub <- expr_circ_val[, keep, drop = FALSE]
pheno_sub <- pheno_circ_val[keep, ]

group_circ <- factor(
  ifelse(
    grepl("Triple-negative", trimws(pheno_sub$characteristics_ch1.1), ignore.case = TRUE),
    "TNBC",
    "Normal"
  ),
  levels = c("Normal", "TNBC")
)

table(group_circ)
dim(expr_sub)

pheno_sub[, c("title", "geo_accession", "characteristics_ch1.1")]
group_circ <- factor(
  ifelse(
    grepl("Triple-negative", trimws(pheno_sub$characteristics_ch1.1), ignore.case = TRUE),
    "TNBC",
    "Normal"
  ),
  levels = c("Normal", "TNBC")
)

table(group_circ)
dim(expr_sub)
grep("ASCRP000064", rownames(expr_sub), value = TRUE)
# Validate ASCRP000064 in GSE101124

circ_id <- "ASCRP000064"

x <- as.numeric(expr_sub[circ_id, ])

validation_circ <- data.frame(
  circRNA = circ_id,
  TNBC_mean = mean(x[group_circ == "TNBC"], na.rm = TRUE),
  Normal_mean = mean(x[group_circ == "Normal"], na.rm = TRUE),
  Direction = ifelse(
    mean(x[group_circ == "TNBC"], na.rm = TRUE) >
      mean(x[group_circ == "Normal"], na.rm = TRUE),
    "Up",
    "Down"
  ),
  P_value = signif(
    wilcox.test(
      x[group_circ == "TNBC"],
      x[group_circ == "Normal"]
    )$p.value,
    3
  )
)

validation_circ
dir.create("GSE101124_circRNA_validation_results", showWarnings = FALSE)

write.csv(
  validation_circ,
  "GSE101124_circRNA_validation_results/Table_GSE101124_ASCRP000064_validation.csv",
  row.names = FALSE
)
library(ggplot2)

plot_df <- data.frame(
  Expression = x,
  Group = group_circ
)

p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("Normal" = "#56B4E9",
                               "TNBC" = "#E64B35")) +
  labs(
    title = "Validation of ASCRP000064 in GSE101124",
    x = "",
    y = "Expression Level"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(plot_df$Expression) * 1.05,
    label = paste0("P = ", signif(validation_circ$P_value, 3))
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p
ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.tiff",
  plot = p,
  width = 5,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.png",
  plot = p,
  width = 5,
  height = 5,
  dpi = 600
)
library(ggplot2)

# Prepare plotting data
plot_df <- data.frame(
  Expression = x,
  Group = group_circ
)

# Generate boxplot
p <- ggplot(plot_df, aes(x = Group,
                         y = Expression,
                         fill = Group)) +
  
  geom_boxplot(
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 3,
    alpha = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "Normal" = "#56B4E9",
      "TNBC"   = "#E64B35"
    )
  ) +
  
  labs(
    title = "Validation of ASCRP000064 in GSE101124",
    x = "",
    y = "Expression Level"
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = max(plot_df$Expression) + 0.25,
    label = paste0(
      "Wilcoxon P = ",
      signif(validation_circ$P_value, 3)
    ),
    size = 6,
    fontface = "bold"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Display plot
p
library(ggplot2)

plot_df <- data.frame(
  Expression = x,
  Group = group_circ
)

y_max <- max(plot_df$Expression, na.rm = TRUE)
y_min <- min(plot_df$Expression, na.rm = TRUE)

p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(
    alpha = 0.85,
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#56B4E9",
      "TNBC" = "#E64B35"
    )
  ) +
  scale_y_continuous(
    limits = c(y_min - 0.2, y_max + 0.7),
    expand = c(0, 0)
  ) +
  labs(
    title = "Validation of ASCRP000064 in GSE101124",
    x = "",
    y = "Expression Level"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = y_max + 0.45,
    label = "Wilcoxon P = 0.057",
    size = 5.5,
    fontface = "bold"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15,
      margin = margin(b = 12)
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

print(p)

ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.tiff",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.png",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 600
)
library(ggplot2)

plot_df <- data.frame(
  Expression = x,
  Group = group_circ
)

y_max <- max(plot_df$Expression, na.rm = TRUE)
y_min <- min(plot_df$Expression, na.rm = TRUE)

p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(
    alpha = 0.85,
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#56B4E9",
      "TNBC" = "#E64B35"
    )
  ) +
  scale_y_continuous(
    limits = c(y_min - 0.2, y_max + 0.8),
    expand = c(0, 0)
  ) +
  labs(
    title = "Validation of ASCRP000064 in GSE101124",
    x = "",
    y = "Expression Level"
  ) +
  geom_segment(aes(x = 1, xend = 2, y = y_max + 0.35, yend = y_max + 0.35),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(aes(x = 1, xend = 1, y = y_max + 0.30, yend = y_max + 0.35),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = y_max + 0.30, yend = y_max + 0.35),
               inherit.aes = FALSE, linewidth = 0.8) +
  annotate(
    "text",
    x = 1.5,
    y = y_max + 0.52,
    label = "P = 0.057",
    size = 6,
    fontface = "bold"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15,
      margin = margin(b = 12)
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

print(p)

ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.tiff",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.png",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 600
)
library(ggplot2)

# Prepare plotting data
plot_df <- data.frame(
  Expression = x,
  Group = group_circ
)

# Define y-axis limits
y_max <- max(plot_df$Expression, na.rm = TRUE)
y_min <- min(plot_df$Expression, na.rm = TRUE)

# Create publication-quality boxplot
p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
  
  geom_boxplot(
    alpha = 0.85,
    width = 0.55,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.85
  ) +
  
  scale_fill_manual(
    values = c(
      "Normal" = "#56B4E9",
      "TNBC" = "#E64B35"
    )
  ) +
  
  scale_y_continuous(
    limits = c(y_min - 0.2, y_max + 0.6),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Validation of ASCRP000064 in GSE101124",
    x = "",
    y = "Expression Level"
  ) +
  
  # Significance bracket
  geom_segment(
    aes(x = 1, xend = 2, y = y_max + 0.20, yend = y_max + 0.20),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  geom_segment(
    aes(x = 1, xend = 1, y = y_max + 0.15, yend = y_max + 0.20),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  geom_segment(
    aes(x = 2, xend = 2, y = y_max + 0.15, yend = y_max + 0.20),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  # P-value label
  annotate(
    "text",
    x = 1.5,
    y = y_max + 0.28,
    label = "P = 0.057",
    size = 5.5,
    fontface = "bold"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15,
      margin = margin(b = 12)
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

# Display plot
print(p)

# Save TIFF
ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.tiff",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# Save PNG
ggsave(
  "GSE101124_circRNA_validation_results/Figure_ASCRP000064_validation_boxplot.png",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 600
)
# Create folder for target prediction
dir.create("06_miRNA_target_prediction_results", showWarnings = FALSE)

# Read TargetScan downloaded table
targetscan_590 <- read.delim(
  "TargetScan_miR590_5p.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# View data
dim(targetscan_590)
head(targetscan_590)
colnames(targetscan_590)

# Extract target genes
targetscan_590_genes <- unique(targetscan_590$`Target gene`)

length(targetscan_590_genes)
head(targetscan_590_genes)

# Save clean gene list
write.csv(
  data.frame(Gene = targetscan_590_genes),
  "06_miRNA_target_prediction_results/TargetScan_miR590_5p_genes.csv",
  row.names = FALSE
)
getwd()

list.files()
targetscan_590 <- read.delim(
  "06_miRNA_target_prediction_results/TargetScan_miR590_5p.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(targetscan_590)
head(targetscan_590)
colnames(targetscan_590)
list.files("06_miRNA_target_prediction_results")
targetscan_590 <- read.delim(
  "06_miRNA_target_prediction_results/TargetScan_miR590_5p.txt.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(targetscan_590)

targetscan_590_genes <- unique(targetscan_590$`Target gene`)

length(targetscan_590_genes)
head(targetscan_590_genes)

write.csv(
  data.frame(Gene = targetscan_590_genes),
  "06_miRNA_target_prediction_results/TargetScan_miR590_5p_genes.csv",
  row.names = FALSE
)
files <- list.files("06_miRNA_target_prediction_results", full.names = TRUE)
file_182 <- files[grepl("182", files)]

targetscan_182 <- read.delim(
  file_182,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dim(targetscan_182)
head(targetscan_182)
colnames(targetscan_182)
targetscan_182_genes <- unique(targetscan_182$`Target gene`)

length(targetscan_182_genes)
head(targetscan_182_genes)

write.csv(
  data.frame(Gene = targetscan_182_genes),
  "06_miRNA_target_prediction_results/TargetScan_miR182_5p_genes.csv",
  row.names = FALSE
)
mirdb_590_genes <- unique(mirdb_590$`Gene Symbol`)
mirdb_182_genes <- unique(mirdb_182$`Gene Symbol`)

length(mirdb_590_genes)
length(mirdb_182_genes)

write.csv(
  data.frame(Gene = mirdb_590_genes),
  "06_miRNA_target_prediction_results/miRDB_miR590_5p_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = mirdb_182_genes),
  "06_miRNA_target_prediction_results/miRDB_miR182_5p_genes.csv",
  row.names = FALSE
)
list.files("06_miRNA_target_prediction_results")
colnames(mirtarbase_590)
colnames(mirdb_590)
mirtarbase_590 <- read.csv(
  "06_miRNA_target_prediction_results/miRTarBase_miR590_5p_validated_targets.csv.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mirtarbase_182 <- read.csv(
  "06_miRNA_target_prediction_results/miRTarBase_miR182_5p_validated_targets.csv.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(mirtarbase_590)
dim(mirtarbase_182)

colnames(mirtarbase_590)
colnames(mirtarbase_182)
mirtarbase_590_genes <- unique(mirtarbase_590$Target)

mirtarbase_182_genes <- unique(mirtarbase_182$Target)

length(mirtarbase_590_genes)
length(mirtarbase_182_genes)

head(mirtarbase_590_genes)
head(mirtarbase_182_genes)
write.csv(
  data.frame(Gene = mirtarbase_590_genes),
  "06_miRNA_target_prediction_results/miRTarBase_miR590_5p_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = mirtarbase_182_genes),
  "06_miRNA_target_prediction_results/miRTarBase_miR182_5p_genes.csv",
  row.names = FALSE
)
common_590 <- Reduce(
  intersect,
  list(
    targetscan_590_genes,
    mirdb_590_genes,
    mirtarbase_590_genes
  )
)

length(common_590)
common_590
grep("mirdb", ls(), value = TRUE)
# Read miRDB files automatically
files <- list.files(
  "06_miRNA_target_prediction_results",
  full.names = TRUE
)

file_mirdb_590 <- files[grepl("miRDB.*590", files, ignore.case = TRUE)]
file_mirdb_182 <- files[grepl("miRDB.*182", files, ignore.case = TRUE)]

file_mirdb_590
file_mirdb_182

library(readxl)

mirdb_590 <- read_excel(file_mirdb_590)
mirdb_182 <- read_excel(file_mirdb_182)

dim(mirdb_590)
dim(mirdb_182)

colnames(mirdb_590)
colnames(mirdb_182)

mirdb_590_genes <- unique(mirdb_590$`Gene Symbol`)
mirdb_182_genes <- unique(mirdb_182$`Gene Symbol`)

length(mirdb_590_genes)
length(mirdb_182_genes)

write.csv(
  data.frame(Gene = mirdb_590_genes),
  "06_miRNA_target_prediction_results/miRDB_miR590_5p_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = mirdb_182_genes),
  "06_miRNA_target_prediction_results/miRDB_miR182_5p_genes.csv",
  row.names = FALSE
)
common_590 <- Reduce(
  intersect,
  list(
    targetscan_590_genes,
    mirdb_590_genes,
    mirtarbase_590_genes
  )
)

length(common_590)
common_590
common_182 <- Reduce(
  intersect,
  list(
    targetscan_182_genes,
    mirdb_182_genes,
    mirtarbase_182_genes
  )
)

length(common_182)
common_182
write.csv(
  data.frame(miRNA = "hsa-miR-590-5p", Gene = common_590),
  "06_miRNA_target_prediction_results/Common_targets_miR590_5p_3databases.csv",
  row.names = FALSE
)

write.csv(
  data.frame(miRNA = "hsa-miR-182-5p", Gene = common_182),
  "06_miRNA_target_prediction_results/Common_targets_miR182_5p_3databases.csv",
  row.names = FALSE
)

all_common_targets <- rbind(
  data.frame(miRNA = "hsa-miR-590-5p", Gene = common_590),
  data.frame(miRNA = "hsa-miR-182-5p", Gene = common_182)
)

write.csv(
  all_common_targets,
  "06_miRNA_target_prediction_results/Common_targets_all_validated_miRNAs_3databases.csv",
  row.names = FALSE
)
sink("06_miRNA_target_prediction_results/README_target_prediction_summary.txt")

cat("miRNA target prediction summary\n")
cat("Databases used: TargetScan, miRDB, miRTarBase\n\n")

cat("hsa-miR-590-5p common targets:", length(common_590), "\n")
cat(common_590, sep = ", ")
cat("\n\n")

cat("hsa-miR-182-5p common targets:", length(common_182), "\n")
cat(common_182, sep = ", ")
cat("\n\n")

cat("Total high-confidence targets:", nrow(all_common_targets), "\n")

sink()