rm(list = ls())
gc()
rm(list = ls())
gc()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("GEOquery", "limma", "ggplot2", "pheatmap", "dplyr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GEOquery", "limma")) {
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

gse_val2 <- getGEO("GSE38167", GSEMatrix = TRUE)

length(gse_val2)

gse_val2

eset_val2 <- gse_val2[[1]]

expr_val2 <- exprs(eset_val2)

dim(expr_val2)

pheno_val2 <- pData(eset_val2)

colnames(pheno_val2)

View(pheno_val2)
table(pheno_val2$`tissue:ch1`)
## Keep only Tumour and Normal samples
keep <- pheno_val2$`tissue:ch1` %in%
  c("Breast tumour", "Normal adjacent breast")

expr_val <- expr_val2[, keep]
pheno_sub <- pheno_val2[keep, ]

## Define groups
group_val <- factor(
  ifelse(pheno_sub$`tissue:ch1`=="Breast tumour",
         "Tumour",
         "Normal")
)

table(group_val)
fdata <- fData(eset_val2)

colnames(fdata)
grep("miR", colnames(fdata), value = TRUE)
target_miRNAs <- c(
  "hsa-miR-590-5p",
  "hsa-miR-182-5p",
  "hsa-miR-3125",
  "hsa-miR-452-5p"
)

fdata[fdata$miRNA_ID %in% target_miRNAs,
      c("ProbeName","miRNA_ID")]
library(dplyr)

mir_info <- fdata %>%
  filter(miRNA_ID %in% target_miRNAs)

mir_info
selected_probes <- data.frame()

for(m in target_miRNAs){
  
  tmp <- mir_info[mir_info$miRNA_ID == m, ]
  
  means <- apply(
    expr_val[tmp$ID, , drop = FALSE],
    1,
    mean,
    na.rm = TRUE
  )
  
  best_probe <- names(which.max(means))
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      miRNA = m,
      Probe = best_probe
    )
  )
}

selected_probes
library(dplyr)

mir_info <- fdata %>%
  filter(miRNA_ID %in% target_miRNAs)

selected_probes <- data.frame()

for(m in target_miRNAs){
  
  tmp <- mir_info[mir_info$miRNA_ID == m, ]
  
  tmp <- tmp[tmp$ProbeName %in% rownames(expr_val), ]
  
  means <- apply(
    expr_val[tmp$ProbeName, , drop = FALSE],
    1,
    mean,
    na.rm = TRUE
  )
  
  best_probe <- names(which.max(means))
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      miRNA = m,
      Probe = best_probe
    )
  )
}

selected_probes
# Check rownames and feature columns
head(rownames(expr_val))
head(fdata$ID)
head(fdata$ProbeName)
head(fdata$SystematicName)
head(fdata$miRNA_ID)

# Which column matches rownames(expr_val)?
sum(fdata$ID %in% rownames(expr_val))
sum(fdata$ProbeName %in% rownames(expr_val))
sum(fdata$SystematicName %in% rownames(expr_val))
library(dplyr)

mir_info <- fdata %>%
  filter(miRNA_ID %in% target_miRNAs)

selected_probes <- data.frame()

for (m in target_miRNAs) {
  
  tmp <- mir_info[mir_info$miRNA_ID == m, ]
  
  tmp <- tmp[tmp$ID %in% rownames(expr_val), ]
  
  means <- apply(
    expr_val[tmp$ID, , drop = FALSE],
    1,
    mean,
    na.rm = TRUE
  )
  
  best_probe <- names(which.max(means))
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      miRNA = m,
      Probe = best_probe
    )
  )
}

selected_probes
library(dplyr)

mir_info <- fdata %>%
  filter(miRNA_ID %in% target_miRNAs)

selected_probes <- data.frame()

for (m in target_miRNAs) {
  
  tmp <- mir_info[mir_info$miRNA_ID == m, ]
  
  tmp <- tmp[tmp$ID %in% rownames(expr_val2), ]
  
  means <- apply(
    expr_val2[tmp$ID, , drop = FALSE],
    1,
    mean,
    na.rm = TRUE
  )
  
  best_probe <- names(which.max(means))
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      miRNA = m,
      Probe = best_probe
    )
  )
}

selected_probes
expr_final <- expr_val2[selected_probes$Probe, keep]

rownames(expr_final) <- selected_probes$miRNA

validation_results <- data.frame()

for (i in 1:nrow(expr_final)) {
  
  x <- as.numeric(expr_final[i, ])
  
  test <- wilcox.test(
    x[group_val == "Tumour"],
    x[group_val == "Normal"]
  )
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      miRNA = rownames(expr_final)[i],
      Tumour_mean = mean(x[group_val == "Tumour"]),
      Normal_mean = mean(x[group_val == "Normal"]),
      Direction = ifelse(
        mean(x[group_val == "Tumour"]) > mean(x[group_val == "Normal"]),
        "Up",
        "Down"
      ),
      P_value = signif(test$p.value, 3)
    )
  )
}

validation_results
library(dplyr)

mir_info <- fdata %>%
  filter(miRNA_ID %in% target_miRNAs)

selected_probes <- data.frame()

for (m in target_miRNAs) {
  
  tmp <- mir_info[mir_info$miRNA_ID == m, ]
  tmp <- tmp[tmp$ID %in% rownames(expr_val2), ]
  
  if (nrow(tmp) == 0) {
    message(m, " not found in GSE38167")
    next
  }
  
  means <- apply(
    expr_val2[tmp$ID, , drop = FALSE],
    1,
    mean,
    na.rm = TRUE
  )
  
  best_probe <- names(which.max(means))
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      miRNA = m,
      Probe = best_probe
    )
  )
}

selected_probes
expr_final <- expr_val2[selected_probes$Probe, keep, drop = FALSE]

rownames(expr_final) <- selected_probes$miRNA

dim(expr_final)
rownames(expr_final)
validation_results <- data.frame()

for (i in seq_len(nrow(expr_final))) {
  
  x <- as.numeric(expr_final[i, ])
  
  test <- wilcox.test(
    x[group_val == "Tumour"],
    x[group_val == "Normal"]
  )
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      miRNA = rownames(expr_final)[i],
      Tumour_mean = mean(x[group_val == "Tumour"]),
      Normal_mean = mean(x[group_val == "Normal"]),
      Direction = ifelse(
        mean(x[group_val == "Tumour"]) > mean(x[group_val == "Normal"]),
        "Up",
        "Down"
      ),
      P_value = signif(test$p.value, 3)
    )
  )
}

validation_results
library(dplyr)

## Target miRNAs from discovery
target_miRNAs <- c(
  "hsa-miR-590-5p",
  "hsa-miR-182-5p",
  "hsa-miR-3125",
  "hsa-miR-452-5p"
)

## Clean miRNA IDs
fdata$miRNA_ID_clean <- trimws(fdata$miRNA_ID)

## Select best row/probe for each miRNA using row index
selected_probes <- data.frame()

for (m in target_miRNAs) {
  
  idx <- which(fdata$miRNA_ID_clean == m)
  
  if (length(idx) == 0) {
    message(m, " not found")
    next
  }
  
  expr_tmp <- expr_val2[idx, keep, drop = FALSE]
  
  means <- rowMeans(expr_tmp, na.rm = TRUE)
  
  best_idx <- idx[which.max(means)]
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      miRNA = m,
      RowIndex = best_idx,
      ProbeID = rownames(expr_val2)[best_idx]
    )
  )
}

selected_probes
expr_final <- expr_val2[
  selected_probes$RowIndex,
  keep,
  drop = FALSE
]

rownames(expr_final) <- selected_probes$miRNA

dim(expr_final)
rownames(expr_final)
validation_results <- data.frame()

for (i in seq_len(nrow(expr_final))) {
  
  x <- as.numeric(expr_final[i, ])
  
  test <- wilcox.test(
    x[group_val == "Tumour"],
    x[group_val == "Normal"]
  )
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      miRNA = rownames(expr_final)[i],
      Tumour_mean = mean(x[group_val == "Tumour"]),
      Normal_mean = mean(x[group_val == "Normal"]),
      Direction = ifelse(
        mean(x[group_val == "Tumour"]) > mean(x[group_val == "Normal"]),
        "Up",
        "Down"
      ),
      P_value = signif(test$p.value, 3)
    )
  )
}

validation_results
target_miRNAs <- c(
  "hsa-miR-590-5p",
  "hsa-miR-182-5p",
  "hsa-miR-3125",
  "hsa-miR-452-5p"
)

sapply(target_miRNAs, function(m) {
  sum(trimws(fdata$miRNA_ID) == m)
})
grep("182", fdata$miRNA_ID, value = TRUE)
grep("3125", fdata$miRNA_ID, value = TRUE)
grep("452", fdata$miRNA_ID, value = TRUE)
target_map <- data.frame(
  discovery_miRNA = c("hsa-miR-590-5p", "hsa-miR-182-5p", "hsa-miR-452-5p"),
  validation_miRNA = c("hsa-miR-590-5p", "hsa-miR-182", "hsa-miR-452")
)

selected_probes <- data.frame()

for (j in 1:nrow(target_map)) {
  
  m_val <- target_map$validation_miRNA[j]
  m_dis <- target_map$discovery_miRNA[j]
  
  idx <- which(trimws(fdata$miRNA_ID) == m_val)
  
  expr_tmp <- expr_val2[idx, keep, drop = FALSE]
  means <- rowMeans(expr_tmp, na.rm = TRUE)
  best_idx <- idx[which.max(means)]
  
  selected_probes <- rbind(
    selected_probes,
    data.frame(
      Discovery_miRNA = m_dis,
      Validation_miRNA = m_val,
      RowIndex = best_idx,
      ProbeID = rownames(expr_val2)[best_idx]
    )
  )
}

selected_probes

expr_final <- expr_val2[selected_probes$RowIndex, keep, drop = FALSE]
rownames(expr_final) <- selected_probes$Discovery_miRNA

validation_results <- data.frame()

for (i in seq_len(nrow(expr_final))) {
  
  x <- as.numeric(expr_final[i, ])
  
  test <- wilcox.test(
    x[group_val == "Tumour"],
    x[group_val == "Normal"]
  )
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      miRNA = rownames(expr_final)[i],
      Validation_ID = selected_probes$Validation_miRNA[i],
      Tumour_mean = mean(x[group_val == "Tumour"]),
      Normal_mean = mean(x[group_val == "Normal"]),
      Direction = ifelse(
        mean(x[group_val == "Tumour"]) > mean(x[group_val == "Normal"),
        "Up",
        "Down"
      ),
      P_value = signif(test$p.value, 3)
    )
  )
}

validation_results
validation_results <- data.frame()

for (i in seq_len(nrow(expr_final))) {
  
  x <- as.numeric(expr_final[i, ])
  
  tumour_mean <- mean(x[group_val == "Tumour"])
  normal_mean <- mean(x[group_val == "Normal"])
  
  test <- wilcox.test(
    x[group_val == "Tumour"],
    x[group_val == "Normal"]
  )
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      miRNA = rownames(expr_final)[i],
      Validation_ID = selected_probes$Validation_miRNA[i],
      Tumour_mean = tumour_mean,
      Normal_mean = normal_mean,
      Direction = ifelse(tumour_mean > normal_mean, "Up", "Down"),
      P_value = signif(test$p.value, 3)
    )
  )
}

validation_results
library(ggplot2)
library(reshape2)

# Select validated miRNAs
validated_mirs <- c("hsa-miR-590-5p", "hsa-miR-182-5p")

# Extract expression values
expr_plot <- expr_final[validated_mirs, ]

# Convert expression matrix to long format
plot_df <- data.frame(t(expr_plot))
plot_df$Group <- group_val

plot_long <- melt(
  plot_df,
  id.vars = "Group",
  variable.name = "miRNA",
  value.name = "Expression"
)

# Create validation boxplots
p <- ggplot(
  plot_long,
  aes(x = Group, y = Expression, fill = Group)
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.8,
    outlier.shape = 16
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.5,
    alpha = 0.7
  ) +
  facet_wrap(
    ~miRNA,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#4DBBD5",
      "Tumour" = "#E64B35"
    )
  ) +
  labs(
    title = "External Validation of Candidate miRNAs in GSE38167",
    x = "",
    y = "Expression Level"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

# Display plot
print(p)

# Save figure
ggsave(
  "Figure5_Validation_Boxplots.tiff",
  plot = p,
  width = 8,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "Figure5_Validation_Boxplots.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 600
)
# Create folder for GSE38167 validation results
dir.create("GSE38167_validation_results", showWarnings = FALSE)

# Save validation results table
write.csv(
  validation_results,
  "GSE38167_validation_results/Table_GSE38167_validation_results.csv",
  row.names = FALSE
)

# Save selected probes table
write.csv(
  selected_probes,
  "GSE38167_validation_results/Table_GSE38167_selected_probes.csv",
  row.names = FALSE
)

# Save expression matrix of validated miRNAs
write.csv(
  expr_final,
  "GSE38167_validation_results/Table_GSE38167_validated_miRNA_expression.csv"
)

# Save validation boxplot again inside the folder
ggsave(
  "GSE38167_validation_results/Figure5_Validation_Boxplots.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 600
)

ggsave(
  "GSE38167_validation_results/Figure5_Validation_Boxplots.tiff",
  plot = p,
  width = 8,
  height = 5,
  dpi = 600,
  compression = "lzw"
)
# Create final validation summary table
validation_summary <- validation_results

validation_summary$Discovery_direction <- c("Up", "Up", "Down")

validation_summary$Validation_status <- ifelse(
  validation_summary$Discovery_direction == validation_summary$Direction &
    validation_summary$P_value < 0.05,
  "Validated",
  "Not validated"
)

validation_summary <- validation_summary[, c(
  "miRNA",
  "Validation_ID",
  "Discovery_direction",
  "Direction",
  "Tumour_mean",
  "Normal_mean",
  "P_value",
  "Validation_status"
)]

colnames(validation_summary) <- c(
  "miRNA",
  "Validation ID",
  "Discovery direction",
  "Validation direction",
  "Tumour mean",
  "Normal mean",
  "P-value",
  "Validation status"
)

validation_summary

write.csv(
  validation_summary,
  "GSE38167_validation_results/Table_GSE38167_final_validation_summary.csv",
  row.names = FALSE
)