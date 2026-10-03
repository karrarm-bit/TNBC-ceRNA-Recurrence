library(GEOquery)

candidate_gse <- c(
  "GSE111504",
  "GSE182471",
  "GSE101123",
  "GSE101124"
)

target_circ <- "hsa_circRNA_000554"

results <- data.frame()

for(gse_id in candidate_gse){
  
  cat("\nChecking:", gse_id, "\n")
  
  gse <- tryCatch(
    getGEO(gse_id, GSEMatrix = TRUE),
    error = function(e) NULL
  )
  
  if(is.null(gse)) next
  
  for(i in seq_along(gse)){
    
    eset <- gse[[i]]
    
    expr <- exprs(eset)
    fdat <- fData(eset)
    
    row_hit <- grep(
      target_circ,
      rownames(expr),
      ignore.case = TRUE,
      value = TRUE
    )
    
    fd_hit <- grep(
      target_circ,
      apply(fdat,1,paste,collapse=" "),
      ignore.case = TRUE,
      value = TRUE
    )
    
    tmp <- data.frame(
      GSE = gse_id,
      Platform = annotation(eset),
      Features = nrow(expr),
      Samples = ncol(expr),
      RowMatch = length(row_hit)>0,
      FdataMatch = length(fd_hit)>0
    )
    
    results <- rbind(results,tmp)
    
    if(length(row_hit)>0) print(row_hit)
    
    if(length(fd_hit)>0) print(fd_hit)
  }
}

print(results)
gse <- getGEO("GSE111504", GSEMatrix = TRUE)
eset <- gse[[1]]

fdat <- fData(eset)

hit <- fdat[
  apply(fdat, 1, function(x)
    any(grepl("hsa_circRNA_000554",
              x,
              ignore.case = TRUE))),
]

hit
expr <- exprs(eset)
pheno <- pData(eset)

pheno$title
library(GEOquery)

gse <- getGEO("GSE182471", GSEMatrix = TRUE)

length(gse)

eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)
fdat <- fData(eset)

dim(expr)
annotation(eset)
pheno$title
probe_id <- rownames(hit)[1]
probe_id

expr <- exprs(eset)

circ_expr <- expr[probe_id, ]

pheno$Group <- ifelse(
  grepl("adjacent", pheno$title, ignore.case = TRUE),
  "Normal",
  "Tumor"
)

df <- data.frame(
  Sample = colnames(expr),
  Expression = as.numeric(circ_expr),
  Group = pheno$Group
)

df$Group <- factor(df$Group,
                   levels = c("Normal", "Tumor"))

print(df)

table(df$Group)

wilcox.test(Expression ~ Group, data = df)
library(ggpubr)

p <- ggboxplot(
  df,
  x = "Group",
  y = "Expression",
  fill = "Group",
  add = "jitter"
) +
  stat_compare_means(method = "wilcox.test") +
  ggtitle("hsa_circRNA_000554 validation in GSE182471")

print(p)
outdir <- "09_GSE182471_External_circRNA_Validation"
dir.create(outdir, showWarnings = FALSE)

result_table <- data.frame(
  circRNA_original = "ASCRP000064",
  circRNA_external = "hsa_circRNA_000554",
  Probe_ID = probe_id,
  Dataset = "GSE182471",
  Platform = annotation(eset),
  Normal_N = sum(df$Group == "Normal"),
  Tumor_N = sum(df$Group == "Tumor"),
  Normal_mean = mean(df$Expression[df$Group == "Normal"]),
  Tumor_mean = mean(df$Expression[df$Group == "Tumor"]),
  Direction = "Downregulated in Tumor",
  Wilcoxon_W = 24,
  P_value = 0.01587
)

write.xlsx(
  list(
    Expression_Data = df,
    Validation_Result = result_table
  ),
  file = file.path(outdir, "GSE182471_hsa_circRNA_000554_external_validation.xlsx"),
  rowNames = FALSE
)

ggsave(
  filename = file.path(outdir, "GSE182471_hsa_circRNA_000554_boxplot.png"),
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "GSE182471_hsa_circRNA_000554_boxplot.tiff"),
  plot = p,
  width = 6,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

print(result_table)
library(ggplot2)
library(ggpubr)

p <- ggplot(df, aes(x = Group,
                    y = Expression,
                    fill = Group)) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  geom_jitter(
    width = 0.08,
    size = 2.8,
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 5
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
    plot.title = element_text(
      hjust = 0.5,
      face = "italic",
      size = 16
    ),
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 13
    ),
    legend.position = "none"
  )

print(p)
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
    label = "p = 0.016",
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
# ==============================
# Save final GSE182471 validation results
# ==============================

library(openxlsx)
library(ggplot2)

outdir <- "09_GSE182471_External_circRNA_Validation"
dir.create(outdir, showWarnings = FALSE)

# Statistical test
stat_test <- wilcox.test(Expression ~ Group, data = df)

# Summary table
summary_table <- df %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = dplyr::n(),
    .groups = "drop"
  )

# Final result table
result_table <- data.frame(
  circRNA_original = "ASCRP000064",
  circRNA_external = "hsa_circRNA_000554",
  Probe_ID = "ASCRP3009378",
  Dataset = "GSE182471",
  Platform = annotation(eset),
  Normal_N = sum(df$Group == "Normal"),
  Tumor_N = sum(df$Group == "Tumor"),
  Normal_mean = mean(df$Expression[df$Group == "Normal"]),
  Tumor_mean = mean(df$Expression[df$Group == "Tumor"]),
  Direction = "Downregulated in Tumor",
  Wilcoxon_W = as.numeric(stat_test$statistic),
  P_value = stat_test$p.value
)

# Save Excel
write.xlsx(
  list(
    Expression_Data = df,
    Summary = summary_table,
    Validation_Result = result_table
  ),
  file = file.path(outdir, "GSE182471_hsa_circRNA_000554_external_validation_FINAL.xlsx"),
  rowNames = FALSE
)

# Save CSV
write.csv(
  result_table,
  file = file.path(outdir, "GSE182471_validation_result_FINAL.csv"),
  row.names = FALSE
)

write.csv(
  df,
  file = file.path(outdir, "GSE182471_expression_data_FINAL.csv"),
  row.names = FALSE
)

# Save figure
ggsave(
  filename = file.path(outdir, "Figure_GSE182471_hsa_circRNA_000554_FINAL.png"),
  plot = p,
  width = 4.5,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "Figure_GSE182471_hsa_circRNA_000554_FINAL.tiff"),
  plot = p,
  width = 4.5,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(outdir, "Figure_GSE182471_hsa_circRNA_000554_FINAL.pdf"),
  plot = p,
  width = 4.5,
  height = 5
)

print(result_table)
list.files(outdir)