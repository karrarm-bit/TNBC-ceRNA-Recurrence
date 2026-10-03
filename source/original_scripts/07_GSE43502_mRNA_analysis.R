############################################################
# GSE43502 mRNA analysis
# TNBC vs Normal
############################################################

# Load packages
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

# Create results folder
dir.create("07_GSE43502_mRNA_results",
           showWarnings = FALSE)

# Download GEO dataset
gse <- getGEO("GSE43502",
              GSEMatrix = TRUE)

length(gse)

# Extract ExpressionSet
eset <- gse[[1]]

# Expression matrix
expr <- exprs(eset)

# Phenotype data
pheno <- pData(eset)

# Check dimensions
dim(expr)

# View phenotype columns
colnames(pheno)

# First samples
head(pheno)
length(gse)
eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)

dim(expr)
colnames(pheno)
head(pheno)
table(pheno$`prognosis:ch1`)
table(pheno$`disease state:ch1`)
group_mrna <- factor(
  ifelse(
    pheno$`prognosis:ch1` == "recurrence",
    "Recurrence",
    "No_rec",
  ),
  levels = c("No_rec", "Recurrence")
)

table(group_mrna)
colnames(pheno)
grep("prognosis", colnames(pheno), value = TRUE)
group_mrna <- factor(
  ifelse(
    pheno[["prognosis:ch1"]] == "recurrence",
    "Recurrence",
    "No_rec"
  ),
  levels = c("No_rec", "Recurrence")
)

table(group_mrna)
# Differential expression: Recurrence vs No recurrence

design <- model.matrix(~ 0 + group_mrna)
colnames(design) <- levels(group_mrna)

contrast_matrix <- makeContrasts(
  Recurrence_vs_NoRec = Recurrence - No_rec,
  levels = design
)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

mrna_results <- topTable(
  fit2,
  coef = "Recurrence_vs_NoRec",
  number = Inf,
  adjust.method = "BH"
)

head(mrna_results)
dim(mrna_results)
annotation(eset)
library(hgu133plus2.db)
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  BiocManager::install("hgu133plus2.db")
}

library(hgu133plus2.db)
library(AnnotationDbi)
probe_ids <- rownames(mrna_results)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

mrna_results$ProbeID <- rownames(mrna_results)
mrna_results$Gene <- gene_symbols

head(mrna_results)
library(hgu133plus2.db)
library(AnnotationDbi)
probe_ids <- rownames(mrna_results)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

mrna_results$ProbeID <- rownames(mrna_results)
mrna_results$Gene <- gene_symbols

head(mrna_results)
# Remove probes without gene symbols
mrna_annot <- mrna_results[!is.na(mrna_results$Gene) & mrna_results$Gene != "", ]

# Keep best probe per gene using smallest adjusted p-value
mrna_annot <- mrna_annot[order(mrna_annot$adj.P.Val), ]

mrna_gene <- mrna_annot[!duplicated(mrna_annot$Gene), ]

dim(mrna_gene)
head(mrna_gene)

# Define recurrence-associated DEGs
sig_mrna <- mrna_gene[
  mrna_gene$adj.P.Val < 0.05 & abs(mrna_gene$logFC) > 1,
]

dim(sig_mrna)
head(sig_mrna)
final_network_genes <- sig_mrna[
  sig_mrna$Gene %in% all_common_targets$Gene,
]

dim(final_network_genes)
final_network_genes
dim(sig_mrna)
dim(final_network_genes)
final_network_genes
# Less strict recurrence-associated genes
sig_mrna_p <- mrna_gene[
  mrna_gene$P.Value < 0.05 & abs(mrna_gene$logFC) > 0.5,
]

dim(sig_mrna_p)

final_network_genes_p <- sig_mrna_p[
  sig_mrna_p$Gene %in% all_common_targets$Gene,
]

dim(final_network_genes_p)
final_network_genes_p
common_590
common_182
# Final recurrence-related network genes
final_ceRNA_genes <- final_network_genes_p

# Add miRNA source
final_ceRNA_genes$miRNA <- ifelse(
  final_ceRNA_genes$Gene %in% common_590,
  "hsa-miR-590-5p",
  "hsa-miR-182-5p"
)

# Arrange columns
final_ceRNA_genes <- final_ceRNA_genes[
  , c("Gene", "miRNA", "logFC", "P.Value", "adj.P.Val", "ProbeID")
]

# Save final table
write.csv(
  final_ceRNA_genes,
  "07_GSE43502_mRNA_results/Table_final_ceRNA_recurrence_genes.csv",
  row.names = FALSE
)

# Cytoscape edge file
cytoscape_edges <- rbind(
  data.frame(Source = "ASCRP000064", Target = "hsa-miR-590-5p", Interaction = "circRNA-miRNA"),
  data.frame(Source = "ASCRP000064", Target = "hsa-miR-182-5p", Interaction = "circRNA-miRNA"),
  data.frame(Source = "hsa-miR-590-5p", Target = final_ceRNA_genes$Gene[final_ceRNA_genes$miRNA == "hsa-miR-590-5p"], Interaction = "miRNA-mRNA"),
  data.frame(Source = "hsa-miR-182-5p", Target = final_ceRNA_genes$Gene[final_ceRNA_genes$miRNA == "hsa-miR-182-5p"], Interaction = "miRNA-mRNA")
)

write.csv(
  cytoscape_edges,
  "07_GSE43502_mRNA_results/Cytoscape_ceRNA_network_edges.csv",
  row.names = FALSE
)

# Cytoscape node file
cytoscape_nodes <- data.frame(
  Node = unique(c(cytoscape_edges$Source, cytoscape_edges$Target))
)

cytoscape_nodes$Type <- ifelse(
  cytoscape_nodes$Node == "ASCRP000064",
  "circRNA",
  ifelse(grepl("miR", cytoscape_nodes$Node), "miRNA", "mRNA")
)

write.csv(
  cytoscape_nodes,
  "07_GSE43502_mRNA_results/Cytoscape_ceRNA_network_nodes.csv",
  row.names = FALSE
)

# Summary
sink("07_GSE43502_mRNA_results/README_GSE43502_ceRNA_summary.txt")

cat("GSE43502 recurrence-associated ceRNA network summary\n\n")
cat("Comparison: Recurrence vs No recurrence\n")
cat("Criteria used: P.Value < 0.05 and |logFC| > 0.5\n\n")
cat("Final recurrence-associated ceRNA genes:", nrow(final_ceRNA_genes), "\n\n")

cat("miR-590-5p genes:\n")
cat(final_ceRNA_genes$Gene[final_ceRNA_genes$miRNA == "hsa-miR-590-5p"], sep = ", ")

cat("\n\nmiR-182-5p genes:\n")
cat(final_ceRNA_genes$Gene[final_ceRNA_genes$miRNA == "hsa-miR-182-5p"], sep = ", ")

sink()
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
genes <- unique(final_ceRNA_genes$Gene)

genes
length(genes)
gene.df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene.df
ego <- enrichGO(
  gene          = gene.df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(ego)
ekegg <- enrichKEGG(
  gene         = gene.df$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

head(ekegg)
kegg_df <- as.data.frame(ekegg)

write.csv(
  kegg_df,
  "07_GSE43502_mRNA_results/KEGG_enrichment.csv",
  row.names = FALSE
)
ego_df <- as.data.frame(ego)

write.csv(
  ego_df,
  "07_GSE43502_mRNA_results/GO_enrichment.csv",
  row.names = FALSE
)
library(enrichplot)
library(ggplot2)

# Convert results
ego_df <- as.data.frame(ego)
kegg_df <- as.data.frame(ekegg)

# Save tables
write.csv(ego_df, "07_GSE43502_mRNA_results/GO_enrichment.csv", row.names = FALSE)
write.csv(kegg_df, "07_GSE43502_mRNA_results/KEGG_enrichment.csv", row.names = FALSE)

# KEGG dotplot
p_kegg <- dotplot(ekegg, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment of Final ceRNA Target Genes") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

p_kegg

ggsave(
  "07_GSE43502_mRNA_results/Figure_KEGG_dotplot.png",
  plot = p_kegg,
  width = 8,
  height = 5,
  dpi = 600
)

ggsave(
  "07_GSE43502_mRNA_results/Figure_KEGG_dotplot.tiff",
  plot = p_kegg,
  width = 8,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

# GO dotplot, only if GO has results
if (nrow(ego_df) > 0) {
  p_go <- dotplot(ego, showCategory = 10) +
    ggtitle("GO Enrichment of Final ceRNA Target Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  p_go
  
  ggsave(
    "07_GSE43502_mRNA_results/Figure_GO_dotplot.png",
    plot = p_go,
    width = 8,
    height = 5,
    dpi = 600
  )
  
  ggsave(
    "07_GSE43502_mRNA_results/Figure_GO_dotplot.tiff",
    plot = p_go,
    width = 8,
    height = 5,
    dpi = 600,
    compression = "lzw"
  )
}
############################################################
# Check prognosis values in phenotype data
############################################################

table(pheno$`prognosis:ch1`, useNA = "ifany")
############################################################
# Define recurrence groups
############################################################

group_mrna <- factor(
  ifelse(
    pheno$`prognosis:ch1` == "recurrence",
    "Recurrence",
    ifelse(
      pheno$`prognosis:ch1` == "non-recurrence",
      "No_recurrence",
      NA
    )
  ),
  levels = c("No_recurrence", "Recurrence")
)

# Check number of samples in each group
table(group_mrna, useNA = "ifany")
############################################################
# Create sample information table
############################################################

sample_info <- data.frame(
  SampleID = rownames(pheno),
  Title = pheno$title,
  Disease_state = pheno$`disease state:ch1`,
  Prognosis = pheno$`prognosis:ch1`,
  Group_used = group_mrna
)

# Keep only samples with defined recurrence status
used_samples <- sample_info[!is.na(sample_info$Group_used), ]

# View samples used in the analysis
used_samples

# Count samples per group
table(used_samples$Group_used)

# Save sample table
write.csv(
  used_samples,
  "07_GSE43502_mRNA_results/Samples_used_for_recurrence_analysis.csv",
  row.names = FALSE
)
############################################################
# Create sample information table
############################################################

sample_info <- data.frame(
  SampleID = rownames(pheno),
  Title = pheno$title,
  Disease_state = pheno$`disease state:ch1`,
  Prognosis = pheno$`prognosis:ch1`,
  Group_used = group_mrna
)

# Keep only samples with defined recurrence status
used_samples <- sample_info[!is.na(sample_info$Group_used), ]

# View samples used in the analysis
used_samples

# Count samples per group
table(used_samples$Group_used)

# Save sample table
write.csv(
  used_samples,
  "07_GSE43502_mRNA_results/Samples_used_for_recurrence_analysis.csv",
  row.names = FALSE
)
############################################################
# Clean prognosis values
############################################################

prognosis_clean <- trimws(tolower(pheno$`prognosis:ch1`))

table(prognosis_clean, useNA = "ifany")
unique(prognosis_clean)
############################################################
# Define recurrence groups using cleaned prognosis values
############################################################

group_mrna <- factor(
  ifelse(
    prognosis_clean %in% c("recurrence", "recurrent", "relapse", "relapsed"),
    "Recurrence",
    ifelse(
      prognosis_clean %in% c("non-recurrence", "no recurrence", "non recurrence", 
                             "no_recurrence", "without recurrence", "disease-free"),
      "No_recurrence",
      NA
    )
  ),
  levels = c("No_recurrence", "Recurrence")
)

table(group_mrna, useNA = "ifany")
############################################################
# Search for phenotype columns related to recurrence/prognosis
############################################################

colnames(pheno)

grep("rec|relapse|prognosis|status|outcome|disease", 
     colnames(pheno), 
     ignore.case = TRUE, 
     value = TRUE)
############################################################
# Check values in status column
############################################################

table(pheno$status, useNA = "ifany")

unique(pheno$status)
############################################################
# Check values in Group column
############################################################

table(pheno$Group, useNA = "ifany")

unique(pheno$Group)
############################################################
# Check all ExpressionSets returned by GEO
############################################################

length(gse)

for (i in seq_along(gse)) {
  cat("\n==============================\n")
  cat("ExpressionSet number:", i, "\n")
  cat("Platform:", annotation(gse[[i]]), "\n")
  cat("Expression dimensions:\n")
  print(dim(exprs(gse[[i]])))
  cat("Phenotype dimensions:\n")
  print(dim(pData(gse[[i]])))
  cat("Phenotype columns:\n")
  print(colnames(pData(gse[[i]])))
}
############################################################
# Check all ExpressionSets returned by GEO
############################################################

library(GEOquery)
library(Biobase)

# Check number of ExpressionSets
length(gse)

for (i in seq_along(gse)) {
  cat("\n==============================\n")
  cat("ExpressionSet number:", i, "\n")
  
  cat("Platform:\n")
  print(annotation(gse[[i]]))
  
  cat("Expression dimensions:\n")
  print(dim(exprs(gse[[i]])))
  
  cat("Phenotype dimensions:\n")
  print(dim(pData(gse[[i]])))
  
  cat("Phenotype columns:\n")
  print(colnames(pData(gse[[i]])))
}
############################################################
# Download full GEO object and inspect all GSM samples
############################################################

library(GEOquery)

# Download full GEO record, not only the expression matrix
gse_full <- getGEO("GSE43502", GSEMatrix = FALSE)

# Extract GSM samples
gsm_list <- GSMList(gse_full)

# Number of GSM samples
length(gsm_list)

# Extract metadata from all GSM samples
gsm_meta <- lapply(names(gsm_list), function(gsm_id) {
  
  gsm <- gsm_list[[gsm_id]]
  meta <- Meta(gsm)
  
  data.frame(
    GSM = gsm_id,
    Title = paste(meta$title, collapse = "; "),
    Source = paste(meta$source_name_ch1, collapse = "; "),
    Characteristics = paste(meta$characteristics_ch1, collapse = " | "),
    stringsAsFactors = FALSE
  )
})

gsm_meta <- do.call(rbind, gsm_meta)

# View metadata
head(gsm_meta)
dim(gsm_meta)

# Save all GSM metadata
write.csv(
  gsm_meta,
  "07_GSE43502_mRNA_results/GSE43502_all_GSM_metadata.csv",
  row.names = FALSE
)
############################################################
# Create a stable download folder
############################################################

dir.create("07_GSE43502_mRNA_results/GEO_full",
           recursive = TRUE,
           showWarnings = FALSE)
############################################################
# Set download method for Windows
############################################################

options(download.file.method.GEOquery = "libcurl")
options(timeout = 600)
############################################################
# Download full GEO record
############################################################

library(GEOquery)

gse_full <- getGEO(
  "GSE43502",
  GSEMatrix = FALSE,
  destdir = "07_GSE43502_mRNA_results/GEO_full"
)
############################################################
# Extract GSM samples from full GEO object
############################################################

gsm_list <- GSMList(gse_full)

# Check number of GSM samples
length(gsm_list)
############################################################
# Extract metadata from all GSM samples
############################################################

gsm_meta <- lapply(names(gsm_list), function(gsm_id) {
  
  gsm <- gsm_list[[gsm_id]]
  meta <- Meta(gsm)
  
  data.frame(
    GSM = gsm_id,
    Title = paste(meta$title, collapse = "; "),
    Source = paste(meta$source_name_ch1, collapse = "; "),
    Characteristics = paste(meta$characteristics_ch1, collapse = " | "),
    stringsAsFactors = FALSE
  )
})

gsm_meta <- do.call(rbind, gsm_meta)

# Check metadata table
dim(gsm_meta)
head(gsm_meta)

# Save metadata
write.csv(
  gsm_meta,
  "07_GSE43502_mRNA_results/GSE43502_all_GSM_metadata.csv",
  row.names = FALSE
)
############################################################
# Search prognosis-related terms in GSM metadata
############################################################

all_text <- apply(gsm_meta, 1, paste, collapse = " ")

hits <- grep(
  "prognosis|recurrence|relapse|response|resistant|chemoresistant|residual|survival|poor|good|pcr|non-pcr",
  all_text,
  ignore.case = TRUE
)

gsm_hits <- gsm_meta[hits, ]

dim(gsm_hits)
gsm_hits
############################################################
# Extract prognosis information from GSM metadata
############################################################

# Clean characteristics text
gsm_meta$Characteristics_clean <- trimws(tolower(gsm_meta$Characteristics))

# Define prognosis group
gsm_meta$Prognosis <- ifelse(
  grepl("prognosis: no recurrence", gsm_meta$Characteristics_clean),
  "No_recurrence",
  ifelse(
    grepl("prognosis: recurrence", gsm_meta$Characteristics_clean),
    "Recurrence",
    NA
  )
)

# Check prognosis groups
table(gsm_meta$Prognosis, useNA = "ifany")

# Keep only samples with prognosis information
prognosis_meta <- gsm_meta[!is.na(gsm_meta$Prognosis), ]

# Check final prognosis samples
dim(prognosis_meta)
prognosis_meta[, c("GSM", "Title", "Prognosis", "Characteristics")]
############################################################
# Save prognosis sample metadata
############################################################

write.csv(
  prognosis_meta,
  "07_GSE43502_mRNA_results/GSE43502_prognosis_samples_metadata.csv",
  row.names = FALSE
)

# Separate recurrence and no-recurrence samples
recurrence_samples <- prognosis_meta[
  prognosis_meta$Prognosis == "Recurrence",
]

no_recurrence_samples <- prognosis_meta[
  prognosis_meta$Prognosis == "No_recurrence",
]

write.csv(
  recurrence_samples,
  "07_GSE43502_mRNA_results/GSE43502_recurrence_samples.csv",
  row.names = FALSE
)

write.csv(
  no_recurrence_samples,
  "07_GSE43502_mRNA_results/GSE43502_no_recurrence_samples.csv",
  row.names = FALSE
)
############################################################
# Extract expression matrix from GSM tables
############################################################

# Get GSM IDs with prognosis information
selected_gsms <- prognosis_meta$GSM

# Extract expression tables
expr_list <- lapply(selected_gsms, function(gsm_id) {
  
  gsm <- gsm_list[[gsm_id]]
  tab <- Table(gsm)
  
  # Check required columns
  if (!all(c("ID_REF", "VALUE") %in% colnames(tab))) {
    stop(paste("ID_REF or VALUE column not found in", gsm_id))
  }
  
  values <- tab$VALUE
  names(values) <- tab$ID_REF
  
  return(values)
})

# Combine into expression matrix
expr_prognosis <- do.call(cbind, expr_list)

# Assign sample names
colnames(expr_prognosis) <- selected_gsms

# Check dimensions
dim(expr_prognosis)

# Check first rows
expr_prognosis[1:5, 1:5]
############################################################
# Match phenotype data with expression matrix
############################################################

pheno_prognosis <- prognosis_meta[
  match(colnames(expr_prognosis), prognosis_meta$GSM),
]

# Check matching
all(pheno_prognosis$GSM == colnames(expr_prognosis))

# Define final group factor
group_mrna <- factor(
  pheno_prognosis$Prognosis,
  levels = c("No_recurrence", "Recurrence")
)

# Check groups
table(group_mrna)
############################################################
# Differential expression analysis
# Comparison: Recurrence vs No recurrence
############################################################

library(limma)

# Convert expression values to numeric matrix
expr_prognosis <- apply(expr_prognosis, 2, as.numeric)
rownames(expr_prognosis) <- rownames(do.call(cbind, expr_list))

# Optional log2 transformation check
qx <- quantile(expr_prognosis, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)

if (qx[5] > 100) {
  expr_prognosis <- log2(expr_prognosis + 1)
}

# Design matrix
design <- model.matrix(~ 0 + group_mrna)
colnames(design) <- levels(group_mrna)

# Contrast
contrast_matrix <- makeContrasts(
  Recurrence_vs_NoRecurrence = Recurrence - No_recurrence,
  levels = design
)

# Fit model
fit <- lmFit(expr_prognosis, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
mrna_results <- topTable(
  fit2,
  coef = "Recurrence_vs_NoRecurrence",
  number = Inf,
  adjust.method = "BH"
)

# View results
head(mrna_results)
dim(mrna_results)
############################################################
# Save recurrence differential expression results
############################################################

mrna_results$ProbeID <- rownames(mrna_results)

write.csv(
  mrna_results,
  "07_GSE43502_mRNA_results/GSE43502_Recurrence_vs_NoRecurrence_limma_results.csv",
  row.names = FALSE
)
############################################################
# Final quality checks
############################################################

dim(expr_prognosis)
dim(pheno_prognosis)
table(group_mrna)

head(pheno_prognosis[, c("GSM", "Prognosis")])
############################################################
# Create folders for final GSE43502 recurrence results
############################################################

dir.create(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_objects",
  recursive = TRUE,
  showWarnings = FALSE
)
############################################################
# Save sample and phenotype information
############################################################

write.csv(
  pheno_prognosis,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_prognosis_sample_metadata.csv",
  row.names = FALSE
)

sample_summary <- data.frame(
  Group = names(table(group_mrna)),
  Sample_number = as.numeric(table(group_mrna))
)

write.csv(
  sample_summary,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_sample_group_summary.csv",
  row.names = FALSE
)

recurrence_samples <- pheno_prognosis[
  pheno_prognosis$Prognosis == "Recurrence",
]

no_recurrence_samples <- pheno_prognosis[
  pheno_prognosis$Prognosis == "No_recurrence",
]

write.csv(
  recurrence_samples,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_recurrence_samples.csv",
  row.names = FALSE
)

write.csv(
  no_recurrence_samples,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_no_recurrence_samples.csv",
  row.names = FALSE
)
############################################################
# Save expression matrix used in recurrence analysis
############################################################

expr_prognosis_df <- as.data.frame(expr_prognosis)

expr_prognosis_df$ProbeID <- rownames(expr_prognosis_df)

expr_prognosis_df <- expr_prognosis_df[
  , c("ProbeID", setdiff(colnames(expr_prognosis_df), "ProbeID"))
]

write.csv(
  expr_prognosis_df,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_expression_matrix_25_prognosis_samples.csv",
  row.names = FALSE
)
############################################################
# Save all limma differential expression results
############################################################

mrna_results$ProbeID <- rownames(mrna_results)

write.csv(
  mrna_results,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_Recurrence_vs_NoRecurrence_all_limma_results.csv",
  row.names = FALSE
)
############################################################
# Save significant recurrence-associated probes
############################################################

sig_mrna_adj <- mrna_results[
  mrna_results$adj.P.Val < 0.05 & abs(mrna_results$logFC) > 1,
]

sig_mrna_p <- mrna_results[
  mrna_results$P.Value < 0.05 & abs(mrna_results$logFC) > 0.5,
]

write.csv(
  sig_mrna_adj,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_significant_adjP005_logFC1.csv",
  row.names = FALSE
)

write.csv(
  sig_mrna_p,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_significant_P005_logFC05.csv",
  row.names = FALSE
)
############################################################
# Save R objects for future analysis
############################################################

saveRDS(
  expr_prognosis,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_objects/expr_prognosis.rds"
)

saveRDS(
  pheno_prognosis,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_objects/pheno_prognosis.rds"
)

saveRDS(
  group_mrna,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_objects/group_mrna.rds"
)

saveRDS(
  mrna_results,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_objects/mrna_results_limma.rds"
)

saveRDS(
  fit2,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_objects/limma_fit2.rds"
)
############################################################
# Save analysis summary README
############################################################

sink(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/README_GSE43502_recurrence_analysis_summary.txt"
)

cat("GSE43502 recurrence analysis summary\n")
cat("====================================\n\n")

cat("Dataset: GSE43502\n")
cat("Cancer type: Triple-negative breast cancer (TNBC)\n")
cat("Clinical setting: Chemoresistant / residual disease samples after neoadjuvant chemotherapy\n")
cat("Comparison: Recurrence vs No recurrence\n\n")

cat("Number of samples:\n")
print(table(group_mrna))

cat("\nExpression matrix dimensions:\n")
print(dim(expr_prognosis))

cat("\nPhenotype table dimensions:\n")
print(dim(pheno_prognosis))

cat("\nDifferential expression method: limma\n")
cat("Contrast: Recurrence - No_recurrence\n")
cat("Multiple testing correction: Benjamini-Hochberg\n\n")

cat("Number of all tested probes:", nrow(mrna_results), "\n")
cat("Significant probes using adj.P.Val < 0.05 and |logFC| > 1:", nrow(sig_mrna_adj), "\n")
cat("Significant probes using P.Value < 0.05 and |logFC| > 0.5:", nrow(sig_mrna_p), "\n\n")

cat("Output files saved in:\n")
cat("07_GSE43502_mRNA_results/Final_Recurrence_Results/\n")

sink()
############################################################
# Check all saved files
############################################################

list.files(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results",
  recursive = TRUE,
  full.names = TRUE
)
############################################################
# Define the 14 candidate genes
############################################################

candidate_genes <- c(
  "FOXO3",
  "PDCD4",
  "RECK",
  "BCL2",
  "PPM1L",
  "CBX4",
  "STAT3",
  "KDM5A",
  "CADM1",
  "CYLD",
  "THBS1",
  "HOXA9",
  "CASP2",
  "LSM14A"
)
############################################################
# Add gene symbols to limma results
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  BiocManager::install("hgu133plus2.db")
}

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}

library(hgu133plus2.db)
library(AnnotationDbi)

# Add ProbeID
mrna_results$ProbeID <- rownames(mrna_results)

# Map ProbeID to Gene Symbol
mrna_results$Gene <- mapIds(
  hgu133plus2.db,
  keys = mrna_results$ProbeID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

head(mrna_results)
############################################################
# Extract the 14 candidate genes from recurrence results
############################################################

candidate_results <- mrna_results[
  mrna_results$Gene %in% candidate_genes,
]

candidate_results <- candidate_results[
  order(candidate_results$Gene, candidate_results$P.Value),
]

candidate_results[
  , c("Gene", "ProbeID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
]
############################################################
# Keep the best probe for each candidate gene
############################################################

candidate_results_best <- candidate_results[
  order(candidate_results$Gene, candidate_results$P.Value),
]

candidate_results_best <- candidate_results_best[
  !duplicated(candidate_results_best$Gene),
]

candidate_results_best <- candidate_results_best[
  order(candidate_results_best$logFC),
]

candidate_results_best[
  , c("Gene", "ProbeID", "logFC", "P.Value", "adj.P.Val")
]
############################################################
# Extract downregulated candidate genes
############################################################

down_candidate_genes <- candidate_results_best[
  candidate_results_best$logFC < 0,
]

down_candidate_genes[
  , c("Gene", "ProbeID", "logFC", "P.Value", "adj.P.Val")
]
############################################################
# Extract significantly downregulated candidate genes
############################################################

significant_down_candidate_genes <- candidate_results_best[
  candidate_results_best$logFC < 0 &
    candidate_results_best$P.Value < 0.05,
]

significant_down_candidate_genes[
  , c("Gene", "ProbeID", "logFC", "P.Value", "adj.P.Val")
]
############################################################
# Extract strongly downregulated candidate genes
############################################################

strong_down_candidate_genes <- candidate_results_best[
  candidate_results_best$logFC < -0.5 &
    candidate_results_best$P.Value < 0.05,
]

strong_down_candidate_genes[
  , c("Gene", "ProbeID", "logFC", "P.Value", "adj.P.Val")
]
############################################################
# Save candidate gene validation results
############################################################

write.csv(
  candidate_results_best,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_14_candidate_genes_recurrence_validation.csv",
  row.names = FALSE
)

write.csv(
  down_candidate_genes,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_14_candidate_genes_downregulated_in_recurrence.csv",
  row.names = FALSE
)

write.csv(
  significant_down_candidate_genes,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_14_candidate_genes_significant_downregulated_P005.csv",
  row.names = FALSE
)

write.csv(
  strong_down_candidate_genes,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_14_candidate_genes_strong_downregulated_logFC05_P005.csv",
  row.names = FALSE
)
############################################################
# Check missing candidate genes
############################################################

found_genes <- unique(candidate_results_best$Gene)

missing_genes <- setdiff(candidate_genes, found_genes)

missing_genes
length(missing_genes)

write.csv(
  data.frame(Missing_gene = missing_genes),
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_missing_candidate_genes.csv",
  row.names = FALSE
)
############################################################
# Show final candidate gene validation table
############################################################

candidate_results_best[
  , c("Gene", "ProbeID", "logFC", "Direction",
      "P.Value", "adj.P.Val",
      "Significant_P005", "Significant_adjP005")
]
############################################################
# Show all 14 candidate genes without filtering
############################################################

candidate_results_best[
  , c("Gene", "ProbeID", "logFC", "P.Value", "adj.P.Val", "Direction")
]
############################################################
# Fix candidate results table and rank by P.Value
############################################################

# Add ProbeID if missing
if (!"ProbeID" %in% colnames(candidate_results_best)) {
  candidate_results_best$ProbeID <- rownames(candidate_results_best)
}

# Add Direction if missing
candidate_results_best$Direction <- ifelse(
  candidate_results_best$logFC < 0,
  "Downregulated_in_Recurrence",
  "Upregulated_in_Recurrence"
)

# Add significance by nominal P.Value
candidate_results_best$Significant_P005 <- ifelse(
  candidate_results_best$P.Value < 0.05,
  "Yes",
  "No"
)

# Add significance by adjusted P.Value
candidate_results_best$Significant_adjP005 <- ifelse(
  candidate_results_best$adj.P.Val < 0.05,
  "Yes",
  "No"
)

# Rank all 14 genes by P.Value
candidate_results_best_Pvalue <- candidate_results_best[
  order(candidate_results_best$P.Value),
]

# Show final table ranked by P.Value
candidate_results_best_Pvalue[
  , c("Gene", "ProbeID", "logFC", "Direction",
      "P.Value", "adj.P.Val",
      "Significant_P005", "Significant_adjP005")
]
############################################################
# Final candidate gene table ranked by P.Value
############################################################

candidate_results_best_Pvalue <- candidate_results_best[
  order(candidate_results_best$P.Value),
]

final_candidate_table <- candidate_results_best_Pvalue[
  , c("Gene", "ProbeID", "logFC", "Direction",
      "P.Value", "adj.P.Val",
      "Significant_P005", "Significant_adjP005")
]

final_candidate_table

write.csv(
  final_candidate_table,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_14_candidate_genes_final_ranked_by_Pvalue.csv",
  row.names = FALSE
)
############################################################
# Final significant downregulated candidate genes
############################################################

final_down_sig_table <- final_candidate_table[
  final_candidate_table$logFC < 0 &
    final_candidate_table$P.Value < 0.05,
]

final_down_sig_table <- final_down_sig_table[
  order(final_down_sig_table$P.Value),
]

final_down_sig_table

write.csv(
  final_down_sig_table,
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_candidate_genes_downregulated_P005_ranked.csv",
  row.names = FALSE
)
############################################################
# Candidate genes not significant by P.Value
############################################################

not_significant_Pvalue <- final_candidate_table[
  final
  ############################################################
  # Safer way to define miRNA-target candidate genes
  ############################################################
  
  mir590_genes <- c(
    "PDCD4",
    "TMEM170A",
    "RBPJ",
    "STAT3",
    "FRS2",
    "CBX4",
    "MBNL1",
    "KLHL15",
    "GID4"
  )
  
  mir182_genes <- c(
    "FLOT1",
    "MITF",
    "FOXF2",
    "ADCY6",
    "RARG",
    "NPTX1",
    "RASA2",
    "RECK",
    "FBXW7",
    "PPM1L",
    "FOXO3",
    "EVI5",
    "FGF9",
    "IGF1R",
    "BCL2",
    "GABRB1",
    "LSM14A",
    "CITED2",
    "CCND2",
    "PFN1",
    "USP5",
    "CREB1",
    "MTSS1",
    "HOXA9",
    "NUP50",
    "CLOCK",
    "THBS1",
    "CASP2",
    "DDAH1",
    "RCC2",
    "SATB2",
    "LRRC4",
    "TMEM170B",
    "NDRG1",
    "FOXO1",
    "EIF4EBP2",
    "PPP1R12A",
    "ACER2",
    "ARRDC3",
    "PRKAA2",
    "CHL1",
    "TP53INP1",
    "SESN2",
    "STK17B",
    "PPP1R11",
    "KDM5A",
    "QSER1",
    "CYLD",
    "NUFIP2",
    "CADM1"
  )
  
  mirna_targets <- rbind(
    data.frame(
      miRNA = "hsa-miR-590-5p",
      Gene = mir590_genes,
      stringsAsFactors = FALSE
    ),
    data.frame(
      miRNA = "hsa-miR-182-5p",
      Gene = mir182_genes,
      stringsAsFactors = FALSE
    )
  )
  
  dim(mirna_targets)
  head(mirna_targets)
  tail(mirna_targets)
  ############################################################
  # Define miRNA-target candidate genes
  ############################################################
  
  mirna_targets <- data.frame(
    miRNA = c(
      rep("hsa-miR-590-5p", 9),
      rep("hsa-miR-182-5p", 50)
    ),
    Gene = c(
      "PDCD4",
      "TMEM170A",
      "RBPJ",
      "STAT3",
      "FRS2",
      "CBX4",
      "MBNL1",
      "KLHL15",
      "GID4",
      
      "FLOT1",
      "MITF",
      "FOXF2",
      "ADCY6",
      "RARG",
      "NPTX1",
      "RASA2",
      "RECK",
      "FBXW7",
      "PPM1L",
      "FOXO3",
      "EVI5",
      "FGF9",
      "IGF1R",
      "BCL2",
      "GABRB1",
      "LSM14A",
      "CITED2",
      "CCND2",
      "PFN1",
      "USP5",
      "CREB1",
      "MTSS1",
      "HOXA9",
      "NUP50",
      "CLOCK",
      "THBS1",
      "CASP2",
      "DDAH1",
      "RCC2",
      "SATB2",
      "LRRC4",
      "TMEM170B",
      "NDRG1",
      "FOXO1",
      "EIF4EBP2",
      "PPP1R12A",
      "ACER2",
      "ARRDC3",
      "PRKAA2",
      "CHL1",
      "TP53INP1",
      "SESN2",
      "STK17B",
      "PPP1R11",
      "KDM5A",
      "QSER1",
      "CYLD",
      "NUFIP2",
      "CADM1"
    ),
    stringsAsFactors = FALSE
  )
  ############################################################
  # Add gene symbols to limma results if needed
  ############################################################
  
  if (!"Gene" %in% colnames(mrna_results)) {
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
      BiocManager::install("hgu133plus2.db")
    }
    
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      BiocManager::install("AnnotationDbi")
    }
    
    library(hgu133plus2.db)
    library(AnnotationDbi)
    
    mrna_results$ProbeID <- rownames(mrna_results)
    
    mrna_results$Gene <- mapIds(
      hgu133plus2.db,
      keys = mrna_results$ProbeID,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
  }
  
  if (!"ProbeID" %in% colnames(mrna_results)) {
    mrna_results$ProbeID <- rownames(mrna_results)
  }
  ############################################################
  # Extract miRNA target genes from GSE43502 recurrence results
  ############################################################
  
  target_results_all_probes <- mrna_results[
    mrna_results$Gene %in% unique(mirna_targets$Gene),
  ]
  
  # Add miRNA information
  target_results_all_probes <- merge(
    target_results_all_probes,
    mirna_targets,
    by = "Gene",
    all.x = TRUE
  )
  
  # Rank by Gene and P.Value
  target_results_all_probes <- target_results_all_probes[
    order(target_results_all_probes$Gene, target_results_all_probes$P.Value),
  ]
  
  # Keep best probe per gene based on lowest P.Value
  target_results_best <- target_results_all_probes[
    !duplicated(target_results_all_probes$Gene),
  ]
  
  # Classify direction
  target_results_best$Direction <- ifelse(
    target_results_best$logFC < 0,
    "Downregulated_in_Recurrence",
    "Upregulated_in_Recurrence"
  )
  
  # Nominal significance
  target_results_best$Nominal_P005 <- ifelse(
    target_results_best$P.Value < 0.05,
    "Yes",
    "No"
  )
  
  # FDR significance
  target_results_best$FDR_adjP005 <- ifelse(
    target_results_best$adj.P.Val < 0.05,
    "Yes",
    "No"
  )
  
  # Rank final table by P.Value
  target_results_best <- target_results_best[
    order(target_results_best$P.Value),
  ]
  
  target_results_best[
    , c("miRNA", "Gene", "ProbeID", "logFC", "Direction",
        "P.Value", "adj.P.Val", "Nominal_P005", "FDR_adjP005")
  ]
  ############################################################
  # Nominally significant miRNA target genes
  ############################################################
  
  target_nominal_sig <- target_results_best[
    target_results_best$P.Value < 0.05,
  ]
  
  target_nominal_sig <- target_nominal_sig[
    order(target_nominal_sig$P.Value),
  ]
  
  target_nominal_sig[
    , c("miRNA", "Gene", "ProbeID", "logFC", "Direction",
        "P.Value", "adj.P.Val")
  ]
  ############################################################
  # Downregulated nominally significant miRNA target genes
  ############################################################
  
  target_down_nominal_sig <- target_results_best[
    target_results_best$logFC < 0 &
      target_results_best$P.Value < 0.05,
  ]
  
  target_down_nominal_sig <- target_down_nominal_sig[
    order(target_down_nominal_sig$P.Value),
  ]
  
  target_down_nominal_sig[
    , c("miRNA", "Gene", "ProbeID", "logFC", "Direction",
        "P.Value", "adj.P.Val")
  ]
  ############################################################
  # Upregulated nominally significant miRNA target genes
  ############################################################
  
  target_up_nominal_sig <- target_results_best[
    target_results_best$logFC > 0 &
      target_results_best$P.Value < 0.05,
  ]
  
  target_up_nominal_sig <- target_up_nominal_sig[
    order(target_up_nominal_sig$P.Value),
  ]
  
  target_up_nominal_sig[
    , c("miRNA", "Gene", "ProbeID", "logFC", "Direction",
        "P.Value", "adj.P.Val")
  ]
  ############################################################
  # Save miRNA target validation results
  ############################################################
  
  write.csv(
    target_results_best,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_all_miRNA_target_genes_ranked_by_Pvalue.csv",
    row.names = FALSE
  )
  
  write.csv(
    target_nominal_sig,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_nominally_significant_miRNA_target_genes_P005.csv",
    row.names = FALSE
  )
  
  write.csv(
    target_down_nominal_sig,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_downregulated_nominally_significant_miRNA_targets.csv",
    row.names = FALSE
  )
  
  write.csv(
    target_up_nominal_sig,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_upregulated_nominally_significant_miRNA_targets.csv",
    row.names = FALSE
  )
  ############################################################
  # Summary counts by miRNA
  ############################################################
  
  table(target_results_best$miRNA)
  
  table(target_nominal_sig$miRNA)
  
  table(target_down_nominal_sig$miRNA)
  
  table(target_up_nominal_sig$miRNA)
  ############################################################
  # Check missing miRNA target genes
  ############################################################
  
  found_target_genes <- unique(target_results_best$Gene)
  
  missing_target_genes <- setdiff(unique(mirna_targets$Gene), found_target_genes)
  
  missing_target_genes
  length(missing_target_genes)
  
  write.csv(
    data.frame(Missing_gene = missing_target_genes),
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_missing_miRNA_target_genes.csv",
    row.names = FALSE
  )
  ############################################################
  # Check missing miRNA target genes
  ############################################################
  
  found_target_genes <- unique(target_results_best$Gene)
  
  missing_target_genes <- setdiff(unique(mirna_targets$Gene), found_target_genes)
  
  missing_target_genes
  length(missing_target_genes)
  
  write.csv(
    data.frame(Missing_gene = missing_target_genes),
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/Tables/GSE43502_missing_miRNA_target_genes.csv",
    row.names = FALSE
  )
  ############################################################
  # Define the 59 miRNA target genes
  ############################################################
  
  mir590_genes <- c(
    "PDCD4",
    "TMEM170A",
    "RBPJ",
    "STAT3",
    "FRS2",
    "CBX4",
    "MBNL1",
    "KLHL15",
    "GID4"
  )
  
  mir182_genes <- c(
    "FLOT1",
    "MITF",
    "FOXF2",
    "ADCY6",
    "RARG",
    "NPTX1",
    "RASA2",
    "RECK",
    "FBXW7",
    "PPM1L",
    "FOXO3",
    "EVI5",
    "FGF9",
    "IGF1R",
    "BCL2",
    "GABRB1",
    "LSM14A",
    "CITED2",
    "CCND2",
    "PFN1",
    "USP5",
    "CREB1",
    "MTSS1",
    "HOXA9",
    "NUP50",
    "CLOCK",
    "THBS1",
    "CASP2",
    "DDAH1",
    "RCC2",
    "SATB2",
    "LRRC4",
    "TMEM170B",
    "NDRG1",
    "FOXO1",
    "EIF4EBP2",
    "PPP1R12A",
    "ACER2",
    "ARRDC3",
    "PRKAA2",
    "CHL1",
    "TP53INP1",
    "SESN2",
    "STK17B",
    "PPP1R11",
    "KDM5A",
    "QSER1",
    "CYLD",
    "NUFIP2",
    "CADM1"
  )
  
  mirna_targets_59 <- rbind(
    data.frame(
      miRNA = "hsa-miR-590-5p",
      Gene = mir590_genes,
      stringsAsFactors = FALSE
    ),
    data.frame(
      miRNA = "hsa-miR-182-5p",
      Gene = mir182_genes,
      stringsAsFactors = FALSE
    )
  )
  
  # Check number of target genes
  dim(mirna_targets_59)
  length(unique(mirna_targets_59$Gene))
  
  write.csv(
    mirna_targets_59,
    file.path(output_59, "01_Input_59_miRNA_target_genes.csv"),
    row.names = FALSE
  )
  ############################################################
  # Start fresh: Filter 59 miRNA target genes
  ############################################################
  
  output_59 <- "07_GSE43502_mRNA_results/Final_Recurrence_Results/Filter_59_miRNA_Target_Genes"
  
  dir.create(
    output_59,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  output_59
  ############################################################
  # Define miR-590-5p and miR-182-5p target genes
  ############################################################
  
  mir590_genes <- c(
    "PDCD4",
    "TMEM170A",
    "RBPJ",
    "STAT3",
    "FRS2",
    "CBX4",
    "MBNL1",
    "KLHL15",
    "GID4"
  )
  
  mir182_genes <- c(
    "FLOT1",
    "MITF",
    "FOXF2",
    "ADCY6",
    "RARG",
    "NPTX1",
    "RASA2",
    "RECK",
    "FBXW7",
    "PPM1L",
    "FOXO3",
    "EVI5",
    "FGF9",
    "IGF1R",
    "BCL2",
    "GABRB1",
    "LSM14A",
    "CITED2",
    "CCND2",
    "PFN1",
    "USP5",
    "CREB1",
    "MTSS1",
    "HOXA9",
    "NUP50",
    "CLOCK",
    "THBS1",
    "CASP2",
    "DDAH1",
    "RCC2",
    "SATB2",
    "LRRC4",
    "TMEM170B",
    "NDRG1",
    "FOXO1",
    "EIF4EBP2",
    "PPP1R12A",
    "ACER2",
    "ARRDC3",
    "PRKAA2",
    "CHL1",
    "TP53INP1",
    "SESN2",
    "STK17B",
    "PPP1R11",
    "KDM5A",
    "QSER1",
    "CYLD",
    "NUFIP2",
    "CADM1"
  )
  
  mirna_targets_59 <- rbind(
    data.frame(
      miRNA = "hsa-miR-590-5p",
      Gene = mir590_genes,
      stringsAsFactors = FALSE
    ),
    data.frame(
      miRNA = "hsa-miR-182-5p",
      Gene = mir182_genes,
      stringsAsFactors = FALSE
    )
  )
  
  dim(mirna_targets_59)
  length(unique(mirna_targets_59$Gene))
  ############################################################
  # Save input 59 miRNA target genes
  ############################################################
  
  write.csv(
    mirna_targets_59,
    file.path(output_59, "01_Input_59_miRNA_target_genes.csv"),
    row.names = FALSE
  )
  ############################################################
  # Prepare limma results columns
  ############################################################
  
  if (!"ProbeID" %in% colnames(mrna_results)) {
    mrna_results$ProbeID <- rownames(mrna_results)
  }
  
  if (!"Gene" %in% colnames(mrna_results)) {
    
    library(hgu133plus2.db)
    library(AnnotationDbi)
    
    mrna_results$Gene <- mapIds(
      hgu133plus2.db,
      keys = mrna_results$ProbeID,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
  }
  
  colnames(mrna_results)
  ############################################################
  # Extract all probes for the 59 miRNA target genes
  ############################################################
  
  target_59_all_probes <- mrna_results[
    mrna_results$Gene %in% unique(mirna_targets_59$Gene),
  ]
  
  target_59_all_probes <- merge(
    target_59_all_probes,
    mirna_targets_59,
    by = "Gene",
    all.x = TRUE
  )
  
  target_59_all_probes$Direction <- ifelse(
    target_59_all_probes$logFC < 0,
    "Downregulated_in_Recurrence",
    "Upregulated_in_Recurrence"
  )
  
  target_59_all_probes$Nominal_P005 <- ifelse(
    target_59_all_probes$P.Value < 0.05,
    "Yes",
    "No"
  )
  
  target_59_all_probes$FDR_adjP005 <- ifelse(
    target_59_all_probes$adj.P.Val < 0.05,
    "Yes",
    "No"
  )
  
  target_59_all_probes <- target_59_all_probes[
    order(target_59_all_probes$miRNA,
          target_59_all_probes$Gene,
          target_59_all_probes$P.Value),
  ]
  
  dim(target_59_all_probes)
  
  write.csv(
    target_59_all_probes,
    file.path(output_59, "02_GSE43502_59_targets_all_probes.csv"),
    row.names = FALSE
  )
  ############################################################
  # Keep best probe per gene based on lowest P.Value
  ############################################################
  
  target_59_best <- target_59_all_probes[
    order(target_59_all_probes$Gene,
          target_59_all_probes$P.Value),
  ]
  
  target_59_best <- target_59_best[
    !duplicated(target_59_best$Gene),
  ]
  
  target_59_best <- target_59_best[
    order(target_59_best$P.Value),
  ]
  
  dim(target_59_best)
  
  write.csv(
    target_59_best,
    file.path(output_59, "03_GSE43502_59_targets_best_probe_ranked_by_Pvalue.csv"),
    row.names = FALSE
  )
  ############################################################
  # Save missing genes
  ############################################################
  
  found_59_genes <- unique(target_59_best$Gene)
  
  missing_59_genes <- setdiff(
    unique(mirna_targets_59$Gene),
    found_59_genes
  )
  
  missing_59_genes
  length(missing_59_genes)
  
  write.csv(
    data.frame(Missing_gene = missing_59_genes),
    file.path(output_59, "04_GSE43502_59_targets_missing_genes.csv"),
    row.names = FALSE
  )
  ############################################################
  # Filter by nominal P.Value < 0.05
  ############################################################
  
  target_59_nominal_P005 <- target_59_best[
    target_59_best$P.Value < 0.05,
  ]
  
  target_59_nominal_P005 <- target_59_nominal_P005[
    order(target_59_nominal_P005$P.Value),
  ]
  
  write.csv(
    target_59_nominal_P005,
    file.path(output_59, "05_GSE43502_59_targets_nominally_significant_P005.csv"),
    row.names = FALSE
  )
  ############################################################
  # Downregulated genes with P.Value < 0.05
  ############################################################
  
  target_59_down_P005 <- target_59_best[
    target_59_best$logFC < 0 &
      target_59_best$P.Value < 0.05,
  ]
  
  target_59_down_P005 <- target_59_down_P005[
    order(target_59_down_P005$P.Value),
  ]
  
  write.csv(
    target_59_down_P005,
    file.path(output_59, "06_GSE43502_59_targets_downregulated_P005.csv"),
    row.names = FALSE
  )
  ############################################################
  # Upregulated genes with P.Value < 0.05
  ############################################################
  
  target_59_up_P005 <- target_59_best[
    target_59_best$logFC > 0 &
      target_59_best$P.Value < 0.05,
  ]
  
  target_59_up_P005 <- target_59_up_P005[
    order(target_59_up_P005$P.Value),
  ]
  
  write.csv(
    target_59_up_P005,
    file.path(output_59, "07_GSE43502_59_targets_upregulated_P005.csv"),
    row.names = FALSE
  )
  ############################################################
  # Save filtering summary
  ############################################################
  
  summary_59 <- data.frame(
    Item = c(
      "Input miRNA-target pairs",
      "Unique input genes",
      "Genes found in GSE43502",
      "Missing genes",
      "Nominally significant genes P.Value < 0.05",
      "Downregulated significant genes",
      "Upregulated significant genes",
      "FDR significant genes adj.P.Val < 0.05"
    ),
    Count = c(
      nrow(mirna_targets_59),
      length(unique(mirna_targets_59$Gene)),
      length(found_59_genes),
      length(missing_59_genes),
      nrow(target_59_nominal_P005),
      nrow(target_59_down_P005),
      nrow(target_59_up_P005),
      sum(target_59_best$adj.P.Val < 0.05)
    )
  )
  
  summary_59
  
  write.csv(
    summary_59,
    file.path(output_59, "08_GSE43502_59_targets_filtering_summary.csv"),
    row.names = FALSE
  )
  ############################################################
  # Check saved files
  ############################################################
  
  list.files(
    output_59,
    full.names = TRUE
  )
  ############################################################
  # Save current GSE43502 recurrence validation workspace
  ############################################################
  
  dir.create(
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  save.image(
    file = "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/GSE43502_recurrence_validation_workspace.RData"
  )
  ############################################################
  # Save important objects separately
  ############################################################
  
  saveRDS(
    mrna_results,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/mrna_results_GSE43502_recurrence.rds"
  )
  
  saveRDS(
    target_59_best,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/target_59_best_ranked_by_Pvalue.rds"
  )
  
  saveRDS(
    target_59_down_P005,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/target_59_downregulated_P005.rds"
  )
  
  saveRDS(
    pheno_prognosis,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/pheno_prognosis_25_samples.rds"
  )
  
  saveRDS(
    expr_prognosis,
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/expr_prognosis_25_samples.rds"
  )
  ############################################################
  # Check saved R workspace files
  ############################################################
  
  list.files(
    "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace",
    full.names = TRUE
  )