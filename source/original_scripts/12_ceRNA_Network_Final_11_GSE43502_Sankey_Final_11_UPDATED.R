############################################################
# Final ceRNA network using 11 validated downregulated genes
# GSE43502 recurrence validation
############################################################

# Set working directory
setwd("D:/breast cancer/breast cancer")

# Load final validated GSE43502 recurrence results
target_59_down_P005 <- readRDS(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/R_workspace/target_59_downregulated_P005.rds"
)

# Define final 11 genes for network
final_11_genes <- c(
  "FOXO3",
  "PDCD4",
  "RECK",
  "BCL2",
  "PPM1L",
  "KDM5A",
  "CADM1",
  "CYLD",
  "HOXA9",
  "CASP2",
  "LSM14A"
)

# Extract final 11 genes
final_11_network_genes <- target_59_down_P005[
  target_59_down_P005$Gene %in% final_11_genes,
]

final_11_network_genes <- final_11_network_genes[
  order(final_11_network_genes$P.Value),
]

final_11_network_genes
############################################################
# Alternative: load final 11 genes from CSV
############################################################

final_11_network_genes <- read.csv(
  "07_GSE43502_mRNA_results/Final_Recurrence_Results/Filter_59_miRNA_Target_Genes/06_GSE43502_59_targets_downregulated_P005.csv"
)

final_11_genes <- c(
  "FOXO3",
  "PDCD4",
  "RECK",
  "BCL2",
  "PPM1L",
  "KDM5A",
  "CADM1",
  "CYLD",
  "HOXA9",
  "CASP2",
  "LSM14A"
)

final_11_network_genes <- final_11_network_genes[
  final_11_network_genes$Gene %in% final_11_genes,
]

final_11_network_genes
############################################################
# Find all target_59 files
############################################################

list.files(
  "07_GSE43502_mRNA_results",
  pattern = "target|59|down",
  recursive = TRUE,
  full.names = TRUE
)
############################################################
# Final ceRNA network using 11 validated downregulated genes
# Manual clean version
############################################################

setwd("D:/breast cancer/breast cancer")

# Create final network folder
output_network <- "12_ceRNA_Network_Final_11_GSE43502"

dir.create(
  output_network,
  recursive = TRUE,
  showWarnings = FALSE
)

# Define final 11 validated downregulated genes
final_11_genes <- c(
  "FOXO3",
  "PDCD4",
  "RECK",
  "BCL2",
  "PPM1L",
  "KDM5A",
  "CADM1",
  "CYLD",
  "HOXA9",
  "CASP2",
  "LSM14A"
)

# Define miRNA-gene relationships based on target filtering
final_11_network_genes <- data.frame(
  Gene = c(
    "FOXO3",
    "PDCD4",
    "RECK",
    "BCL2",
    "PPM1L",
    "KDM5A",
    "CADM1",
    "CYLD",
    "HOXA9",
    "CASP2",
    "LSM14A"
  ),
  miRNA = c(
    "hsa-miR-182-5p",
    "hsa-miR-590-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p"
  ),
  Direction = "Downregulated_in_Recurrence",
  stringsAsFactors = FALSE
)

# Save final 11 gene table
write.csv(
  final_11_network_genes,
  file.path(output_network, "Final_11_downregulated_genes_for_network.csv"),
  row.names = FALSE
)

final_11_network_genes
############################################################
# Final ceRNA network after removing upregulated genes
# circRNA: hsa_circRNA_000554
############################################################

setwd("D:/breast cancer/breast cancer")

output_network <- "12_ceRNA_Network_Final_11_GSE43502"

dir.create(
  output_network,
  recursive = TRUE,
  showWarnings = FALSE
)

# Correct circRNA name
circRNA_id <- "hsa_circRNA_000554"

# Final downregulated mRNAs only
final_nodes <- data.frame(
  name = c(
    "hsa_circRNA_000554",
    "hsa-miR-182-5p",
    "hsa-miR-590-5p",
    "FOXO3",
    "RECK",
    "BCL2",
    "PPM1L",
    "KDM5A",
    "CADM1",
    "CYLD",
    "HOXA9",
    "CASP2",
    "LSM14A",
    "PDCD4"
  ),
  type = c(
    "circRNA",
    "miRNA",
    "miRNA",
    rep("mRNA", 11)
  ),
  x = c(
    0,
    0,
    0,
    -3.262708033,
    -2.836838713,
    -2.265553408,
    -1.578136128,
    -0.809823796,
    -6.61E-16,
    0.809823796,
    2.265553408,
    2.836838713,
    3.262708033,
    2.539371337
  ),
  y = c(
    0,
    -2,
    2,
    -3.521425742,
    -4.216381311,
    -4.797725461,
    -5.235658567,
    -5.507732233,
    -5.6,
    -5.507732233,
    -4.797725461,
    -4.216381311,
    -3.521425742,
    3.778086953
  ),
  stringsAsFactors = FALSE
)

# Save node file
write.csv(
  final_nodes,
  file.path(output_network, "Cytoscape_nodes_final_downregulated_ceRNA_network.csv"),
  row.names = FALSE
)

final_nodes
############################################################
# Create final ceRNA edge file
############################################################

# circRNA-miRNA interactions
edges_circ_mirna <- data.frame(
  source = c(
    "hsa_circRNA_000554",
    "hsa_circRNA_000554"
  ),
  target = c(
    "hsa-miR-182-5p",
    "hsa-miR-590-5p"
  ),
  interaction = "circRNA-miRNA",
  stringsAsFactors = FALSE
)

# miRNA-mRNA interactions
edges_mirna_mrna <- data.frame(
  source = c(
    rep("hsa-miR-182-5p", 10),
    "hsa-miR-590-5p"
  ),
  target = c(
    "FOXO3",
    "RECK",
    "BCL2",
    "PPM1L",
    "KDM5A",
    "CADM1",
    "CYLD",
    "HOXA9",
    "CASP2",
    "LSM14A",
    "PDCD4"
  ),
  interaction = "miRNA-mRNA",
  stringsAsFactors = FALSE
)

final_edges <- rbind(
  edges_circ_mirna,
  edges_mirna_mrna
)

write.csv(
  final_edges,
  file.path(output_network, "Cytoscape_edges_final_downregulated_ceRNA_network.csv"),
  row.names = FALSE
)

final_edges
############################################################
# Save final downregulated mRNA list for PPI
############################################################

final_downregulated_genes <- c(
  "FOXO3",
  "RECK",
  "BCL2",
  "PPM1L",
  "KDM5A",
  "CADM1",
  "CYLD",
  "HOXA9",
  "CASP2",
  "LSM14A",
  "PDCD4"
)

write.table(
  final_downregulated_genes,
  file.path(output_network, "STRING_PPI_final_downregulated_11_genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
list.files(
  output_network,
  full.names = TRUE
)
############################################################
# Load required packages
############################################################

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}

if (!requireNamespace("tidygraph", quietly = TRUE)) {
  install.packages("tidygraph")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(openxlsx)
############################################################
# Read final ceRNA network files
############################################################

setwd("D:/breast cancer/breast cancer")

output_network <- "12_ceRNA_Network_Final_11_GSE43502"

final_edges <- read.csv(
  file.path(output_network, "Cytoscape_edges_final_downregulated_ceRNA_network.csv"),
  stringsAsFactors = FALSE
)

final_nodes <- read.csv(
  file.path(output_network, "Cytoscape_nodes_final_downregulated_ceRNA_network.csv"),
  stringsAsFactors = FALSE
)

final_genes <- read.csv(
  file.path(output_network, "Final_11_downregulated_genes_for_network.csv"),
  stringsAsFactors = FALSE
)

final_edges
final_nodes
final_genes
############################################################
# Draw final ceRNA network
############################################################

# Build graph
g <- graph_from_data_frame(
  d = final_edges,
  vertices = final_nodes,
  directed = TRUE
)

# Convert to tidygraph
tg <- as_tbl_graph(g)

# Node colors
node_colors <- c(
  "circRNA" = "#E64B35",
  "miRNA"   = "#4DBBD5",
  "mRNA"    = "#00A087"
)

# Node shapes
node_shapes <- c(
  "circRNA" = 21,
  "miRNA"   = 22,
  "mRNA"    = 24
)

# Plot
p_ceRNA <- ggraph(tg, layout = "manual", x = final_nodes$x, y = final_nodes$y) +
  geom_edge_link(
    aes(label = interaction),
    arrow = arrow(length = unit(3, "mm")),
    end_cap = circle(4, "mm"),
    color = "grey55",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(fill = type, shape = type),
    size = 8,
    color = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = node_colors) +
  scale_shape_manual(values = node_shapes) +
  theme_void() +
  ggtitle("Final hsa_circRNA_000554-centered ceRNA Network") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_ceRNA
############################################################
# Save ceRNA network figure
############################################################

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_11_downregulated_genes.png"),
  plot = p_ceRNA,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_11_downregulated_genes.pdf"),
  plot = p_ceRNA,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_11_downregulated_genes.tiff"),
  plot = p_ceRNA,
  width = 10,
  height = 8,
  dpi = 600,
  compression = "lzw"
)
############################################################
# Save important results into one Excel workbook
############################################################

excel_file <- file.path(
  output_network,
  "Final_ceRNA_Network_11_downregulated_genes_GSE43502.xlsx"
)

wb <- createWorkbook()

# Sheet 1: final genes
addWorksheet(wb, "Final_11_genes")
writeData(wb, "Final_11_genes", final_genes)

# Sheet 2: Cytoscape edges
addWorksheet(wb, "Cytoscape_edges")
writeData(wb, "Cytoscape_edges", final_edges)

# Sheet 3: Cytoscape nodes
addWorksheet(wb, "Cytoscape_nodes")
writeData(wb, "Cytoscape_nodes", final_nodes)

# Sheet 4: excluded upregulated genes
excluded_upregulated_genes <- data.frame(
  Gene = c("CBX4", "STAT3", "THBS1"),
  Reason = "Excluded from final ceRNA network because logFC > 0 in Recurrence vs No_recurrence",
  Direction = "Upregulated_in_Recurrence",
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Excluded_upregulated")
writeData(wb, "Excluded_upregulated", excluded_upregulated_genes)

# Sheet 5: summary
summary_table <- data.frame(
  Item = c(
    "circRNA",
    "miRNAs",
    "Final mRNAs",
    "Excluded genes",
    "Validation dataset",
    "Comparison",
    "Selection rule"
  ),
  Description = c(
    "hsa_circRNA_000554",
    "hsa-miR-182-5p; hsa-miR-590-5p",
    "11 downregulated target genes",
    "CBX4, STAT3, THBS1",
    "GSE43502",
    "Recurrence vs No_recurrence",
    "P.Value < 0.05 and logFC < 0"
  ),
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_table)

# Format workbook
sheets <- names(wb)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAD3",
  border = "Bottom"
)

for (sh in sheets) {
  addStyle(
    wb,
    sheet = sh,
    style = header_style,
    rows = 1,
    cols = 1:50,
    gridExpand = TRUE
  )
  freezePane(wb, sheet = sh, firstRow = TRUE)
  setColWidths(wb, sheet = sh, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  excel_file,
  overwrite = TRUE
)
############################################################
# Save README summary
############################################################

sink(
  file.path(output_network, "README_final_ceRNA_network_11_downregulated_genes.txt")
)

cat("Final ceRNA network summary\n")
cat("====================================\n\n")

cat("circRNA: hsa_circRNA_000554\n")
cat("miRNAs: hsa-miR-182-5p, hsa-miR-590-5p\n")
cat("Validation dataset: GSE43502\n")
cat("Comparison: Recurrence vs No_recurrence\n")
cat("Selection criteria: P.Value < 0.05 and logFC < 0\n\n")

cat("Final 11 downregulated mRNA genes:\n")
cat(final_downregulated_genes, sep = ", ")

cat("\n\nExcluded upregulated genes:\n")
cat("CBX4, STAT3, THBS1\n")

cat("\n\nReason for exclusion:\n")
cat("These genes showed logFC > 0, indicating upregulation in recurrence samples.\n")
cat("Therefore, they were excluded from the final downregulated ceRNA network.\n\n")

cat("Saved files:\n")
print(list.files(output_network, full.names = FALSE))

sink()
############################################################
# Check all final files
############################################################

list.files(
  output_network,
  full.names = TRUE
)
############################################################
# Open final output folder
############################################################

shell.exec(
  normalizePath(output_network)
)
############################################################
# Draw final ceRNA network without edge labels
############################################################

g <- graph_from_data_frame(
  d = final_edges,
  vertices = final_nodes,
  directed = TRUE
)

tg <- as_tbl_graph(g)

node_colors <- c(
  "circRNA" = "#E64B35",
  "miRNA"   = "#4DBBD5",
  "mRNA"    = "#00A087"
)

node_shapes <- c(
  "circRNA" = 21,
  "miRNA"   = 22,
  "mRNA"    = 24
)

p_ceRNA_clean <- ggraph(tg, layout = "manual", x = final_nodes$x, y = final_nodes$y) +
  geom_edge_link(
    arrow = arrow(length = unit(3, "mm")),
    end_cap = circle(4, "mm"),
    color = "grey45",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(fill = type, shape = type),
    size = 8,
    color = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = node_colors) +
  scale_shape_manual(values = node_shapes) +
  theme_void() +
  ggtitle("Final hsa_circRNA_000554-centered ceRNA Network") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_ceRNA_clean
############################################################
# Save clean ceRNA network figure without edge labels
############################################################

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_11_downregulated_genes_clean.png"),
  plot = p_ceRNA_clean,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_11_downregulated_genes_clean.pdf"),
  plot = p_ceRNA_clean,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_11_downregulated_genes_clean.tiff"),
  plot = p_ceRNA_clean,
  width = 10,
  height = 8,
  dpi = 600,
  compression = "lzw"
)
############################################################
# Save final clean ceRNA network outputs
############################################################

# Make sure output folder is defined
output_network <- "12_ceRNA_Network_Final_11_GSE43502"

dir.create(
  output_network,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save clean ceRNA figure
ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_publication_clean.png"),
  plot = p_ceRNA_clean,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_publication_clean.pdf"),
  plot = p_ceRNA_clean,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_network, "Figure_final_ceRNA_network_publication_clean.tiff"),
  plot = p_ceRNA_clean,
  width = 10,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

# Save workspace for ceRNA network
save.image(
  file = file.path(output_network, "Final_ceRNA_network_workspace.RData")
)

# Check saved files
list.files(
  output_network,
  full.names = TRUE
)
############################################################
# Save final ceRNA network Excel workbook
############################################################

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(openxlsx)

excel_file <- file.path(
  output_network,
  "Final_ceRNA_Network_11_downregulated_genes_publication_tables.xlsx"
)

wb <- createWorkbook()

addWorksheet(wb, "Final_11_genes")
writeData(wb, "Final_11_genes", final_11_network_genes)

addWorksheet(wb, "Cytoscape_edges")
writeData(wb, "Cytoscape_edges", final_edges)

addWorksheet(wb, "Cytoscape_nodes")
writeData(wb, "Cytoscape_nodes", final_nodes)

excluded_genes <- data.frame(
  Gene = c("CBX4", "STAT3", "THBS1"),
  Direction = "Upregulated_in_Recurrence",
  Reason = "Excluded from final ceRNA network because logFC > 0",
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Excluded_upregulated")
writeData(wb, "Excluded_upregulated", excluded_genes)

summary_ceRNA <- data.frame(
  Item = c(
    "circRNA",
    "miRNAs",
    "Final mRNAs",
    "Excluded genes",
    "Dataset",
    "Comparison",
    "Selection criteria"
  ),
  Description = c(
    "hsa_circRNA_000554",
    "hsa-miR-182-5p; hsa-miR-590-5p",
    "11 downregulated mRNAs",
    "CBX4, STAT3, THBS1",
    "GSE43502",
    "Recurrence vs No_recurrence",
    "P.Value < 0.05 and logFC < 0"
  ),
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_ceRNA)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAD3",
  border = "Bottom"
)

for (sh in names(wb)) {
  addStyle(
    wb,
    sheet = sh,
    style = header_style,
    rows = 1,
    cols = 1:50,
    gridExpand = TRUE
  )
  freezePane(wb, sheet = sh, firstRow = TRUE)
  setColWidths(wb, sheet = sh, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  excel_file,
  overwrite = TRUE
)
############################################################
# Prepare final 11 genes for PPI / STRING analysis
############################################################

output_ppi <- "13_PPI_Final_11_GSE43502"

dir.create(
  output_ppi,
  recursive = TRUE,
  showWarnings = FALSE
)

ppi_genes <- c(
  "FOXO3",
  "RECK",
  "BCL2",
  "PPM1L",
  "KDM5A",
  "CADM1",
  "CYLD",
  "HOXA9",
  "CASP2",
  "LSM14A",
  "PDCD4"
)

# Save TXT for STRING website
write.table(
  ppi_genes,
  file.path(output_ppi, "STRING_input_final_11_downregulated_genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Save CSV
write.csv(
  data.frame(Gene = ppi_genes),
  file.path(output_ppi, "STRING_input_final_11_downregulated_genes.csv"),
  row.names = FALSE
)

# Check files
list.files(
  output_ppi,
  full.names = TRUE
)
############################################################
# PPI analysis for final 11 downregulated genes
# STRING interaction files
############################################################

setwd("D:/breast cancer/breast cancer")

output_ppi <- "13_PPI_Final_11_GSE43502"

# Read STRING interaction file
ppi_edges <- read.delim(
  file.path(output_ppi, "string_interactions_short.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Read STRING node degree file
ppi_degrees <- read.delim(
  file.path(output_ppi, "string_node_degrees.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check files
head(ppi_edges)
head(ppi_degrees)

colnames(ppi_edges)
colnames(ppi_degrees)
############################################################
# Save STRING files as CSV
############################################################

write.csv(
  ppi_edges,
  file.path(output_ppi, "PPI_STRING_interactions_edges.csv"),
  row.names = FALSE
)

write.csv(
  ppi_degrees,
  file.path(output_ppi, "PPI_STRING_node_degrees.csv"),
  row.names = FALSE
)
############################################################
# Draw PPI network
############################################################

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}

if (!requireNamespace("tidygraph", quietly = TRUE)) {
  install.packages("tidygraph")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

# Build graph using node1 and node2
ppi_graph <- graph_from_data_frame(
  d = ppi_edges[, c("node1", "node2")],
  directed = FALSE
)

# Add degree
V(ppi_graph)$degree <- degree(ppi_graph)

# Plot PPI
p_ppi <- ggraph(ppi_graph, layout = "fr") +
  geom_edge_link(
    color = "grey60",
    linewidth = 0.8
  ) +
  geom_node_point(
    aes(size = degree),
    shape = 21,
    fill = "#4DBBD5",
    color = "black",
    stroke = 0.7
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 4,
    fontface = "bold"
  ) +
  scale_size_continuous(range = c(5, 12)) +
  theme_void() +
  ggtitle("PPI Network of 11 Downregulated Recurrence-Associated Genes") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

p_ppi
############################################################
# Check STRING PPI interaction column names
############################################################

colnames(ppi_edges)

head(ppi_edges)
############################################################
# Build PPI network using STRING node1 and node2 columns
############################################################

setwd("D:/breast cancer/breast cancer")

output_ppi <- "13_PPI_Final_11_GSE43502"

# Read STRING interaction file
ppi_edges <- read.delim(
  file.path(output_ppi, "string_interactions_short.tsv"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check
colnames(ppi_edges)
head(ppi_edges)

# Keep only node1, node2, and score
ppi_edge_table <- ppi_edges[
  , c("node1", "node2", "combined_score")
]

# Remove duplicate edges
ppi_edge_table <- unique(ppi_edge_table)

ppi_edge_table
############################################################
# Fix STRING PPI edge table safely
############################################################

setwd("D:/breast cancer/breast cancer")

output_ppi <- "13_PPI_Final_11_GSE43502"

# Read STRING interaction file
ppi_edges <- read.delim(
  file.path(output_ppi, "string_interactions_short.tsv"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Clean column names
colnames(ppi_edges) <- trimws(colnames(ppi_edges))

# Show column names
print(colnames(ppi_edges))

# Use first two columns as interacting genes
ppi_edge_table <- data.frame(
  node1 = ppi_edges[[1]],
  node2 = ppi_edges[[2]],
  combined_score = ppi_edges[[ncol(ppi_edges)]],
  stringsAsFactors = FALSE
)

# Remove duplicated interactions
ppi_edge_table <- unique(ppi_edge_table)

# Check final edge table
ppi_edge_table
############################################################
# Draw PPI network
############################################################

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

ppi_graph <- graph_from_data_frame(
  d = ppi_edge_table[, c("node1", "node2")],
  directed = FALSE
)

V(ppi_graph)$degree <- degree(ppi_graph)

p_ppi <- ggraph(ppi_graph, layout = "fr") +
  geom_edge_link(
    color = "grey55",
    linewidth = 1,
    alpha = 0.8
  ) +
  geom_node_point(
    aes(size = degree),
    shape = 21,
    fill = "#4DBBD5",
    color = "black",
    stroke = 0.8
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_size_continuous(range = c(5, 12)) +
  theme_void() +
  ggtitle("PPI Network of 11 Downregulated Recurrence-Associated Genes") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

p_ppi
############################################################
# Identify hub gene(s) and draw PPI network
############################################################

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

# Build graph
ppi_graph <- graph_from_data_frame(
  d = ppi_edge_table[, c("node1", "node2")],
  directed = FALSE
)

# Calculate degree
V(ppi_graph)$degree <- degree(ppi_graph)

# Define hub gene(s): highest degree
max_degree <- max(V(ppi_graph)$degree)

V(ppi_graph)$hub_status <- ifelse(
  V(ppi_graph)$degree == max_degree,
  "Hub_gene",
  "Other_genes"
)

# Check hub gene(s)
hub_genes <- V(ppi_graph)$name[V(ppi_graph)$hub_status == "Hub_gene"]
hub_genes

# Draw PPI network
p_ppi <- ggraph(ppi_graph, layout = "fr") +
  geom_edge_link(
    color = "grey60",
    linewidth = 1,
    alpha = 0.8
  ) +
  geom_node_point(
    aes(size = degree, fill = hub_status),
    shape = 21,
    color = "black",
    stroke = 0.8
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Hub_gene" = "red",
      "Other_genes" = "#4DBBD5"
    )
  ) +
  scale_size_continuous(range = c(5, 12)) +
  theme_void() +
  ggtitle("PPI Network of 11 Downregulated Recurrence-Associated Genes") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_ppi
############################################################
# Create degree table with hub status
############################################################

ppi_degree_table <- data.frame(
  Gene = V(ppi_graph)$name,
  Degree = V(ppi_graph)$degree,
  Hub_status = V(ppi_graph)$hub_status,
  stringsAsFactors = FALSE
)

ppi_degree_table <- ppi_degree_table[
  order(ppi_degree_table$Degree, decreasing = TRUE),
]

ppi_degree_table
############################################################
# Save final PPI figure with hub gene(s) in red
############################################################

ggsave(
  file.path(output_ppi, "Figure_PPI_network_hub_gene_red.png"),
  plot = p_ppi,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(output_ppi, "Figure_PPI_network_hub_gene_red.pdf"),
  plot = p_ppi,
  width = 8,
  height = 6
)

ggsave(
  file.path(output_ppi, "Figure_PPI_network_hub_gene_red.tiff"),
  plot = p_ppi,
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)
write.csv(
  ppi_degree_table,
  file.path(output_ppi, "PPI_hub_genes_ranked_by_degree.csv"),
  row.names = FALSE
)
############################################################
# Publication-quality PPI network with hub gene(s) in red
############################################################

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

# Build graph
ppi_graph <- graph_from_data_frame(
  d = ppi_edge_table[, c("node1", "node2")],
  directed = FALSE
)

# Calculate degree
V(ppi_graph)$Degree <- degree(ppi_graph)

# Define hub gene(s): highest degree
max_degree <- max(V(ppi_graph)$Degree)

V(ppi_graph)$Hub_status <- ifelse(
  V(ppi_graph)$Degree == max_degree,
  "Hub_gene",
  "Other_genes"
)

# Create node table for checking
ppi_degree_table <- data.frame(
  Gene = V(ppi_graph)$name,
  Degree = V(ppi_graph)$Degree,
  Hub_status = V(ppi_graph)$Hub_status,
  stringsAsFactors = FALSE
)

ppi_degree_table <- ppi_degree_table[
  order(ppi_degree_table$Degree, decreasing = TRUE),
]

ppi_degree_table

# Extract hub gene names
hub_genes <- ppi_degree_table$Gene[ppi_degree_table$Hub_status == "Hub_gene"]
hub_genes

# Draw publication-quality PPI network
p_ppi_pub <- ggraph(ppi_graph, layout = "fr") +
  geom_edge_link(
    color = "grey70",
    linewidth = 1,
    alpha = 0.9
  ) +
  geom_node_point(
    aes(size = Degree, fill = Hub_status),
    shape = 21,
    color = "black",
    stroke = 0.9
  ) +
  geom_node_text(
    aes(
      label = name,
      color = Hub_status,
      fontface = Hub_status
    ),
    repel = TRUE,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Hub_gene" = "#D73027",
      "Other_genes" = "#4DBBD5"
    )
  ) +
  scale_color_manual(
    values = c(
      "Hub_gene" = "#B22222",
      "Other_genes" = "black"
    )
  ) +
  scale_size_continuous(range = c(6, 14)) +
  scale_fontface_manual(
    values = c(
      "Hub_gene" = "bold",
      "Other_genes" = "plain"
    )
  ) +
  theme_void() +
  ggtitle("Protein-Protein Interaction Network of 11 Downregulated Genes") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.margin = margin(15, 15, 15, 15)
  )

p_ppi_pub
############################################################
# Safe publication-quality PPI network
############################################################

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

ppi_graph <- graph_from_data_frame(
  d = ppi_edge_table[, c("node1", "node2")],
  directed = FALSE
)

V(ppi_graph)$Degree <- degree(ppi_graph)

max_degree <- max(V(ppi_graph)$Degree)

V(ppi_graph)$Hub_status <- ifelse(
  V(ppi_graph)$Degree == max_degree,
  "Hub_gene",
  "Other_genes"
)

node_data <- data.frame(
  name = V(ppi_graph)$name,
  Degree = V(ppi_graph)$Degree,
  Hub_status = V(ppi_graph)$Hub_status,
  stringsAsFactors = FALSE
)

p_ppi_pub <- ggraph(ppi_graph, layout = "fr") +
  geom_edge_link(
    color = "grey70",
    linewidth = 1,
    alpha = 0.9
  ) +
  geom_node_point(
    aes(size = Degree, fill = Hub_status),
    shape = 21,
    color = "black",
    stroke = 0.9
  ) +
  geom_node_text(
    data = node_data[node_data$Hub_status == "Other_genes", ],
    aes(x = x, y = y, label = name),
    color = "black",
    size = 4.5,
    fontface = "plain",
    repel = TRUE
  ) +
  geom_node_text(
    data = node_data[node_data$Hub_status == "Hub_gene", ],
    aes(x = x, y = y, label = name),
    color = "#B22222",
    size = 5.2,
    fontface = "bold",
    repel = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Hub_gene" = "#D73027",
      "Other_genes" = "#4DBBD5"
    )
  ) +
  scale_size_continuous(range = c(6, 14)) +
  theme_void() +
  ggtitle("Protein-Protein Interaction Network of 11 Downregulated Genes") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.margin = margin(15, 15, 15, 15)
  )

p_ppi_pub
############################################################
# Save current PPI work before closing or changing script
############################################################

setwd("D:/breast cancer/breast cancer")

output_ppi <- "13_PPI_Final_11_GSE43502"

dir.create(
  output_ppi,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save current workspace
save.image(
  file = file.path(output_ppi, "PPI_current_workspace_backup.RData")
)
############################################################
# Save available PPI tables if they exist
############################################################

if (exists("ppi_edges")) {
  write.csv(
    ppi_edges,
    file.path(output_ppi, "PPI_STRING_raw_edges_backup.csv"),
    row.names = FALSE
  )
}

if (exists("ppi_edge_table")) {
  write.csv(
    ppi_edge_table,
    file.path(output_ppi, "PPI_edge_table_backup.csv"),
    row.names = FALSE
  )
}

if (exists("ppi_degree_table")) {
  write.csv(
    ppi_degree_table,
    file.path(output_ppi, "PPI_hub_genes_degree_backup.csv"),
    row.names = FALSE
  )
}

if (exists("ppi_graph")) {
  saveRDS(
    ppi_graph,
    file.path(output_ppi, "PPI_graph_backup.rds")
  )
}

if (exists("p_ppi")) {
  ggsave(
    file.path(output_ppi, "Figure_PPI_network_old_version_backup.png"),
    plot = p_ppi,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  ggsave(
    file.path(output_ppi, "Figure_PPI_network_old_version_backup.pdf"),
    plot = p_ppi,
    width = 8,
    height = 6
  )
}
############################################################
# Save PPI README backup
############################################################

sink(
  file.path(output_ppi, "README_PPI_backup_summary.txt")
)

cat("PPI analysis backup\n")
cat("===================\n\n")

cat("Folder: 13_PPI_Final_11_GSE43502\n")
cat("Input genes: 11 downregulated recurrence-associated target genes\n")
cat("Source: STRING database\n")
cat("Analysis: Protein-protein interaction network\n\n")

cat("Note:\n")
cat("This backup was saved before finalizing the publication-style PPI figure.\n")
cat("The old PPI plot version and available tables were saved to prevent data loss.\n\n")

cat("Files in this folder:\n")
print(list.files(output_ppi, full.names = FALSE))

sink()
load("D:/breast cancer/breast cancer/13_PPI_Final_11_GSE43502/PPI_current_workspace_backup.RData")
############################################################
# Final save before closing R
############################################################

setwd("D:/breast cancer/breast cancer")

save.image(
  file = "13_PPI_Final_11_GSE43502/PPI_final_workspace_saved_before_closing.RData"
)

list.files(
  "13_PPI_Final_11_GSE43502",
  full.names = TRUE
)
############################################################
# GO and KEGG enrichment for final 11 downregulated genes
############################################################

final_11_genes <- c(
  "FOXO3",
  "PDCD4",
  "RECK",
  "BCL2",
  "PPM1L",
  "KDM5A",
  "CADM1",
  "CYLD",
  "HOXA9",
  "CASP2",
  "LSM14A"
)

# Create output folder
output_enrich_11 <- "14_GO_KEGG_Final_11_GSE43502"

dir.create(
  output_enrich_11,
  recursive = TRUE,
  showWarnings = FALSE
)
############################################################
# Load packages
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(openxlsx)
############################################################
# Convert gene symbols to Entrez IDs
############################################################

gene_df_11 <- bitr(
  final_11_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_df_11

write.csv(
  gene_df_11,
  file.path(output_enrich_11, "Final_11_genes_EntrezIDs.csv"),
  row.names = FALSE
)
############################################################
# GO enrichment analysis
############################################################

ego_11 <- enrichGO(
  gene          = gene_df_11$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ego_11_df <- as.data.frame(ego_11)

write.csv(
  ego_11_df,
  file.path(output_enrich_11, "GO_enrichment_final_11_genes.csv"),
  row.names = FALSE
)

ego_11_df
############################################################
# KEGG enrichment analysis
############################################################

ekegg_11 <- enrichKEGG(
  gene         = gene_df_11$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg_11_df <- as.data.frame(ekegg_11)

write.csv(
  ekegg_11_df,
  file.path(output_enrich_11, "KEGG_enrichment_final_11_genes.csv"),
  row.names = FALSE
)

ekegg_11_df
############################################################
# Save GO dotplot if results exist
############################################################

if (nrow(ego_11_df) > 0) {
  
  p_go_11 <- dotplot(
    ego_11,
    showCategory = min(10, nrow(ego_11_df))
  ) +
    ggtitle("GO Enrichment of Final 11 Downregulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  p_go_11
  
  ggsave(
    file.path(output_enrich_11, "Figure_GO_dotplot_final_11_genes.png"),
    plot = p_go_11,
    width = 8,
    height = 5,
    dpi = 600
  )
  
  ggsave(
    file.path(output_enrich_11, "Figure_GO_dotplot_final_11_genes.pdf"),
    plot = p_go_11,
    width = 8,
    height = 5
  )
  
  ggsave(
    file.path(output_enrich_11, "Figure_GO_dotplot_final_11_genes.tiff"),
    plot = p_go_11,
    width = 8,
    height = 5,
    dpi = 600,
    compression = "lzw"
  )
}
############################################################
# Save KEGG dotplot if results exist
############################################################

if (nrow(ekegg_11_df) > 0) {
  
  p_kegg_11 <- dotplot(
    ekegg_11,
    showCategory = min(10, nrow(ekegg_11_df))
  ) +
    ggtitle("KEGG Enrichment of Final 11 Downregulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  p_kegg_11
  
  ggsave(
    file.path(output_enrich_11, "Figure_KEGG_dotplot_final_11_genes.png"),
    plot = p_kegg_11,
    width = 8,
    height = 5,
    dpi = 600
  )
  
  ggsave(
    file.path(output_enrich_11, "Figure_KEGG_dotplot_final_11_genes.pdf"),
    plot = p_kegg_11,
    width = 8,
    height = 5
  )
  
  ggsave(
    file.path(output_enrich_11, "Figure_KEGG_dotplot_final_11_genes.tiff"),
    plot = p_kegg_11,
    width = 8,
    height = 5,
    dpi = 600,
    compression = "lzw"
  )
}
############################################################
# Save KEGG dotplot if results exist
############################################################

if (nrow(ekegg_11_df) > 0) {
  
  p_kegg_11 <- dotplot(
    ekegg_11,
    showCategory = min(10, nrow(ekegg_11_df))
  ) +
    ggtitle("KEGG Enrichment of Final 11 Downregulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  p_kegg_11
  
  ggsave(
    file.path(output_enrich_11, "Figure_KEGG_dotplot_final_11_genes.png"),
    plot = p_kegg_11,
    width = 8,
    height = 5,
    dpi = 600
  )
  
  ggsave(
    file.path(output_enrich_11, "Figure_KEGG_dotplot_final_11_genes.pdf"),
    plot = p_kegg_11,
    width = 8,
    height = 5
  )
  
  ggsave(
    file.path(output_enrich_11, "Figure_KEGG_dotplot_final_11_genes.tiff"),
    plot = p_kegg_11,
    width = 8,
    height = 5,
    dpi = 600,
    compression = "lzw"
  )
}
############################################################
# Save GO/KEGG results in Excel workbook
############################################################

wb <- createWorkbook()

addWorksheet(wb, "Final_11_genes")
writeData(wb, "Final_11_genes", data.frame(Gene = final_11_genes))

addWorksheet(wb, "Entrez_IDs")
writeData(wb, "Entrez_IDs", gene_df_11)

addWorksheet(wb, "GO_enrichment")
writeData(wb, "GO_enrichment", ego_11_df)

addWorksheet(wb, "KEGG_enrichment")
writeData(wb, "KEGG_enrichment", ekegg_11_df)

summary_enrich_11 <- data.frame(
  Item = c(
    "Input genes",
    "GO enriched terms",
    "KEGG enriched pathways",
    "Gene set"
  ),
  Count_or_description = c(
    length(final_11_genes),
    nrow(ego_11_df),
    nrow(ekegg_11_df),
    "Final 11 downregulated recurrence-associated ceRNA genes"
  )
)

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_enrich_11)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAD3",
  border = "Bottom"
)

for (sh in names(wb)) {
  addStyle(
    wb,
    sheet = sh,
    style = header_style,
    rows = 1,
    cols = 1:50,
    gridExpand = TRUE
  )
  freezePane(wb, sheet = sh, firstRow = TRUE)
  setColWidths(wb, sheet = sh, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  file.path(output_enrich_11, "GO_KEGG_final_11_genes_results.xlsx"),
  overwrite = TRUE
)
############################################################
# Check KEGG enrichment results
############################################################

dim(ekegg_11_df)
ekegg_11_df
############################################################
# Save KEGG no significant result note
############################################################

sink(file.path(output_enrich_11, "README_KEGG_no_significant_pathways.txt"))

cat("KEGG enrichment analysis for final 11 downregulated genes\n")
cat("=========================================================\n\n")

cat("Input genes:\n")
cat(final_11_genes, sep = ", ")

cat("\n\nResult:\n")
cat("No significantly enriched KEGG pathways were detected for the final 11 genes.\n")
cat("The KEGG result table contained 0 enriched pathways.\n\n")

cat("Interpretation:\n")
cat("This may be due to the small number of input genes and the limited overlap with KEGG pathway gene sets.\n")
cat("Therefore, GO enrichment results were retained for functional interpretation.\n")

sink()
############################################################
# FINAL PROJECT BACKUP
# Save all current work and important outputs
############################################################

setwd("D:/breast cancer/breast cancer")

# Main backup folder
final_backup <- "99_FINAL_PROJECT_BACKUP_TNBC_ceRNA_Recurrence"

dir.create(
  final_backup,
  recursive = TRUE,
  showWarnings = FALSE
)

# Subfolders
dir.create(file.path(final_backup, "R_workspace"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final_backup, "Scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final_backup, "Final_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final_backup, "Final_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final_backup, "Network_files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final_backup, "PPI_files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final_backup, "GO_KEGG_files"), recursive = TRUE, showWarnings = FALSE)

# Save complete R workspace
save.image(
  file = file.path(final_backup, "R_workspace", "TNBC_ceRNA_recurrence_FINAL_workspace.RData")
)
############################################################
# Save important R objects separately if they exist
############################################################

important_objects <- c(
  "mrna_results",
  "pheno_prognosis",
  "expr_prognosis",
  "group_mrna",
  "final_nodes",
  "final_edges",
  "final_11_network_genes",
  "final_downregulated_genes",
  "ppi_edges",
  "ppi_edge_table",
  "ppi_degree_table",
  "ppi_graph",
  "ego_11",
  "ego_11_df",
  "ekegg_11",
  "ekegg_11_df",
  "gene_df_11"
)

for (obj in important_objects) {
  if (exists(obj)) {
    saveRDS(
      get(obj),
      file = file.path(final_backup, "R_workspace", paste0(obj, ".rds"))
    )
  }
}
############################################################
# Copy final output files into backup folder
############################################################

folders_to_backup <- c(
  "07_GSE43502_mRNA_results",
  "12_ceRNA_Network_Final_11_GSE43502",
  "13_PPI_Final_11_GSE43502"
)

for (folder in folders_to_backup) {
  
  if (dir.exists(folder)) {
    
    dest_folder <- file.path(final_backup, folder)
    
    dir.create(
      dest_folder,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    files_to_copy <- list.files(
      folder,
      recursive = TRUE,
      full.names = TRUE
    )
    
    for (f in files_to_copy) {
      
      rel_path <- sub(
        paste0("^", folder, "/?"),
        "",
        f
      )
      
      dest_file <- file.path(dest_folder, rel_path)
      
      dir.create(
        dirname(dest_file),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      file.copy(
        from = f,
        to = dest_file,
        overwrite = TRUE
      )
    }
  }
}
############################################################
# Save manifest of all files in the project
############################################################

project_files <- data.frame(
  File = list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE
  ),
  stringsAsFactors = FALSE
)

project_files$Size_MB <- round(file.info(project_files$File)$size / 1024^2, 3)
project_files$Modified <- file.info(project_files$File)$mtime

write.csv(
  project_files,
  file.path(final_backup, "PROJECT_FILE_MANIFEST.csv"),
  row.names = FALSE
)
############################################################
# Save final README
############################################################

sink(file.path(final_backup, "README_FINAL_PROJECT_BACKUP.txt"))

cat("TNBC ceRNA recurrence project - FINAL BACKUP\n")
cat("============================================\n\n")

cat("Main analysis:\n")
cat("- hsa_circRNA_000554-centered ceRNA network\n")
cat("- miRNAs: hsa-miR-182-5p and hsa-miR-590-5p\n")
cat("- Validation dataset: GSE43502\n")
cat("- Comparison: Recurrence vs No_recurrence\n\n")

cat("Final gene sets:\n")
cat("- 59 common miRNA target genes from TargetScan, miRDB, and miRTarBase\n")
cat("- 11 validated downregulated recurrence-associated ceRNA target genes\n")
cat("- 8 STRING-supported PPI core genes\n")
cat("- Hub gene in PPI: BCL2\n\n")

cat("Final 11 downregulated genes:\n")
if (exists("final_downregulated_genes")) {
  cat(final_downregulated_genes, sep = ", ")
} else {
  cat("FOXO3, RECK, BCL2, PPM1L, KDM5A, CADM1, CYLD, HOXA9, CASP2, LSM14A, PDCD4")
}

cat("\n\nPPI core genes:\n")
cat("BCL2, FOXO3, PDCD4, RECK, CADM1, CASP2, HOXA9, KDM5A\n\n")

cat("Excluded from final downregulated ceRNA network due to upregulation:\n")
cat("CBX4, STAT3, THBS1\n\n")

cat("Genes not included in PPI core due to lack of STRING interactions above score 0.4:\n")
cat("PPM1L, CYLD, LSM14A\n\n")

cat("GO/KEGG:\n")
cat("- GO enrichment was performed for the final 11 genes.\n")
cat("- KEGG enrichment did not identify significant pathways for the final 11 genes.\n\n")

cat("Backup contains:\n")
cat("- R workspace\n")
cat("- RDS objects\n")
cat("- ceRNA network files\n")
cat("- PPI files\n")
cat("- GO/KEGG files\n")
cat("- figures\n")
cat("- Excel/CSV tables\n")
cat("- project file manifest\n\n")

cat("Backup folder:\n")
cat(final_backup, "\n\n")

cat("Created on:\n")
cat(as.character(Sys.time()), "\n")

sink()
############################################################
# Create ZIP archive of final backup
############################################################

zip_file <- paste0(final_backup, ".zip")

if (file.exists(zip_file)) {
  file.remove(zip_file)
}

zip(
  zipfile = zip_file,
  files = final_backup
)

zip_file
############################################################
# Check final backup files
############################################################

list.files(
  final_backup,
  recursive = FALSE,
  full.names = TRUE
)

file.exists(paste0(final_backup, ".zip"))
############################################################
# Sankey plot for final 11 downregulated ceRNA network
############################################################

# Install packages if needed
if (!requireNamespace("networkD3", quietly = TRUE)) {
  install.packages("networkD3")
}

if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}

library(networkD3)
library(htmlwidgets)

############################################################
# Output folder
############################################################

output_sankey <- "12_ceRNA_Network_Final_11_GSE43502/Sankey_Final_11"

dir.create(
  output_sankey,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# Define final nodes
############################################################

nodes <- data.frame(
  name = c(
    "hsa_circRNA_000554",
    "hsa-miR-182-5p",
    "hsa-miR-590-5p",
    "FOXO3",
    "RECK",
    "BCL2",
    "PPM1L",
    "KDM5A",
    "CADM1",
    "CYLD",
    "HOXA9",
    "CASP2",
    "LSM14A",
    "PDCD4"
  ),
  group = c(
    "circRNA",
    "miRNA",
    "miRNA",
    rep("Downregulated_mRNA", 11)
  ),
  stringsAsFactors = FALSE
)

############################################################
# Define links
############################################################

links <- data.frame(
  source_name = c(
    "hsa_circRNA_000554",
    "hsa_circRNA_000554",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-182-5p",
    "hsa-miR-590-5p"
  ),
  target_name = c(
    "hsa-miR-182-5p",
    "hsa-miR-590-5p",
    "FOXO3",
    "RECK",
    "BCL2",
    "PPM1L",
    "KDM5A",
    "CADM1",
    "CYLD",
    "HOXA9",
    "CASP2",
    "LSM14A",
    "PDCD4"
  ),
  value = c(
    10,
    1,
    rep(1, 11)
  ),
  regulation = c(
    "circRNA-miRNA",
    "circRNA-miRNA",
    rep("Downregulated_mRNA", 11)
  ),
  stringsAsFactors = FALSE
)

############################################################
# Convert node names to zero-based indices for networkD3
############################################################

links$source <- match(links$source_name, nodes$name) - 1
links$target <- match(links$target_name, nodes$name) - 1

############################################################
# Color scale
############################################################

colourScale <- JS(
  'd3.scaleOrdinal()
    .domain(["circRNA", "miRNA", "Downregulated_mRNA"])
    .range(["#D9A300", "#B2182B", "#67C5D8"])'
)

############################################################
# Draw Sankey plot
############################################################

sankey_final_11 <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "group",
  LinkGroup = "regulation",
  colourScale = colourScale,
  fontSize = 13,
  nodeWidth = 28,
  nodePadding = 14,
  sinksRight = TRUE
)

sankey_final_11

############################################################
# Save as HTML
############################################################

saveWidget(
  sankey_final_11,
  file = file.path(output_sankey, "Sankey_final_11_downregulated_ceRNA_network.html"),
  selfcontained = TRUE
)
############################################################
# Save Sankey plot in publication-ready formats
############################################################

# Install required packages if needed
if (!requireNamespace("webshot2", quietly = TRUE)) {
  install.packages("webshot2")
}

if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}

library(htmlwidgets)
library(webshot2)

############################################################
# Define output folder
############################################################

output_sankey <- "12_ceRNA_Network_Final_11_GSE43502/Sankey_Final_11"

dir.create(
  output_sankey,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# Save interactive HTML
############################################################

sankey_html <- file.path(
  output_sankey,
  "Figure_4B_Sankey_final_11_downregulated_ceRNA_network.html"
)

saveWidget(
  sankey_final_11,
  file = sankey_html,
  selfcontained = TRUE
)

############################################################
# Save static PNG and PDF from HTML
############################################################

webshot(
  url = sankey_html,
  file = file.path(
    output_sankey,
    "Figure_4B_Sankey_final_11_downregulated_ceRNA_network.png"
  ),
  vwidth = 1200,
  vheight = 900,
  zoom = 2
)

webshot(
  url = sankey_html,
  file = file.path(
    output_sankey,
    "Figure_4B_Sankey_final_11_downregulated_ceRNA_network.pdf"
  ),
  vwidth = 1200,
  vheight = 900,
  zoom = 2
)
############################################################
# Save Sankey plot as TIFF
############################################################

# Install magick if needed
if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages("magick")
}

library(magick)

# Define PNG and TIFF paths
sankey_png <- file.path(
  output_sankey,
  "Figure_4B_Sankey_final_11_downregulated_ceRNA_network.png"
)

sankey_tiff <- file.path(
  output_sankey,
  "Figure_4B_Sankey_final_11_downregulated_ceRNA_network.tiff"
)

# Convert PNG to high-resolution TIFF
img <- image_read(sankey_png)

image_write(
  image = img,
  path = sankey_tiff,
  format = "tiff",
  density = "600x600"
)

# Check TIFF file
file.exists(sankey_tiff)
############################################################
# FINAL SAVE FOR TODAY'S SANKEY WORK
# hsa_circRNA_000554-centered final 11 ceRNA Sankey network
############################################################

# Set main project directory
setwd("D:/breast cancer/breast cancer")

############################################################
# Create backup folder for today's Sankey work
############################################################

sankey_backup <- "99_BACKUP_TODAY_Sankey_Final_11_ceRNA"

dir.create(
  sankey_backup,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(file.path(sankey_backup, "R_workspace"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sankey_backup, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sankey_backup, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sankey_backup, "R_objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sankey_backup, "README"), recursive = TRUE, showWarnings = FALSE)
############################################################
# Save complete R workspace
############################################################

save.image(
  file = file.path(
    sankey_backup,
    "R_workspace",
    "Workspace_today_Sankey_final_11_ceRNA.RData"
  )
)
############################################################
# Save important Sankey R objects separately
############################################################

sankey_objects <- c(
  "nodes",
  "links",
  "sankey_final_11",
  "output_sankey",
  "sankey_html",
  "sankey_png",
  "sankey_pdf",
  "sankey_tiff"
)

for (obj in sankey_objects) {
  if (exists(obj)) {
    saveRDS(
      get(obj),
      file = file.path(
        sankey_backup,
        "R_objects",
        paste0(obj, ".rds")
      )
    )
  }
}
############################################################
# Copy Sankey output files into backup folder
############################################################

# Your current Sankey output folder
output_sankey <- "12_ceRNA_Network_Final_11_GSE43502/Sankey_Final_11"

if (dir.exists(output_sankey)) {
  
  sankey_files <- list.files(
    output_sankey,
    recursive = TRUE,
    full.names = TRUE
  )
  
  for (f in sankey_files) {
    
    dest_file <- file.path(
      sankey_backup,
      "Figures",
      basename(f)
    )
    
    file.copy(
      from = f,
      to = dest_file,
      overwrite = TRUE
    )
  }
}
############################################################
# Save Sankey nodes and links tables
############################################################

if (exists("nodes")) {
  write.csv(
    nodes,
    file.path(
      sankey_backup,
      "Tables",
      "Sankey_nodes_final_11_ceRNA_network.csv"
    ),
    row.names = FALSE
  )
}

if (exists("links")) {
  write.csv(
    links,
    file.path(
      sankey_backup,
      "Tables",
      "Sankey_links_final_11_ceRNA_network.csv"
    ),
    row.names = FALSE
  )
}
############################################################
# Save README for today's Sankey backup
############################################################

sink(
  file.path(
    sankey_backup,
    "README",
    "README_TODAY_Sankey_Final_11_ceRNA.txt"
  )
)

cat("Backup for today's Sankey work\n")
cat("================================\n\n")

cat("Project:\n")
cat("TNBC recurrence-relevant ceRNA network\n\n")

cat("Sankey figure:\n")
cat("Final hsa_circRNA_000554-centered downregulated ceRNA network\n\n")

cat("Network structure:\n")
cat("hsa_circRNA_000554 -> hsa-miR-182-5p -> FOXO3, RECK, BCL2, PPM1L, KDM5A, CADM1, CYLD, HOXA9, CASP2, LSM14A\n")
cat("hsa_circRNA_000554 -> hsa-miR-590-5p -> PDCD4\n\n")

cat("Final 11 downregulated recurrence-associated mRNA targets:\n")
cat("FOXO3, PDCD4, RECK, BCL2, PPM1L, KDM5A, CADM1, CYLD, HOXA9, CASP2, LSM14A\n\n")

cat("Excluded genes from old 14-gene network:\n")
cat("CBX4, STAT3, THBS1 were excluded because they were upregulated in recurrence and were not retained in the final downregulated ceRNA network.\n\n")

cat("Important interpretation:\n")
cat("This Sankey diagram is Figure 4B or a subfigure within Figure 4.\n")
cat("It represents the final 11 downregulated ceRNA network only.\n\n")

cat("Saved items:\n")
cat("- Complete R workspace\n")
cat("- Sankey R objects\n")
cat("- HTML interactive Sankey figure\n")
cat("- PNG Sankey figure\n")
cat("- PDF Sankey figure\n")
cat("- TIFF Sankey figure\n")
cat("- Nodes CSV table\n")
cat("- Links CSV table\n\n")

cat("Suggested figure caption:\n")
cat("Figure 4B. Sankey diagram of the final hsa_circRNA_000554-centered downregulated ceRNA network. ")
cat("The Sankey diagram illustrates the final recurrence-associated ceRNA regulatory structure after exclusion of the upregulated mRNAs CBX4, STAT3, and THBS1. ")
cat("The network includes one circRNA (hsa_circRNA_000554), two validated miRNAs (hsa-miR-182-5p and hsa-miR-590-5p), and 11 downregulated recurrence-associated mRNA targets. ")
cat("hsa-miR-590-5p was linked to PDCD4, whereas hsa-miR-182-5p was linked to FOXO3, RECK, BCL2, PPM1L, KDM5A, CADM1, CYLD, HOXA9, CASP2, and LSM14A.\n\n")

cat("Backup created on:\n")
cat(as.character(Sys.time()), "\n")

sink()
############################################################
# Create ZIP backup
############################################################

zip_file <- paste0(sankey_backup, ".zip")

if (file.exists(zip_file)) {
  file.remove(zip_file)
}

zip(
  zipfile = zip_file,
  files = sankey_backup
)

file.exists(zip_file)
############################################################
# Check backup content
############################################################

list.files(
  sankey_backup,
  recursive = TRUE,
  full.names = TRUE
)

file.exists(
  file.path(
    sankey_backup,
    "R_workspace",
    "Workspace_today_Sankey_final_11_ceRNA.RData"
  )
)

file.exists(
  paste0(sankey_backup, ".zip")
)