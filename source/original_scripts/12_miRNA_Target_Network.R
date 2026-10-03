library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

outdir <- "12_miRNA_Target_Network"
dir.create(outdir, showWarnings = FALSE)

targets <- data.frame(
  miRNA = c(
    rep("hsa-miR-590-5p", 9),
    rep("hsa-miR-182-5p", 51)
  ),
  Gene = c(
    "PDCD4","TMEM170A","RBPJ","STAT3","FRS2","CBX4","MBNL1","KLHL15","GID4",
    "FLOT1","MITF","FOXF2","ADCY6","RARG","NPTX1","RASA2","RECK","FBXW7",
    "PPM1L","FOXO3","EVI5","FGF9","IGF1R","BCL2","GABRB1","LSM14A","CITED2",
    "CCND2","PFN1","USP5","CREB1","MTSS1","HOXA9","NUP50","CLOCK","THBS1",
    "CASP2","DDAH1","RCC2","SATB2","LRRC4","TMEM170B","NDRG1","FOXO1",
    "EIF4EBP2","PPP1R12A","ACER2","ARRDC3","PRKAA2","CHL1","TP53INP1",
    "SESN2","STK17B","PPP1R11","KDM5A","QSER1","CYLD","NUFIP2","CADM1"
  ),
  stringsAsFactors = FALSE
)

edges <- targets %>%
  rename(from = miRNA, to = Gene)

nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  stringsAsFactors = FALSE
)

nodes$type <- ifelse(grepl("^hsa-miR", nodes$name), "miRNA", "Target gene")
nodes$degree <- degree(graph_from_data_frame(edges, vertices = nodes, directed = FALSE))

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

set.seed(123)

p <- ggraph(g, layout = "fr") +
  geom_edge_link(
    color = "grey70",
    linewidth = 0.5,
    alpha = 0.8
  ) +
  geom_node_point(
    aes(size = degree, fill = type, shape = type),
    color = "black",
    stroke = 0.8
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3.3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "miRNA" = "#E64B35",
    "Target gene" = "#4DBBD5"
  )) +
  scale_shape_manual(values = c(
    "miRNA" = 22,
    "Target gene" = 21
  )) +
  scale_size_continuous(range = c(3, 8)) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  ) +
  labs(
    title = "Consensus miRNA–target gene network",
    fill = "",
    shape = "",
    size = "Degree"
  )

print(p)

ggsave(
  file.path(outdir, "Figure_miRNA_Target_Gene_Network_FINAL.png"),
  p,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(outdir, "Figure_miRNA_Target_Gene_Network_FINAL.tiff"),
  p,
  width = 10,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Figure_miRNA_Target_Gene_Network_FINAL.pdf"),
  p,
  width = 10,
  height = 8
)

write.csv(
  targets,
  file.path(outdir, "Table_miRNA_Target_Genes_FINAL.csv"),
  row.names = FALSE
)

list.files(outdir)
#################################################
## miRNA–Target Gene Network Analysis
#################################################

# Create output folder
outdir <- "12_miRNA_Target_Network"
dir.create(outdir, showWarnings = FALSE)

# Load packages
library(igraph)
library(ggraph)
library(ggplot2)

#################################################
# Consensus target genes
#################################################

targets <- data.frame(
  miRNA = c(
    rep("hsa-miR-590-5p", 9),
    rep("hsa-miR-182-5p", 50)
  ),
  
  Gene = c(
    
    # hsa-miR-590-5p
    "PDCD4","TMEM170A","RBPJ","STAT3",
    "FRS2","CBX4","MBNL1","KLHL15","GID4",
    
    # hsa-miR-182-5p
    "FLOT1","MITF","FOXF2","ADCY6","RARG",
    "NPTX1","RASA2","RECK","FBXW7","PPM1L",
    "FOXO3","EVI5","FGF9","IGF1R","BCL2",
    "GABRB1","LSM14A","CITED2","CCND2","PFN1",
    "USP5","CREB1","MTSS1","HOXA9","NUP50",
    "CLOCK","THBS1","CASP2","DDAH1","RCC2",
    "SATB2","LRRC4","TMEM170B","NDRG1","FOXO1",
    "EIF4EBP2","PPP1R12A","ACER2","ARRDC3","PRKAA2",
    "CHL1","TP53INP1","SESN2","STK17B","PPP1R11",
    "KDM5A","QSER1","CYLD","NUFIP2","CADM1"
    
  ),
  
  stringsAsFactors = FALSE
)

#################################################
# Check counts
#################################################

table(targets$miRNA)

# Expected:
# hsa-miR-182-5p 50
# hsa-miR-590-5p 9

#################################################
# Save target table
#################################################

write.csv(
  targets,
  file.path(outdir,
            "Table_miRNA_Target_Genes_FINAL.csv"),
  row.names = FALSE
)

#################################################
# Build network
#################################################

net <- graph_from_data_frame(
  targets,
  directed = TRUE
)

#################################################
# Node information
#################################################

V(net)$type <- ifelse(
  V(net)$name %in% unique(targets$miRNA),
  "miRNA",
  "Gene"
)

V(net)$degree <- degree(net)

#################################################
# Publication-quality network
#################################################

p <- ggraph(
  net,
  layout = "fr"
) +
  
  geom_edge_link(
    colour = "grey75",
    alpha = 0.7
  ) +
  
  geom_node_point(
    aes(
      color = type,
      size = degree
    )
  ) +
  
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 3
  ) +
  
  scale_color_manual(
    values = c(
      "miRNA" = "#E64B35",
      "Gene"  = "#4DBBD5"
    )
  ) +
  
  scale_size(
    range = c(3,10)
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.title = element_blank(),
    
    legend.position = "right"
  ) +
  
  labs(
    title = "Consensus miRNA–Target Gene Network"
  )

#################################################
# Display
#################################################

print(p)

#################################################
# Save figures
#################################################

ggsave(
  file.path(
    outdir,
    "Figure_miRNA_Target_Gene_Network_FINAL.png"
  ),
  p,
  width = 12,
  height = 9,
  dpi = 600
)

ggsave(
  file.path(
    outdir,
    "Figure_miRNA_Target_Gene_Network_FINAL.tiff"
  ),
  p,
  width = 12,
  height = 9,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    outdir,
    "Figure_miRNA_Target_Gene_Network_FINAL.pdf"
  ),
  p,
  width = 12,
  height = 9
)

#################################################
# Save node information
#################################################

node_info <- data.frame(
  Node = V(net)$name,
  Type = V(net)$type,
  Degree = V(net)$degree
)

write.csv(
  node_info,
  file.path(
    outdir,
    "Table_Network_Node_Degree_FINAL.csv"
  ),
  row.names = FALSE
)

#################################################
# Finished
#################################################

cat("\nNetwork analysis completed successfully.\n")
cat("Results saved in:", outdir, "\n")