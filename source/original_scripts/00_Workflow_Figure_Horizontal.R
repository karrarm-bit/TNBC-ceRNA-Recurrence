# =========================================================
# Horizontal workflow figure with selection criteria
# TNBC recurrence-associated ceRNA study
# =========================================================

if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
  install.packages("DiagrammeR")
}
if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
  install.packages("DiagrammeRsvg")
}
if (!requireNamespace("rsvg", quietly = TRUE)) {
  install.packages("rsvg")
}
if (!requireNamespace("png", quietly = TRUE)) {
  install.packages("png")
}

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(png)
library(grid)

outdir <- "00_Workflow_Figure"
dir.create(outdir, showWarnings = FALSE)

workflow <- grViz("
digraph workflow {

  graph [
    layout = dot,
    rankdir = LR,
    bgcolor = white,
    splines = ortho,
    nodesep = 0.55,
    ranksep = 0.55
  ]

  node [
    shape = rectangle,
    style = 'rounded,filled',
    fontname = Helvetica,
    fontsize = 16,
    penwidth = 1.8,
    color = '#2F3A45',
    margin = 0.15
  ]

  edge [
    color = '#4A5568',
    penwidth = 1.5,
    arrowsize = 0.8
  ]

  A [
    label = 'GEO dataset collection\\nGSE154255 | GSE38167\\nGSE101123 | GSE182471\\nGSE43502',
    fillcolor = '#E8F1FA'
  ]

  B [
    label = 'Differential expression analysis\\nmiRNA + circRNA + mRNA\\nCriteria: P < 0.05, |log2FC| > 1',
    fillcolor = '#EAF7EA'
  ]

  C [
    label = 'External validation\\nmiRNA validation: GSE38167\\ncircRNA validation: GSE182471\\nCriteria: independent dataset + matched annotation',
    fillcolor = '#FFF4E6'
  ]

  D [
    label = 'miRNA target prediction\\nTargetScan ∩ miRDB ∩ miRTarBase\\nCriteria: consensus targets',
    fillcolor = '#F3EAF7'
  ]

  E [
    label = 'Recurrence mRNA integration\\nOverlap with GSE43502\\nCriteria: recurrence-associated mRNAs',
    fillcolor = '#EAF7EA'
  ]

  F [
    label = 'Network and prognostic validation\\nPPI network | Kaplan–Meier\\nForest plot | UALCAN | ROC',
    fillcolor = '#FCE8E6'
  ]

  G [
    label = 'Final ceRNA regulatory axis\\nhsa_circRNA_000554\\nhsa-miR-182-5p / hsa-miR-590-5p\\nrecurrence-associated mRNAs',
    fillcolor = '#DFF3F0',
    penwidth = 2.4
  ]

  A -> B -> C -> D -> E -> F -> G
}
")

workflow

# Save SVG
svg_file <- file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Horizontal.svg")
svg_txt <- export_svg(workflow)
writeLines(svg_txt, svg_file)

# Save PNG
png_file <- file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Horizontal.png")
rsvg_png(
  svg_file,
  png_file,
  width = 5200,
  height = 1600
)

# Save PDF
rsvg_pdf(
  svg_file,
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Horizontal.pdf")
)

# Save TIFF
temp_png <- file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Horizontal_temp.png")
rsvg_png(
  svg_file,
  temp_png,
  width = 5200,
  height = 1600
)

img <- png::readPNG(temp_png)

tiff(
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Horizontal.tiff"),
  width = 5200,
  height = 1600,
  res = 300,
  compression = "lzw"
)

grid::grid.raster(img)
dev.off()

file.remove(temp_png)

list.files(outdir)