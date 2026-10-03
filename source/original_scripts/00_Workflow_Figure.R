install.packages("png")
# =========================================================
# Workflow figure for TNBC recurrence-associated ceRNA study
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

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

outdir <- "00_Workflow_Figure"
dir.create(outdir, showWarnings = FALSE)

workflow <- grViz("
digraph workflow {

  graph [
    layout = dot,
    rankdir = TB,
    bgcolor = white,
    splines = ortho,
    nodesep = 0.55,
    ranksep = 0.65
  ]

  node [
    shape = rectangle,
    style = 'rounded,filled',
    fontname = Helvetica,
    fontsize = 18,
    penwidth = 1.6,
    color = '#2F3A45',
    fillcolor = '#F7F9FB',
    margin = 0.16
  ]

  edge [
    color = '#4A5568',
    penwidth = 1.4,
    arrowsize = 0.8
  ]

  A [
    label = 'GEO dataset collection\\nGSE154255, GSE38167, GSE101123,\\nGSE182471, GSE43502',
    fillcolor = '#E8F1FA'
  ]

  B [
    label = 'Differential expression analysis\\nmiRNAs, circRNAs, and mRNAs',
    fillcolor = '#EAF7EA'
  ]

  C [
    label = 'External validation\\nhsa-miR-182-5p / hsa-miR-590-5p\\nhsa_circRNA_000554 + ROC analysis',
    fillcolor = '#FFF4E6'
  ]

  D [
    label = 'miRNA target prediction\\nTargetScan + miRDB + miRTarBase',
    fillcolor = '#F3EAF7'
  ]

  E [
    label = 'Consensus target gene selection\\nIntersection of three databases',
    fillcolor = '#F3EAF7'
  ]

  F [
    label = 'Recurrence-associated mRNA integration\\nGSE43502 recurrence-related genes',
    fillcolor = '#EAF7EA'
  ]

  G [
    label = 'Network and prognostic analyses\\nPPI network, Kaplan-Meier, Forest plot, UALCAN',
    fillcolor = '#FCE8E6'
  ]

  H [
    label = 'Final ceRNA network construction\\nhsa_circRNA_000554–miRNA–mRNA axis',
    fillcolor = '#E8F1FA'
  ]

  I [
    label = 'Outcome\\nRecurrence-associated diagnostic and prognostic\\nceRNA biomarkers in TNBC',
    fillcolor = '#DFF3F0',
    penwidth = 2.2
  ]

  A -> B -> C -> D -> E -> F -> G -> H -> I
}
")

# Display in RStudio Viewer
workflow

# Save as SVG
svg_file <- file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study.svg")
svg_txt <- export_svg(workflow)
writeLines(svg_txt, svg_file)

# Save as PNG and PDF
rsvg_png(
  svg_file,
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study.png"),
  width = 3000,
  height = 4200
)

rsvg_pdf(
  svg_file,
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study.pdf")
)

# Save as TIFF
rsvg_png(
  svg_file,
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study_temp.png"),
  width = 3000,
  height = 4200
)

png_img <- png::readPNG(
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study_temp.png")
)

tiff(
  file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study.tiff"),
  width = 3000,
  height = 4200,
  res = 300,
  compression = "lzw"
)

grid::grid.raster(png_img)
dev.off()

file.remove(file.path(outdir, "Figure1_Workflow_TNBC_ceRNA_Study_temp.png"))

list.files(outdir)