# Reconstruct a schematic from supplied final node/edge tables using base R only.
# Run from repository root: Rscript code/plot_final_network.R
# Newly prepared plotting helper; R execution was unavailable during packaging.
nodes <- read.csv('results/network/ceRNA_nodes.csv', stringsAsFactors = FALSE)
edges <- read.csv('results/network/ceRNA_edges.csv', stringsAsFactors = FALSE)
stopifnot(nrow(nodes) == 14, nrow(edges) == 13,
          !anyDuplicated(nodes$name),
          all(edges$source %in% nodes$name), all(edges$target %in% nodes$name))
colors <- c(circRNA = '#E64B35', miRNA = '#4DBBD5', mRNA = '#00A087')
shapes <- c(circRNA = 21, miRNA = 22, mRNA = 24)
pdf('figures/reconstructed_final_ceRNA_network.pdf', width = 12, height = 10, bg = 'white')
par(mar = c(2, 2, 3, 2))
plot(nodes$x, nodes$y, type = 'n', xlim = c(-5, 5), ylim = c(-7, 5),
     axes = FALSE, xlab = '', ylab = '', main = 'Final candidate ceRNA network')
for (i in seq_len(nrow(edges))) {
  a <- match(edges$source[i], nodes$name); b <- match(edges$target[i], nodes$name)
  arrows(nodes$x[a], nodes$y[a], nodes$x[b], nodes$y[b],
         length = 0.08, col = '#888888', lwd = 1.2)
}
points(nodes$x, nodes$y, pch = shapes[nodes$type], bg = colors[nodes$type],
       col = colors[nodes$type], cex = 2)
text(nodes$x, nodes$y, labels = nodes$name, pos = 3, offset = 0.7, cex = 0.65)
legend('topright', legend = names(colors), pch = shapes, pt.bg = colors,
       col = colors, bty = 'n')
dev.off()
# Arrows show candidate regulatory links, not proof of molecular causation.
