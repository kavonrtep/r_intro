# Visualization of Genomic Data
# =============================
# (Converted from R Markdown; all original content retained as comments)

# --- YAML header ---
# title: "Visualization of Genomic Data"
# format: revealjs (self-contained)
# editor: visual

# --- Custom CSS for reveal.js slides ---
# .reveal {
#   font-size: 160%;
# }

# ------------------------------------------------------------------------------
# Packages for Genomic Data Visualization
# ------------------------------------------------------------------------------
# - GGviz
# - trackViewer
# - ggbio: Extends ggplot2 for genomic data visualization

# ------------------------------------------------------------------------------
# Install necessary Bioconductor packages
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

# ------------------------------------------------------------------------------
# Load libraries (suppress messages)
# ------------------------------------------------------------------------------
library(Gviz)            # Genomic data visualization
library(trackViewer)     # Interactive genomic track visualization
library(GenomicRanges)   # Representation of genomic intervals
library(rtracklayer)     # Import/export of genomic data formats

# ------------------------------------------------------------------------------
# Example Data: CpG Islands
# ------------------------------------------------------------------------------
data(cpgIslands)          # Load built-in CpG islands dataset
class(cpgIslands)         # Inspect object class
print(cpgIslands)         # Show contents

# ------------------------------------------------------------------------------
# Plotting with Gviz
# ------------------------------------------------------------------------------
# AnnotationTrack() creates a track object for plotting
# plotTracks() is the main interface; similar to UCSC Genome Browser output
gen <- genome(cpgIslands)     # Get genome identifier
print(gen)
chr <- seqlevels(cpgIslands)   # Get chromosome names
print(chr)

atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)             # Plot CpG islands track

# ------------------------------------------------------------------------------
# Adding genomic coordinates axis
# ------------------------------------------------------------------------------
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))         # Combine axis + annotation track

# ------------------------------------------------------------------------------
# Combining multiple tracks (custom GRanges)
# ------------------------------------------------------------------------------
# Define two sets of ranges on chr1 and chr2 with strand information
gr1 <- GRanges(
  seqnames = Rle(c("chr1", "chr2"), c(3, 2)),
  ranges   = IRanges(start = c(1, 5, 10, 15, 20),
                     end   = c(4, 8, 12, 18, 25)),
  strand   = Rle(strand(c("+", "-", "+", "-", "+")))
)
gr2 <- GRanges(
  seqnames = Rle(c("chr1", "chr2"), c(3, 2)),
  ranges   = IRanges(start = c(2, 6, 11, 16, 21),
                     end   = c(5, 9, 13, 19, 26)),
  strand   = Rle(strand(c("+", "-", "+", "-", "+")))
)

track1 <- AnnotationTrack(gr1, name = "Track 1")
track2 <- AnnotationTrack(gr2, name = "Track 2")
plotTracks(list(gtrack, track1, track2))
# Note: only chr1 region is shown by default (first chromosome)

# ------------------------------------------------------------------------------
# Specifying the region for plotting
# ------------------------------------------------------------------------------
# Use 'from', 'to', and 'chromosome' arguments to zoom in
plotTracks(
  list(gtrack, track1, track2),
  from       = 1,
  to         = 50,
  chromosome = "chr2"
)

# ------------------------------------------------------------------------------
# Adding ideogram track
# ------------------------------------------------------------------------------
# IdeogramTrack shows entire chromosome with red box indicating current region
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack))

# ------------------------------------------------------------------------------
# Adding gene models
# ------------------------------------------------------------------------------
# Load example gene models (data frame for human genome)
data(geneModels)
print(geneModels)

# Create a GeneRegionTrack for gene models
grtrack <- GeneRegionTrack(
  geneModels,
  genome     = gen,
  chromosome = chr,
  name       = "Gene Model"
)
plotTracks(list(itrack, gtrack, atrack, grtrack))

# ------------------------------------------------------------------------------
# Setting custom coordinates for gene + annotation + ideogram
# ------------------------------------------------------------------------------
plotTracks(
  list(itrack, gtrack, atrack, grtrack),
  from = 26700000,
  to   = 26750000
)

# ------------------------------------------------------------------------------
# Zooming to the base-pair level (add sequence track)
# ------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)   # Load human genome sequence
strack <- SequenceTrack(
  BSgenome.Hsapiens.UCSC.hg19,
  chromosome = chr
)
plotTracks(
  list(itrack, gtrack, atrack, grtrack, strack),
  from = 26591822,
  to   = 26591852,
  cex  = 0.8
)

# ------------------------------------------------------------------------------
# Adding a quantitative track
# ------------------------------------------------------------------------------
set.seed(123)
lim <- c(26700000, 26750000)
coords <- round(seq(lim[1], lim[2], length.out = 101))
dat    <- runif(100, min = -10, max = 10)

dtrack <- DataTrack(
  data       = dat,
  start      = coords[-length(coords)],
  end        = coords[-1],
  chromosome = chr,
  genome     = gen,
  name       = "Uniform"
)
plotTracks(
  list(itrack, gtrack, atrack, grtrack, dtrack),
  from = lim[1],
  to   = lim[2]
)

# ------------------------------------------------------------------------------
# Changing type of quantitative track
# ------------------------------------------------------------------------------
# Supported types: "histogram", "density", "heatmap", etc.
dtrack <- DataTrack(
  data       = dat,
  start      = coords[-length(coords)],
  end        = coords[-1],
  chromosome = chr,
  genome     = gen,
  name       = "Uniform",
  type       = "histogram"
)
plotTracks(
  list(itrack, gtrack, atrack, grtrack, dtrack),
  from = lim[1],
  to   = lim[2]
)

# ------------------------------------------------------------------------------
# TrackViewer: interactive genomic data visualization (alternative to Gviz)
# ------------------------------------------------------------------------------
# Plot interaction data (GInteractions / InteractionSet)
library(trackViewer)
library(InteractionSet)

# Load example chromatin interaction data
gi   <- readRDS(
  system.file("extdata", "nij.chr6.51120000.53200000.gi.rds",
              package = "trackViewer")
)
head(gi)

# Define border colors for highlighting
gi$border_color <- NA
gi$border_color[sample(seq_along(gi), 20)] <- sample(1:7, 20, replace=TRUE)

# Define TADs (topologically associated domains) as GInteractions
tads <- GInteractions(
  GRanges("chr6",
          IRanges(c(51130001, 51130001, 51450001, 52210001),
                  width = 20000)),
  GRanges("chr6",
          IRanges(c(51530001, 52170001, 52210001, 53210001),
                  width = 20000))
)
range <- GRanges("chr6", IRanges(51120000, 53200000))

# Convert gi to a track object (heatmap)
heatmap <- gi2track(gi)

# Load CTCF binding sites (ChIP-seq sample)
ctcf <- readRDS(
  system.file("extdata", "ctcf.sample.rds", package="trackViewer")
)
ctcf

# ------------------------------------------------------------------------------
# Plotting interaction data with viewTracks()
# ------------------------------------------------------------------------------
viewTracks(
  trackList(ctcf, heatmap, heightDist = c(1, 3)),
  gr             = range,
  autoOptimizeStyle = TRUE
)

# ------------------------------------------------------------------------------
# Annotation of interaction data (add TAD lines)
# ------------------------------------------------------------------------------
viewTracks(
  trackList(ctcf, heatmap, heightDist = c(1, 3)),
  gr             = range,
  autoOptimizeStyle = TRUE
)
addInteractionAnnotation(
  tads,
  "heatmap",
  grid.lines,
  gp = gpar(col = "#E69F00", lwd = 3, lty = 3)
)

# ------------------------------------------------------------------------------
# Highlighting regions with high interaction score
# ------------------------------------------------------------------------------
# Select top 5 interactions by score with distance > 200kb
gi_sub <- gi[order(gi$score, decreasing = TRUE)]
gi_sub <- head(gi_sub[distance(first(gi_sub), second(gi_sub)) > 200000], n=5)
start(regions(gi_sub)) <- start(regions(gi_sub)) - 40000
end(regions(gi_sub))   <- end(regions(gi_sub))   + 40000

# Re-plot with TADs and highlight selected regions
viewTracks(
  trackList(ctcf, heatmap, heightDist = c(1, 3)),
  gr             = range,
  autoOptimizeStyle = TRUE
)
addInteractionAnnotation(
  tads,
  "heatmap",
  grid.lines,
  gp = gpar(col = "#E69F00", lwd = 3, lty = 3)
)
addInteractionAnnotation(
  gi_sub,
  "heatmap",
  grid.polygon,
  gp = gpar(col = "red", lwd = 2, lty = 2, fill = NA)
)
