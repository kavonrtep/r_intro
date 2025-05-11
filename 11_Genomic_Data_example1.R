#############################################
## BRCA1 Gene Structure and Repeat Analysis
##
## This script demonstrates how to:
## 1. Retrieve BRCA1 gene annotation (coordinates, exons, introns)
##    using TxDb.Hsapiens.UCSC.hg38.knownGene.
## 2. Query AnnotationHub for UCSC RepeatMasker annotations on hg38.
## 3. Identify repeat elements overlapping the BRCA1 gene region,
##    and classify them by repeat class.
## 4. Determine which repeats fall in exons versus introns.
## 5. Visualize the BRCA1 gene structure alongside overlapping repeats,
##    coloring repeats by class and adding a manual legend.
#############################################

# Install and load AnnotationHub (if not already installed)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
if (!requireNamespace("AnnotationHub", quietly=TRUE))
    BiocManager::install("AnnotationHub")
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly=TRUE))
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
if (!requireNamespace("org.Hs.eg.db", quietly=TRUE))
    BiocManager::install("org.Hs.eg.db")

library(AnnotationHub)
# AnnotationHub provides programmatic access to genomic annotation resources,
# including UCSC RepeatMasker, Ensembl files, and more.

# Load required Bioconductor packages
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Gene/transcript annotations for hg38
library(org.Hs.eg.db)                       # Human gene ID mappings
library(GenomicRanges)                      # Data structures for genomic ranges
library(Gviz)                               # Genomic visualization

# Assign the TxDb object for convenience
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 1. Look up Entrez ID for BRCA1
brca_id <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "BRCA1",
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
brca_id  # should be 672 for BRCA1

# 2. Retrieve BRCA1 gene, exons, introns, and transcripts
brca_gene <- genes(txdb, filter = list(gene_id = brca_id))
brca_exons <- exons(txdb, filter = list(gene_id = brca_id))
brca_transcripts <- transcripts(txdb, filter = list(gene_id = brca_id))

brca_introns <- intronsByTranscript(txdb, use.names = TRUE)[brca_transcripts$tx_name]

# Examine exon statistics
length(brca_exons)         # total exons across all transcripts
unique(brca_exons$exon_id) # unique exon IDs

# Exon length summary
exon_lengths <- width(brca_exons)
summary(exon_lengths)      # summary stats of exon lengths
longest_exon_idx <- which.max(exon_lengths)
brca_exons[longest_exon_idx]  # inspect longest exon

# 3. Query AnnotationHub for RepeatMasker annotations
ah <- AnnotationHub()
query(ah, c("RepeatMasker", "Homo sapiens", "hg38"))
rmsk_hg38 <- ah[["AH99003"]]  # UCSC RepeatMasker annotations (GRanges)

# Subset repeats overlapping the BRCA1 gene region
brca_region <- brca_gene  # GRanges of full BRCA1 span
repeats_in_brca <- subsetByOverlaps(rmsk_hg38, brca_region)
length(repeats_in_brca)               # number of repeats in BRCA1 locus
unique(repeats_in_brca$repClass)      # repeat classes present

# 4. Identify repeats in exons vs introns
repeats_in_exons <- subsetByOverlaps(repeats_in_brca, brca_exons)
length(repeats_in_exons)              # repeats overlapping exons
table(repeats_in_exons$repClass)      # breakdown by class

# Derive intronic regions by subtracting exonic coverage from the gene span
exonic_coverage <- reduce(brca_exons)
intronic_regions <- setdiff(brca_region, exonic_coverage)
repeats_in_introns <- subsetByOverlaps(repeats_in_brca, intronic_regions)
length(repeats_in_introns)            # repeats overlapping introns

# 5. Visualization with Gviz
chrom        <- as.character(seqnames(brca_gene))
start_region <- start(brca_gene) - 5000  # pad 5kb upstream
end_region   <- end(brca_gene)   + 5000  # pad 5kb downstream

# Ideogram and axis tracks
ideo_track   <- IdeogramTrack(genome = "hg38", chromosome = chrom)
axis_track   <- GenomeAxisTrack()

# Gene model track for BRCA1
gene_track   <- GeneRegionTrack(
  txdb,
  genome = "hg38",
  chromosome = chrom,
  start = start_region,
  end = end_region,
  name = "BRCA1 Gene",
  transcriptAnnotation = "symbol"
)

# Prepare repeat classes and colors
rep_classes <- repeats_in_brca$repClass
rep_colors  <- c(
  "SINE"       = "#FF0000",
  "LINE"       = "#00FF00",
  "Retroposon" = "#0000FF",
  "LTR"        = "#FFFF00",
  "DNA"        = "#FF00FF"
)

# AnnotationTrack for repeats, grouped by class, with legend enabled
repeat_track <- AnnotationTrack(
  repeats_in_brca,
  genome       = "hg38",
  chromosome   = chrom,
  name         = "Repeats",
  group        = rep_classes,
  feature      = rep_classes,
  col          = rep_colors[rep_classes],
  fill         = rep_colors[rep_classes],
  legend       = TRUE,
  legendTitle  = "Repeat Class",
  stacking     = "dense"
)

# Plot all tracks together
plotTracks(
  list(ideo_track, axis_track, gene_track, repeat_track),
  from       = start_region,
  to         = end_region,
  chromosome = chrom,
  main       = "BRCA1 Gene Structure and Overlapping Repeats"
)

# Manually add a legend using grid graphics (base legend() won't work with grid)
library(grid)
pushViewport(viewport(x = unit(0.8, "npc"), y = unit(0.85, "npc"), width = unit(0.15, "npc"), height = unit(0.2, "npc")))
legend_labels <- names(rep_colors)
legend_colors <- rep_colors
# Draw legend entries
for (i in seq_along(legend_labels)) {
  y <- unit(1, "npc") - unit(i, "lines")
  grid.rect(x = unit(0, "npc"), y = y, width = unit(0.8, "lines"), height = unit(0.8, "lines"),
            just = c("left","center"), gp = gpar(fill = legend_colors[i], col = NA))
  grid.text(label = legend_labels[i], x = unit(1, "lines"), y = y, just = c("left","center"))
}
popViewport()

