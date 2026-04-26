################################################################################
# R SCRIPT FOR SESSION 10: GENOMICRANGES
#
# This script covers:
# 1. Concept of genomic intervals and the GRanges class
# 2. A small helper function to visualize ranges as interval tracks
# 3. Creating, inspecting, and subsetting GRanges objects
# 4. Range operations: shift, resize, narrow, flank, restrict
# 5. Overlap operations: findOverlaps, countOverlaps, subsetByOverlaps
# 6. Nearest-neighbor operations: nearest, distance, precede, follow
# 7. Set operations on ranges: union, intersect, setdiff
# 8. GRangesList for hierarchical data (e.g., exons by transcript)
# 9. Integration with Biostrings and BSgenome
#
# DATA SOURCES USED:
# - Inline GRanges objects (created with GRanges() / IRanges() constructors)
# - BSgenome.Hsapiens.UCSC.hg38 (optional, large)
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################

library(GenomicRanges)
library(Biostrings)
library(Gviz)


################################################################################
# SECTION 1: WHAT ARE GENOMIC RANGES?
################################################################################

# ---------- SLIDE: What are genomic ranges? ----------

# A genomic range represents an interval on a genome and is defined by:
#   - seqname  : chromosome or contig name (e.g., "chr1")
#   - start    : 1-based start position
#   - end      : 1-based inclusive end position
#   - strand   : "+", "-", or "*" (unstranded)
#   - metadata : optional per-range columns (gene name, score, GC, ...)
#
# Genes, exons, peaks, variants, primer hits — all natural fits for ranges.

# ---------- SLIDE: Core classes in GenomicRanges ----------

# GRanges      : a vector of intervals (the workhorse)
# GRangesList  : a list of GRanges (e.g., exons grouped by transcript)
# GPos         : single-position version (rarely used in this lesson)

# ---------- SLIDE: Relationship with Biostrings ----------

# Biostrings handles SEQUENCE content (DNAStringSet, ...).
# GenomicRanges handles COORDINATES.
# Together: extract DNA at an interval, count GC, scan for motifs.


################################################################################
# SECTION 2: A HELPER FUNCTION FOR PLOTTING RANGES
################################################################################

# ---------- SLIDE: Why plot ranges? ----------

# Range operations are easier to *see* than to read off a printed table.
# We will use a small helper, plotGenomicTracks(), throughout the lesson to
# draw any named list of GRanges objects as Gviz interval tracks.
# Gviz is the Bioconductor track-plotting package — extensible, and we will
# build on it in later lessons.

# ---------- SLIDE: The plotGenomicTracks() helper ----------

# plotGenomicTracks() wraps Gviz::AnnotationTrack + plotTracks() so we can
# pass a simple named list of GRanges and get a clean multi-track figure.
#
# Arguments:
#   tracks_list : named list of GRanges objects (names = track titles)
#   colors      : optional fill color per track (default: gray)
#   chromosome  : optional chromosome to plot (default: auto)
#   from / to   : optional x-axis limits (default: auto)
#
# If a track's GRanges has an `id` metadata column, the labels are drawn
# on each interval.

plotGenomicTracks <- function(tracks_list,
                              colors     = NULL,
                              chromosome = NULL,
                              from       = NULL,
                              to         = NULL) {
  n <- length(tracks_list)
  if (is.null(names(tracks_list)) || any(names(tracks_list) == "")) {
    stop("tracks_list must be a named list (names used as track titles).")
  }
  if (!is.null(colors) && length(colors) != n) {
    stop("If colors is provided, its length must equal number of tracks.")
  }
  if (is.null(colors)) colors <- rep("gray", n)

  annotation_tracks <- mapply(function(gr, col, nm) {
    args <- list(range = gr, name = nm, fill = col,
                 col = "black", stacking = "dense")
    if ("id" %in% names(mcols(gr))) {
      ids <- mcols(gr)$id
      if (!all(is.na(ids)) && length(ids) == length(gr)) {
        args$showFeatureId     <- TRUE
        args$featureAnnotation <- "id"
        args$fontcolor.feature <- "black"
      }
    }
    do.call(AnnotationTrack, args)
  },
  gr = tracks_list, col = colors, nm = names(tracks_list),
  SIMPLIFY = FALSE)

  axis_track <- GenomeAxisTrack()
  all_tracks <- c(list(axis_track), annotation_tracks)
  plotTracks(all_tracks,
             chromosome       = chromosome,
             from             = from,
             to               = to,
             background.panel = "#FFFEDB",
             col.axis         = "black",
             cex.title        = 0.8,
             extend.left      = 5,
             extend.right     = 5)
  invisible(annotation_tracks)
}

# ---------- SLIDE: How to call plotGenomicTracks() ----------

# - Input: a *named* list of GRanges (names become track titles).
# - Output: a Gviz multi-track figure with a genome axis on top.
# - Add an `id` metadata column to label individual intervals.
#   plotGenomicTracks(list(gr1 = gr1, gr2 = gr2))

# ---------- SLIDE: Trying plotGenomicTracks() ----------

demo <- GRanges("chr1", IRanges(start = c(5, 30, 60), end = c(20, 45, 80)),
                strand = c("+", "-", "*"), id = c("a", "b", "c"))
plotGenomicTracks(list(demo = demo))


################################################################################
# SECTION 3: CREATING AND INSPECTING GRanges
################################################################################

# ---------- SLIDE: Creating a GRanges object ----------

gr <- GRanges(
  seqnames = c("chr1", "chr2", "chr2", "chr1", "chr3"),
  ranges   = IRanges(start = c(10, 15, 20, 25, 30),
                     end   = c(20, 25, 30, 35, 40)),
  strand   = c("+", "-", "+", "*", "-"),
  score    = c(5, 10, 15, 20, 25),
  GC       = c(0.45, 0.55, 0.60, 0.70, 0.80)
)
gr

# ---------- SLIDE: Visualizing our GRanges ----------

# Gviz draws one chromosome at a time — show the chr1 ranges only.
plotGenomicTracks(list(gr_chr1 = gr[seqnames(gr) == "chr1"]))

# ---------- SLIDE: Accessing GRanges components ----------

seqnames(gr)
ranges(gr)
strand(gr)
mcols(gr)
gr$score
width(gr)
length(gr)

# TASK 1:
# Create your own GRanges of three intervals on "chr1" with strands of your
# choice and a metadata column called `name`. Print it and plot it with
# plotGenomicTracks().


################################################################################
# SECTION 4: SUBSETTING GRanges
################################################################################

# ---------- SLIDE: Subsetting GRanges ----------

gr[1:3]                                         # by index
gr[width(gr) > 10]                              # by width
gr[seqnames(gr) == "chr1"]                      # by chromosome
gr[strand(gr) == "+"]                           # by strand
gr[gr$score > 10]                               # by metadata
gr[seqnames(gr) == "chr2" & strand(gr) == "+"]  # combined

# TASK 2:
# From `gr`, keep only ranges on chr1 OR chr2 with score >= 10. How many
# ranges remain?


################################################################################
# SECTION 5: BASIC RANGE OPERATIONS
################################################################################

# ---------- SLIDE: A toy GRanges to play with ----------

# A small toy GRanges on a single chromosome — easier to plot.
g <- GRanges("chr1", IRanges(start = c(10, 30, 60), end = c(20, 45, 75)),
             strand = c("+", "-", "+"))

# ---------- SLIDE: shift() ----------
# Move every range by a fixed offset.

plotGenomicTracks(list(original = g, shifted = shift(g, 5)))

# ---------- SLIDE: resize() ----------
# Set every range to a fixed width. `fix` chooses what stays put.

plotGenomicTracks(list(original = g,
                 resize_start = resize(g, width = 10, fix = "start"),
                 resize_end   = resize(g, width = 10, fix = "end"),
                 resize_center = resize(g, width = 10, fix = "center")))

# ---------- SLIDE: narrow() ----------
# Trim from either end (positive = from start, negative = from end).

plotGenomicTracks(list(original = g,
                 narrow = narrow(g, start = 2, end = -2)))

# ---------- SLIDE: flank() ----------
# Get a window adjacent to each range — useful for promoters / 3' UTRs.

plotGenomicTracks(list(original     = g,
                 flank_5prime = flank(g, width = 5, start = TRUE),
                 flank_3prime = flank(g, width = 5, start = FALSE),
                 flank_both   = flank(g, width = 5, both = TRUE)))

# ---------- SLIDE: restrict() ----------
# Clip ranges to a window; ranges outside the window are dropped.

plotGenomicTracks(list(original = g,
                 restricted = restrict(g, start = 15, end = 50)))

# TASK 3:
# Build a GRanges of 3 "genes" on chr1, then extract a 200 bp promoter for
# each (use flank() with start = TRUE and width = 200). Plot both tracks.


################################################################################
# SECTION 6: FINDING OVERLAPS
################################################################################

# ---------- SLIDE: Why overlap operations? ----------

# The single most common operation in genomics:
#   "Which of my peaks fall inside a gene?"
#   "Which variants hit an exon?"
# All use the same overlap machinery.

# ---------- SLIDE: findOverlaps() ----------

gr1 <- GRanges("chr1", IRanges(start = c(5, 80, 130), width = 20))
gr2 <- GRanges("chr1", IRanges(start = c(10, 100, 125), width = 15))
plotGenomicTracks(list(gr1 = gr1, gr2 = gr2))

hits <- findOverlaps(gr1, gr2)
hits

# ---------- SLIDE: Reading a Hits object ----------
# Each pair (queryHits[i], subjectHits[i]) is one overlapping pair.
queryHits(hits)     # indices into gr1
subjectHits(hits)   # indices into gr2

# ---------- SLIDE: countOverlaps() ----------

countOverlaps(gr1, gr2)

# ---------- SLIDE: subsetByOverlaps() ----------

gr1_overlapping <- subsetByOverlaps(gr1, gr2)
plotGenomicTracks(list(gr1 = gr1, gr2 = gr2,
                 gr1_overlapping = gr1_overlapping))

# ---------- SLIDE: overlapsAny() ----------

overlapsAny(gr1, gr2)

# TASK 4:
# Build two GRanges on chr1 of your choice and use subsetByOverlaps() to keep
# only the ranges of the first that overlap any range of the second.
# Plot all three (input, reference, result).


################################################################################
# SECTION 7: NEAREST-NEIGHBOR OPERATIONS
################################################################################

# ---------- SLIDE: nearest() and distance() ----------
# nearest(): index of closest range in subject for each query range.
# distance(): pairwise distance (0 means overlap).

query   <- GRanges("chr1", IRanges(c(10, 25, 55), width = 5))
subject <- GRanges("chr1", IRanges(c(1, 18, 40, 70), width = 5))
plotGenomicTracks(list(query = query, subject = subject))

nearest(query, subject)
distance(query, subject[1:3])

# ---------- SLIDE: precede() and follow() — what they answer ----------
#
# Unlike nearest(), which can return a subject that overlaps the query,
# precede() and follow() always look OUTSIDE the query and pick a direction:
#
#   precede(x, s) : index of the FIRST range in s that comes AFTER  x along
#                   the 5' -> 3' direction  ==> the next downstream feature.
#   follow(x, s)  : index of the LAST  range in s that came  BEFORE x along
#                   the 5' -> 3' direction  ==> the previous upstream feature.
#
# If no such range exists, the result is NA.

# ---------- SLIDE: Why "precede" feels backwards at first ----------
#
# Read the names from the SUBJECT's point of view:
#   "Which range in s does x precede?" -> the one AFTER x.
#   "Which range in s does x follow?"  -> the one BEFORE x.
#
# So precede(x, s) returns the DOWNSTREAM neighbor of x, and
#    follow(x, s)  returns the UPSTREAM   neighbor of x.

# ---------- SLIDE: A worked example ----------

q <- GRanges("chr1", IRanges(c(25), width = 5), strand = "+", id = "q")
s <- GRanges("chr1", IRanges(c(5, 15, 50, 70), width = 5),
             strand = "+", id = paste0("s", 1:4))
plotGenomicTracks(list(query = q, subject = s))

precede(q, s)   # which subject comes AFTER  q?  -> s3 (index 3)
follow(q, s)    # which subject came  BEFORE q?  -> s2 (index 2)
# So q precedes s[3] and follows s[2].

# ---------- SLIDE: Distance to the precede/follow neighbor ----------
#
# Common pattern: find upstream/downstream neighbor and how far away it is.

i_down <- precede(q, s)
i_up   <- follow(q, s)
distance(q, s[i_down])   # bp to next downstream feature
distance(q, s[i_up])     # bp to previous upstream feature
# This is exactly how you'd annotate "distance to next gene" for ChIP-seq
# peaks or variants.

# ---------- SLIDE: Why the strand matters ----------
#
# The 5' -> 3' direction depends on STRAND, but plot positions go
# left-to-right. So on a -strand gene the upstream neighbor is on the right
# of the plot, and follow() returns it.

q_minus <- GRanges("chr1", IRanges(25, width = 5), strand = "-", id = "q-")
s_plus  <- GRanges("chr1", IRanges(c(5, 50), width = 5),
                   strand = "+", id = c("s1", "s2"))
precede(q_minus, s_plus)   # NA: opposite strands

# ---------- SLIDE: Strand effects on precede/follow ----------
#
# The strand of x and subject determines whether a comparison is even
# defined. Set ignore.strand = TRUE to disable this check.
#
#   x  | subject | orientation
#   -- | ------- | -----------
#   +  |   +     | --->
#   +  |   -     | NA
#   +  |   *     | --->
#   -  |   -     | <---
#   -  |   *     | <---
#   *  |   +     | --->
#   *  |   *     | --->  (only case where * arbitrarily means +)


################################################################################
# SECTION 8: SET OPERATIONS
################################################################################

# ---------- SLIDE: Setup ----------
gr_a <- GRanges("chr1", IRanges(start = c(10, 20, 30), end = c(15, 25, 35)))
gr_b <- GRanges("chr1", IRanges(start = c(12, 22, 40), end = c(17, 27, 45)))

# ---------- SLIDE: union() ----------
plotGenomicTracks(list(gr_a = gr_a, gr_b = gr_b,
                 union = union(gr_a, gr_b)))

# ---------- SLIDE: intersect() ----------
plotGenomicTracks(list(gr_a = gr_a, gr_b = gr_b,
                 intersect = intersect(gr_a, gr_b)))

# ---------- SLIDE: setdiff() ----------
plotGenomicTracks(list(gr_a = gr_a, gr_b = gr_b,
                 setdiff = setdiff(gr_a, gr_b)))

# TASK 5:
# Using gr_a and gr_b above, compute the symmetric difference
# (parts in either but not in both). Hint: setdiff(union(...), intersect(...)).


################################################################################
# SECTION 9: GRangesList
################################################################################

# ---------- SLIDE: Why GRangesList? ----------

# Use a GRangesList when intervals come in groups: exons of a transcript,
# peaks per sample, binding sites per TF.

# ---------- SLIDE: Creating a GRangesList ----------

grl <- GRangesList(
  sample1 = GRanges(seqnames = c("chr1", "chr2"),
                    ranges   = IRanges(start = c(10, 20), end = c(20, 30))),
  sample2 = GRanges(seqnames = c("chr1", "chr3"),
                    ranges   = IRanges(start = c(15, 40), end = c(25, 50)))
)
grl

# ---------- SLIDE: Working with a GRangesList ----------

grl[[1]]              # first element
grl$sample2           # by name
elementNROWS(grl)     # length of each element
unlist(grl)           # flatten into one GRanges
endoapply(grl, shift, 5)  # apply shift() to each element


################################################################################
# SECTION 10: INTEGRATION WITH BIOSTRINGS
################################################################################

# ---------- SLIDE: Sequences from coordinates ----------
# getSeq() extracts DNA from a BSgenome at the given GRanges coordinates.

# (Optional — requires the human genome BSgenome package, ~1 GB)
# library(BSgenome.Hsapiens.UCSC.hg38)
#
# genes <- GRanges(
#   seqnames = c("chr1", "chr2"),
#   ranges   = IRanges(start = c(10000, 20000), end = c(10500, 20500)),
#   strand   = c("+", "-"),
#   gene_id  = c("gene1", "gene2")
# )
#
# gene_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, genes)
# names(gene_seqs) <- genes$gene_id
# gene_seqs
# letterFrequency(gene_seqs, "GC", as.prob = TRUE)
