# GenomicRanges: Bioconductor Package for Genomic Data
#
# Introduction to GenomicRanges
# Genomic ranges represent intervals or regions in a genome:
# - Sequence name (chromosome or contig)
# - Start and end positions
# - Strand information (+ for forward, - for reverse, * for unstranded)
# - Additional metadata (e.g., gene names, scores)
#
# Core classes in GenomicRanges:
# GRanges: collection of genomic ranges with start and end locations
# GRangesList: groups of genomic ranges (e.g., exons by transcript)
# GPos: genomic positions (single points)

# Relationship with Biostrings:
# - Extract sequences from coordinates
# - Analyze sequence content within ranges
# - Connect features with annotations

# Basic GenomicRanges Operations

# Load required packages
library(GenomicRanges)
library(Biostrings)

# Create a simple GRanges object
gr <- GRanges(
  seqnames = c("chr1", "chr2", "chr2", "chr1", "chr3"),
  ranges = IRanges(start = c(10, 15, 20, 25, 30), end = c(20, 25, 30, 35, 40)),
  strand = c("+", "-", "+", "*", "-"),
  score = c(5, 10, 15, 20, 25),
  GC = c(0.45, 0.55, 0.60, 0.70, 0.80)
)

# Display and inspect the GRanges object
print(gr)
str(gr)

# Access GRanges components
seqnames(gr)
ranges(gr)
strand(gr)
mcols(gr)
gr$score
width(gr)
length(gr)

# Subsetting GRanges objects
gr[1:3]
gr[width(gr) > 10]
gr[seqnames(gr) == "chr1"]
gr[strand(gr) == "+"]
gr[gr$score > 10]
gr[seqnames(gr) == "chr2" & strand(gr) == "+"]

# Basic range operations
shift(gr, 5)            # shift by 5 bp
resize(gr, width = 10)   # resize to width 10 (anchor start)
resize(gr, width = 10, fix = "start")
resize(gr, width = 10, fix = "end")
narrow(gr, start = 2, end = -2)
flank(gr, width = 5, both = TRUE)
flank(gr, width = 5, start = TRUE)  # 5' flanking
flank(gr, width = 5, start = FALSE) # 3' flanking

# Interval operations: overlaps
# Create two GRanges for overlap demonstrations
gr1 <- GRanges("chr1", IRanges(c(5, 80, 329), width = 20))
gr2 <- GRanges("chr1", IRanges(c(10, 100, 300), width = 15))

# Find overlaps
hits <- findOverlaps(gr1, gr2)
hits
q_idx <- queryHits(hits)
s_idx <- subjectHits(hits)
gr1[q_idx]
gr2[s_idx]

# Count and subset overlaps
countOverlaps(gr1, gr2)
subsetByOverlaps(gr1, gr2)
overlapsAny(gr1, gr2)

# Nearest-neighbor operations
query   <- GRanges("chr1", IRanges(c(10, 20, 50), width = 5))
subject <- GRanges("chr1", IRanges(c(1, 15, 40, 60), width = 5))
hits_nearest <- nearest(query, subject)
hits_nearest
distance(query, subject[1:3])
hits_precede <- precede(query, subject)
hits_precede
hits_follow  <- follow(query, subject)
hits_follow

# Set operations on ranges
gr_a <- GRanges("chr1", IRanges(start = c(10, 20, 30), end = c(15, 25, 35)))
gr_b <- GRanges("chr1", IRanges(start = c(12, 22, 40), end = c(17, 27, 45)))
union(gr_a, gr_b)
intersect(gr_a, gr_b)
setdiff(gr_a, gr_b)

# Working with GRangesList
library(GenomicRanges)
grl <- GRangesList(
  sample1 = GRanges(
    seqnames = c("chr1", "chr2"),
    ranges    = IRanges(start = c(10, 20), end = c(20, 30))
  ),
  sample2 = GRanges(
    seqnames = c("chr1", "chr3"),
    ranges    = IRanges(start = c(15, 40), end = c(25, 50))
  )
)

# Inspect GRangesList
print(grl)
grlist_element <- grl[[1]]
gr$sample2
lengths(grl)
unlist(grl)
endoapply(grl, shift, 5)

# Integration with Biostrings
library(BSgenome.Hsapiens.UCSC.hg38)
genes <- GRanges(
  seqnames = c("chr1", "chr2"),
  ranges   = IRanges(start = c(10000, 20000), end = c(10500, 20500)),
  strand   = c("+", "-"),
  gene_id  = c("gene1", "gene2")
)

# Extract sequences and compute GC content
gene_seqs  <- getSeq(BSgenome.Hsapiens.UCSC.hg38, genes)
names(gene_seqs) <- genes$gene_id
print(gene_seqs)
gc_content <- letterFrequency(gene_seqs, "GC", as.prob = TRUE)
print(gc_content)
