################################################################################
# R SCRIPT FOR SESSION 9 — EXAMPLE 2:
# GC CONTENT AROUND TRANSCRIPTION START SITES (TSS)
#
# This script is a case study that combines everything from the main lesson:
# 1. Install and load BSgenome genome packages (hg38, danRer11)
# 2. Import TSS coordinates from BED files with rtracklayer
# 3. Align seqnames between BED and BSgenome (the "chr" prefix problem)
# 4. Extract sequences around each TSS with getSeq()
# 5. Compute GC content in a sliding window with letterFrequencyInSlidingView()
# 6. Turn the analysis into a reusable function
# 7. Compare the human and zebrafish GC profiles, and contrast them with the
#    genome-wide GC baseline
# 8. Polish the final figure with ggplot2
#
# BIOLOGICAL MOTIVATION:
# In vertebrates, most protein-coding genes have a peak of GC content near
# their 5' transcriptional start site (TSS). This feature is thought to
# promote efficient nuclear export and translation of the resulting mRNA.
# We will visualize this peak in human (Homo sapiens) and zebrafish
# (Danio rerio) by extracting sequences in a ±1 kb window around each TSS,
# computing mean GC content per position in a sliding window, and comparing
# the resulting profile to the genome-wide GC baseline.
#
# DATA SOURCES USED:
# - BSgenome.Hsapiens.UCSC.hg38        — human genome (installed via BiocManager)
# - BSgenome.Drerio.UCSC.danRer11      — zebrafish genome (installed via BiocManager)
# - data/seq_data/tss_hs.bed           — human TSS coordinates (±1 kb windows)
# - data/seq_data/tss_drerio.bed       — zebrafish TSS coordinates (±1 kb windows)
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################


################################################################################
# SECTION 0: INSTALLATION
################################################################################

# ---------- SLIDE: Installing the genome packages ----------

# Both BSgenome packages used below are large (~1 GB each) and are distributed
# through Bioconductor. They are typically already installed on the classroom
# workstations, but if you run this script on your own machine, install them
# once with BiocManager:
#
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
#   BiocManager::install(c(
#     "BSgenome.Hsapiens.UCSC.hg38",     # human, ~900 MB
#     "BSgenome.Drerio.UCSC.danRer11",   # zebrafish, ~450 MB
#     "rtracklayer"                      # used here for import() of BED files
#   ))
#
# Check what you have installed at any time:
#   requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
#   installed.genomes()


# ---------- SLIDE: Load libraries ----------

library(Biostrings)
library(BSgenome)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Drerio.UCSC.danRer11)
library(ggplot2)


################################################################################
# SECTION 1: BIOLOGICAL QUESTION
################################################################################

# ---------- SLIDE: The biological question ----------

# Protein-coding genes in vertebrates tend to be *embedded* in regions that are
# more GC-rich than the rest of the genome. That GC peak sits right at the
# transcription start site (TSS) and has functional consequences:
#   - it overlaps with CpG islands,
#   - it influences chromatin accessibility,
#   - it is associated with better mRNA export and translation efficiency.
#
# We will:
#   1. extract ±1 kb windows around TSSs in two species,
#   2. compute GC content in a 50-bp sliding window along each region,
#   3. average across thousands of TSSs,
#   4. compare the resulting profile to the whole-chromosome GC baseline.


################################################################################
# SECTION 2: LOADING GENOME SEQUENCES FROM BSGENOME
################################################################################

# ---------- SLIDE: Loading genomes as BSgenome objects ----------

# Unlike Example 1 (where we read FASTA files directly), we use precompiled
# BSgenome packages here. Each package is a complete, lazily-loaded reference
# genome — no FASTA parsing, no re-downloading, no manual indexing.

BSgenome.Hsapiens.UCSC.hg38
BSgenome.Drerio.UCSC.danRer11

# A BSgenome object behaves like a named collection of chromosomes.
# Chromosomes are loaded from disk only when you access them.
seqnames(BSgenome.Hsapiens.UCSC.hg38) |> head()
seqlengths(BSgenome.Hsapiens.UCSC.hg38)["chr1"]


################################################################################
# SECTION 3: IMPORTING TSS COORDINATES FROM BED
################################################################################

# ---------- SLIDE: Importing BED files with rtracklayer ----------

# The TSS coordinates were pre-computed and saved as BED files. Each row is a
# ±1 kb window centered on a transcription start site (width = 2001 bp).
# rtracklayer::import() reads a BED file into a GRanges object.

hsapiens_bed_path <- "data/seq_data/tss_hs.bed"
drerio_bed_path   <- "data/seq_data/tss_drerio.bed"

hs_bed     <- import(hsapiens_bed_path, format = "BED")
drerio_bed <- import(drerio_bed_path,   format = "BED")

# Inspect the GRanges object
hs_bed
length(hs_bed)           # how many TSSs
width(hs_bed)[1:5]       # all 2001 bp by construction

# ---------- SLIDE: Aligning seqnames — the "chr" prefix problem ----------

# A very common pitfall: the BED file uses Ensembl-style seqnames ("1", "2",
# "X", ...) while the UCSC BSgenome uses "chr1", "chr2", "chrX", ...
# getSeq() will fail silently or error if the names do not match exactly.

seqlevels(hs_bed) |> head()             # e.g. "Y", "1", "2", ...
seqlevels(BSgenome.Hsapiens.UCSC.hg38) |> head()   # "chr1", "chr2", ...

# Fix: prepend "chr" to the seqlevels of the BED-based GRanges.
seqlevels(hs_bed)     <- paste0("chr", seqlevels(hs_bed))
seqlevels(drerio_bed) <- paste0("chr", seqlevels(drerio_bed))

seqlevels(hs_bed) |> head()


################################################################################
# SECTION 4: EXTRACTING SEQUENCES AROUND THE TSS
################################################################################

# ---------- SLIDE: Extracting sequences with getSeq() ----------

# getSeq(genome, ranges) returns a DNAStringSet with one sequence per range.
# Strand is respected — ranges on the minus strand are reverse-complemented.

hs_sequences     <- getSeq(BSgenome.Hsapiens.UCSC.hg38, hs_bed)
drerio_sequences <- getSeq(BSgenome.Drerio.UCSC.danRer11, drerio_bed)

hs_sequences


################################################################################
# SECTION 5: GC CONTENT IN A SLIDING WINDOW
################################################################################

# ---------- SLIDE: letterFrequencyInSlidingView() ----------

# letterFrequencyInSlidingView(seq, view.width = W, letters = "GC")
# slides a window of width W across a single sequence and returns the count
# of the requested letters in each window.
#
# The result is a matrix with (length(seq) - W + 1) rows and one column per
# letter set. We apply this to every TSS window, stack the results, and
# average across windows to get the mean GC profile.

# The full human set has >300k TSSs — sampling 5000 is plenty for a smooth
# profile and keeps the loop fast.
n <- 5000
x <- sapply(sample(hs_sequences, n), function(seq) {
  letterFrequencyInSlidingView(seq, view.width = 50, letters = "GC")
})
dim(x)   # rows = window positions, cols = sampled TSSs

# ---------- SLIDE: Average GC content across TSSs ----------

gc_mean_content <- rowMeans(x, na.rm = TRUE)
plot(gc_mean_content,
     xlab = "Window start (bp from window 1)",
     ylab = "GC count in 50 bp window",
     main = "Mean GC content around human TSS")

# The input windows were 2001 bp centered on the TSS, and the sliding window
# is 50 bp, so the TSS itself sits near the middle of the x-axis.
n_windows <- length(gc_mean_content)
abline(v = n_windows / 2, col = "red", lty = 2)


################################################################################
# SECTION 6: TURN THE ANALYSIS INTO A FUNCTION
################################################################################

# ---------- SLIDE: A reusable profile function ----------

# The same code works for any letter set (GC, AT, just G, just CpG as "CG"),
# any window width, and any DNAStringSet. Wrapping it makes the comparison
# between species trivial.

get_N_profile <- function(s, W = 50, letters = "GC", sample_size = 5000) {
  # s            - DNAStringSet of equal-width sequences
  # W            - sliding-window width
  # letters      - letters to count (character vector or single string)
  # sample_size  - number of sequences to sample (NULL for all)

  x <- sapply(sample(s, sample_size), function(seq) {
    letterFrequencyInSlidingView(seq, view.width = W, letters = letters)
  })
  mean_content <- rowMeans(x, na.rm = TRUE)
  mean_content / W   # normalize to a per-base proportion in [0, 1]
}

hs_gc_profile <- get_N_profile(hs_sequences,     W = 50, letters = "GC", sample_size = 5000)
dr_gc_profile <- get_N_profile(drerio_sequences, W = 50, letters = "GC", sample_size = 5000)

# TASK 1:
# Use get_N_profile() to compute separate A, T, G, and C profiles for the
# human set. Plot all four on the same axes. Which bases contribute most to
# the peak at the TSS?
# Hint: call the function once per letter, then plot(..., type="l") and
# points() to overlay.


################################################################################
# SECTION 7: COMPARING SPECIES AGAINST THE GENOME-WIDE BASELINE
################################################################################

# ---------- SLIDE: Comparing species with base R ----------

plot(dr_gc_profile, type = "l", col = "blue", ylim = c(0, 1),
     xlab = "Window start (bp)", ylab = "GC content",
     main = "GC profile around TSS")
points(hs_gc_profile, type = "l", col = "red")

# ---------- SLIDE: Genome-wide GC baseline ----------

# To judge whether the TSS peak is unusual we need a reference: the
# chromosome-wide GC content. Computed once on chr1 of each species.

gc_hs_chr1 <- letterFrequency(BSgenome.Hsapiens.UCSC.hg38$chr1,
                              letters = "GC", as.prob = TRUE)
gc_dr_chr1 <- letterFrequency(BSgenome.Drerio.UCSC.danRer11$chr1,
                              letters = "GC", as.prob = TRUE)

cat("Human   chr1 GC:", round(gc_hs_chr1, 3), "\n")
cat("Zebrafish chr1 GC:", round(gc_dr_chr1, 3), "\n")

# Add the baselines to the base-R plot as horizontal lines
abline(h = gc_hs_chr1, col = "red",  lty = 2)
abline(h = gc_dr_chr1, col = "blue", lty = 2)


################################################################################
# SECTION 8: A POLISHED FIGURE WITH GGPLOT2
################################################################################

# ---------- SLIDE: The same figure with ggplot2 ----------

gc_profile_df <- data.frame(
  Position   = 1:length(hs_gc_profile),
  GC_Content = c(hs_gc_profile, dr_gc_profile),
  Species    = rep(c("Human", "Zebrafish"), each = length(hs_gc_profile))
)

ggplot(gc_profile_df, aes(x = Position, y = GC_Content, color = Species)) +
  geom_line() +
  scale_color_manual(values = c("Human" = "red", "Zebrafish" = "blue")) +
  geom_hline(yintercept = gc_hs_chr1, color = "red",  linetype = "dashed") +
  geom_hline(yintercept = gc_dr_chr1, color = "blue", linetype = "dashed") +
  labs(title = "GC content profile around TSS",
       subtitle = "Dashed lines: chromosome-1 GC baseline per species",
       x = "Position relative to TSS (bp)",
       y = "Mean GC content") +
  theme_minimal() +
  theme(legend.position = "top")

# TASK 2:
# The x-axis currently shows "window position" starting at 1. Rescale it so
# that the TSS sits at 0 and the axis runs from roughly -1000 to +1000 bp.
# Hint: the windows were 2001 bp wide, centered on the TSS, and the sliding
# view shortens the vector by (W - 1) positions.


################################################################################
# FINAL EXERCISE
################################################################################

# Repeat the analysis for CpG dinucleotides using letters = "CG" in the
# profile function, and compare the CpG profile to the GC profile around the
# human TSS. Do they peak at the same location? Hint: CpG depletion in
# vertebrate genomes is strong outside of CpG islands — expect the CpG peak
# to be much sharper than the GC peak.
