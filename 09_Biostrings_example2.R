#!/usr/bin/env Rscript
library(Biostrings)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Drerio.UCSC.danRer11)

# In vertebrates, most protein-coding genes have a peak of GC-content near their 5′
# transcriptional start site (TSS). This feature promotes both the efficient nuclear
# export and translation of mRNAs. Despite the importance of GC-content for RNA
# metabolism, its general features, origin, and maintenance remain uncleat We will
# investigate the GC-content around the TSS of two vertebrate species: human and zebrafish.
# The goal of this exercise is to extract the sequences around the TSS and calculate mean GC-content
# in a sliding window of 50 bp. We will use the Biostrings package to read the genome sequences and
# the rtracklayer package to import the BED files containing the TSS coordinates. The analysis will
# be performed on a sample of 5000 sequences to speed up the computation. The results will be plotted
# to visualize the GC-content profile around the TSS.



# Paths to the BED files that contain coordinates for the +/- TSS regions.
hsapiens_bed_path <- "data/seq_data/tss_hs.bed"             # e.g., "data/H_sapiens_TSS.bed"
drerio_bed_path   <- "data/seq_data/tss_drerio.bed"               # e.g., "data/D_rerio_TSS.bed"


# instead of reading the genome sequences from the FASTA files, we will use the BSgenome package
# BSgenome packages are pre-packaged genome sequences for various model organisms

# -------------------------------------------------------------------
# LOAD GENOME SEQUENCES as BSgenome objects
# -------------------------------------------------------------------
# Read the FASTA files.
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Drerio.UCSC.danRer11)

# inspect data:
BSgenome.Hsapiens.UCSC.hg38
BSgenome.Drerio.UCSC.danRer11



# -------------------------------------------------------------------
# IMPORT BED FILES
# -------------------------------------------------------------------
# Read the BED files into GRanges objects.
hs_bed <- import(hsapiens_bed_path, format = "BED")
drerio_bed <- import(drerio_bed_path, format = "BED")

# inspect granges object
hs_bed
# grange and BSgenome object should have the same seqnames but they are not
# we need to adjust the seqnames in the granges object
seqlevels(hs_bed)
# append 'chr' to the seqnames
seqlevels(hs_bed) <- paste0("chr", seqlevels(hs_bed))
# same for zebrafish
seqlevels(drerio_bed)
# append 'chr' to the seqnames
seqlevels(drerio_bed) <- paste0("chr", seqlevels(drerio_bed))



# -------------------------------------------------------------------
# EXTRACT SUBSEQUENCES AT TSS REGIONS
# -------------------------------------------------------------------
# Use getSeq to extract regions from the genome based on the BED coordinates.
hs_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, hs_bed)
drerio_sequences <- getSeq(BSgenome.Drerio.UCSC.danRer11, drerio_bed)




# function letterFrequencyInSlidingView computes the frequency of specified letters in a sliding window.
# the number of sequences is quite large - use only a sample of 5000 sequences
# to speed up the computation.

n <- 5000
x <- sapply(sample(hs_sequences, n), function(seq) {
  letterFrequencyInSlidingView(seq, view.width = 50, letters = c("GC"))
})
dim(x)

# calculate mean on each position
gc_mean_content <- rowMeans(x, na.rm = TRUE)
plot(gc_mean_content)
n_windows <- length(gc_mean_content)
# we were using 50 window so TSS was in the middle
abline(v = n_windows / 2, col = "red", lty = 2)


# convert above code to a function
get_N_profile <- function(s, W=50, letters="CG", sample_size = 5000) {
  # s - DNAStringSet object
  # W - window size
  # letters - letters to count
  # sample_size - number of sequences to sample
  x <- sapply(sample(s, sample_size), function(seq) {
    letterFrequencyInSlidingView(seq, view.width = W, letters = letters)
  })
  dim(x)

  # calculate mean on each position
  mean_content <- rowMeans(x, na.rm = TRUE)
  # normalize to size of the window
  mean_content <- mean_content / W
  return(mean_content)
}

hs_gc_profile <- get_N_profile(hs_sequences, W=50, letters="GC", sample_size = 5000)
dr_gc_profile <- get_N_profile(drerio_sequences, W=50, letters="GC", sample_size = 5000)

plot(dr_gc_profile, type = "l", col = "blue", ylim = c(0,1))
points(hs_gc_profile, col = "red", type = "l", ylim = c(0,1))

# Calculate avarage nucleotide content in the genomes
gc_hs_chr1 <- letterFrequency(BSgenome.Hsapiens.UCSC.hg38$chr1, letters = "GC", as.prob = TRUE)
abline(h = gc_hs_chr1, col = "red", lty = 2)

gc_drerio_chr1 <- letterFrequency(BSgenome.Drerio.UCSC.danRer11$chr1, letters = "GC", as.prob = TRUE)
abline(h = gc_drerio_chr1, col = "blue", lty = 2)

letterFrequency(hs_sequences, "GC", as.prob = TRUE)