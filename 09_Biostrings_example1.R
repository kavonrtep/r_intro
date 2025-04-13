#!/usr/bin/env Rscript
library(Biostrings)
library(rtracklayer)
library(BSgenome)

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



# data - genomic data can be downloaded with get_data.sh script located in data/seq_data directory!

hsapiens_genome_path <- "data/seq_data/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"   # e.g., "data/H_sapiens.fasta"
drerio_genome_path   <- "data/seq_data/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz"      # e.g., "data/D_rerio.fasta"

# Paths to the BED files that contain coordinates for the +/- TSS regions.
hsapiens_bed_path <- "data/seq_data/tss_hs.bed"             # e.g., "data/H_sapiens_TSS.bed"
drerio_bed_path   <- "data/seq_data/tss_drerio.bed"               # e.g., "data/D_rerio_TSS.bed"



# -------------------------------------------------------------------
# LOAD GENOME SEQUENCES
# -------------------------------------------------------------------
# Read the FASTA files.
hs_genome <- readDNAStringSet(hsapiens_genome_path)
drerio_genome <- readDNAStringSet(drerio_genome_path)

# inspect data:
hs_genome
drerio_genome

# Ensure that the names are clean. (Assumes the first word in the FASTA header is the chromosome name)
names(hs_genome) <- gsub(" .*", "", names(hs_genome))
names(drerio_genome) <- gsub(" .*", "", names(drerio_genome))

# -------------------------------------------------------------------
# IMPORT BED FILES
# -------------------------------------------------------------------
# Read the BED files into GRanges objects.
hs_bed <- import(hsapiens_bed_path, format = "BED")
drerio_bed <- import(drerio_bed_path, format = "BED")

# Optional: Check that the seqnames in the BED file match the names in the genome FASTA.
# For example:
#   unique(seqnames(hs_bed))
#   names(hs_genome)

# -------------------------------------------------------------------
# EXTRACT SUBSEQUENCES AT TSS REGIONS
# -------------------------------------------------------------------
# Use getSeq to extract regions from the genome based on the BED coordinates.
hs_sequences <- getSeq(hs_genome, hs_bed)
drerio_sequences <- getSeq(drerio_genome, drerio_bed)

export(drerio_bed, drerio_bed_path, format = "BED")

# inspect extracted sequences:
hs_sequences

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

hs_content <- alphabetFrequency(hs_genome)
hs_content_total <- colSums(hs_content)/sum(hs_content)

dr_content <- alphabetFrequency(drerio_genome)
dr_content_total <- colSums(dr_content)/sum(dr_content)

## compare GC content in the genome with GC content in the TSS region - use abline to add
## horizontal line to the plot



## Task Use the function get_N_profile to calculate profiles for C,G,A,T separately and the results
