################################################################################
# R SCRIPT FOR SESSION 9 — EXAMPLE 1 (FALLBACK, NO BSGENOME):
# GC CONTENT AROUND TRANSCRIPTION START SITES (TSS)
#
# This is a FALLBACK version of Example 2: if the BSgenome packages are not
# installed (and you can't install them — they are ~1.4 GB combined), you can
# still run exactly the same analysis by downloading two Ensembl FASTA files.
# The only Bioconductor packages required are Biostrings and rtracklayer,
# which are typically already installed.
#
# This script mirrors 09_Biostrings_example2.R slide-for-slide; only the
# "load genome" section is different.
#
# This script covers:
# 1. Download the genome FASTA files (one-time, ~1.4 GB total)
# 2. Read genome sequences with readDNAStringSet()
# 3. Clean FASTA header names
# 4. Import TSS coordinates from BED files with rtracklayer
# 5. Extract sequences around each TSS with getSeq()
# 6. Compute GC content in a sliding window with letterFrequencyInSlidingView()
# 7. Wrap the analysis into a reusable function
# 8. Compare the human and zebrafish profiles against the genome-wide baseline
# 9. Polish the final figure with ggplot2
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
# - data/seq_data/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz  — human genome (~1 GB)
# - data/seq_data/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz   — zebrafish genome (~500 MB)
# - data/seq_data/tss_hs.bed                                 — human TSS coordinates
# - data/seq_data/tss_drerio.bed                             — zebrafish TSS coordinates
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################


################################################################################
# SECTION 0: DOWNLOADING THE GENOME DATA
################################################################################

# ---------- SLIDE: Downloading the genome FASTA files ----------

# This script does NOT use BSgenome packages, so you don't need to install
# ~1 GB of R packages. Instead, it reads the genome straight from gzipped
# FASTA files downloaded once from Ensembl release 113:
#
#   Human (GRCh38):     ~1.0 GB
#   Zebrafish (GRCz11): ~0.5 GB
#
# The shell script data/seq_data/get_data.sh runs these two downloads:
#
#   wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
#   wget https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz
#
# From the terminal, at the repository root:
#
#   cd data/seq_data
#   bash get_data.sh
#
# Or from within R (slower, but convenient if you can't open a terminal):
#
#   dir.create("data/seq_data", recursive = TRUE, showWarnings = FALSE)
#   download.file(
#     "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz",
#     "data/seq_data/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz")
#   download.file(
#     "https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz",
#     "data/seq_data/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz")


# ---------- SLIDE: Load libraries ----------

library(Biostrings)
library(BSgenome)      # for getSeq(); no species package needed
library(rtracklayer)
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
#   4. compare the resulting profile to the whole-genome GC baseline.


################################################################################
# SECTION 2: LOADING GENOME SEQUENCES FROM FASTA
################################################################################

# ---------- SLIDE: Loading genomes from FASTA ----------

# readDNAStringSet() reads the whole (gzipped) FASTA file into a named
# DNAStringSet — one sequence per chromosome / contig / scaffold.
#
# Unlike BSgenome, nothing is loaded lazily here: the full genome is parsed
# into memory. For the human FASTA this is ~3 GB of decompressed sequence
# plus R overhead — make sure you have RAM to spare.

hsapiens_genome_path <- "data/seq_data/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
drerio_genome_path   <- "data/seq_data/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz"

hs_genome     <- readDNAStringSet(hsapiens_genome_path)
drerio_genome <- readDNAStringSet(drerio_genome_path)

hs_genome
drerio_genome

# ---------- SLIDE: Cleaning FASTA header names ----------

# Ensembl FASTA headers carry extra metadata after the chromosome name:
#   "1 dna_sm:chromosome chromosome:GRCh38:1:1:248956422:1 REF"
# getSeq() will only match names that equal the BED seqnames exactly, so we
# keep just the first word (everything before the first space).

head(names(hs_genome), 3)

names(hs_genome)     <- gsub(" .*", "", names(hs_genome))
names(drerio_genome) <- gsub(" .*", "", names(drerio_genome))

head(names(hs_genome), 5)


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

# ---------- SLIDE: Seqname conventions match here ----------

# In Example 2 (BSgenome) we had to prepend "chr" because UCSC genomes use
# "chr1", "chr2", ... while the BED file uses Ensembl-style "1", "2", ...
# With Ensembl FASTA we are already in the Ensembl world — no rename needed.

seqlevels(hs_bed) |> head()
names(hs_genome)  |> head()

# Quick sanity check — every seqlevel in the BED should appear in the FASTA
all(seqlevels(hs_bed) %in% names(hs_genome))


################################################################################
# SECTION 4: EXTRACTING SEQUENCES AROUND THE TSS
################################################################################

# ---------- SLIDE: Extracting sequences with getSeq() ----------

# getSeq() also works on a plain DNAStringSet — the set just needs named
# sequences, and the GRanges seqlevels must match those names. Strand is
# respected: minus-strand ranges are reverse-complemented automatically.

hs_sequences     <- getSeq(hs_genome,     hs_bed)
drerio_sequences <- getSeq(drerio_genome, drerio_bed)

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
# genome-wide GC content. With a plain DNAStringSet we compute it across
# ALL sequences in one go using alphabetFrequency() + colSums().

hs_content       <- alphabetFrequency(hs_genome)
hs_content_total <- colSums(hs_content) / sum(hs_content)

dr_content       <- alphabetFrequency(drerio_genome)
dr_content_total <- colSums(dr_content) / sum(dr_content)

# GC proportion = sum of the C and G columns
gc_hs_genome <- sum(hs_content_total[c("C", "G")])
gc_dr_genome <- sum(dr_content_total[c("C", "G")])

cat("Human   genome GC:", round(gc_hs_genome, 3), "\n")
cat("Zebrafish genome GC:", round(gc_dr_genome, 3), "\n")

# Add the baselines to the base-R plot as horizontal lines
abline(h = gc_hs_genome, col = "red",  lty = 2)
abline(h = gc_dr_genome, col = "blue", lty = 2)


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
  geom_hline(yintercept = gc_hs_genome, color = "red",  linetype = "dashed") +
  geom_hline(yintercept = gc_dr_genome, color = "blue", linetype = "dashed") +
  labs(title = "GC content profile around TSS",
       subtitle = "Dashed lines: whole-genome GC baseline per species",
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
