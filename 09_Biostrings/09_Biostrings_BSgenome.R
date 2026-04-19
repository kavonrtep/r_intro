################################################################################
# R SCRIPT FOR SESSION 9: BIOSTRINGS, BSGENOME, AND BIOCONDUCTOR
#
# This script covers:
# 0. Installation and Bioconductor-version compatibility (Biostrings / pwalign)
# 1. Biostrings — sequence objects (DNAString, RNAString, AAString, DNAStringSet)
# 2. Basic sequence operations (width, subseq, reverseComplement, translate)
# 3. Pattern matching (matchPattern, vmatchPattern, countPattern, mismatches)
# 4. Reading and writing FASTA files
# 5. Multiple sequence alignments and consensus
# 6. Pairwise alignment (+ BLOSUM / PAM substitution matrices)
# 7. Finding palindromes — restriction sites and hairpins
# 8. Indexing DNAStringSet objects
# 9. Sequence analysis (alphabetFrequency, oligonucleotideFrequency)
# 10. BSgenome — prepackaged reference genomes
# 11. TxDb — transcript databases
# 12. Bioconductor — the shared S4 infrastructure behind all of the above
#
# DATA SOURCES USED:
# - Inline sequences (created with DNAString / DNAStringSet constructors)
# - data/seq_data/someORF.fa         (ORF DNA sequences)
# - data/seq_data/dm3_upstream2000.fa.gz  (Drosophila upstream regions)
# - data/seq_data/UP000000625.fasta.gz    (E. coli proteome)
# - data/seq_data/example_msa.fasta       (example multiple sequence alignment)
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################

library(Biostrings)
library(BSgenome)
library(GenomicRanges)


################################################################################
# SECTION 0: INSTALLATION AND BIOCONDUCTOR-VERSION COMPATIBILITY
################################################################################

# ---------- SLIDE: Installing Biostrings, BSgenome, pwalign ----------

# All packages in this session come from Bioconductor, installed via BiocManager:
#
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
#   BiocManager::install(c(
#     "Biostrings",
#     "BSgenome",
#     "GenomicRanges",
#     "pwalign",                         # Bioc 3.20+ only (see note below)
#     "BSgenome.Hsapiens.UCSC.hg38",     # human genome — optional, large (~1 GB)
#     "BSgenome.Drerio.UCSC.danRer11"    # zebrafish genome — optional
#   ))
#
# Check your Bioconductor version with:
BiocManager::version()

# ---------- SLIDE: Heads-up — the pwalign split ----------

# In Bioconductor 3.20 (Oct 2024) and later, pairwiseAlignment() and the
# substitution-matrix data (BLOSUM*, PAM*) were moved OUT of Biostrings into
# a new package called 'pwalign'. This affects how you load them:
#
#   - Bioc ≤ 3.19  (e.g. this host, 3.18): they still live in Biostrings —
#                   no extra library() call needed.
#   - Bioc ≥ 3.20  (e.g. classroom, 3.22): install pwalign and call
#                   library(pwalign) BEFORE using pairwiseAlignment / BLOSUM62.
#
# The small helper below loads pwalign if it's installed and otherwise falls
# back to Biostrings — so the same code runs on both versions.

if (requireNamespace("pwalign", quietly = TRUE)) {
  library(pwalign)
  message("Using pwalign (Bioc >= 3.20).")
} else {
  message("pwalign not available — using pairwiseAlignment from Biostrings ",
          "(Bioc <= 3.19). If you are on Bioc >= 3.20, install pwalign:\n",
          "  BiocManager::install(\"pwalign\")")
}

# Robust loader for substitution matrices: same call works on both versions.
load_sub_matrix <- function(name) {
  if (requireNamespace("pwalign", quietly = TRUE)) {
    get(name, envir = asNamespace("pwalign"))
  } else {
    env <- new.env()
    data(list = name, package = "Biostrings", envir = env)
    env[[name]]
  }
}


################################################################################
# SECTION 1: BIOSTRINGS — CREATING SEQUENCE OBJECTS
################################################################################

# ---------- SLIDE: Biostrings: Sequence Analysis in R ----------

# Biostrings provides high-performance classes for biological sequences:
#   DNAString, RNAString, AAString   — single sequences
#   DNAStringSet, RNAStringSet, AAStringSet — collections
# Built for fast pattern matching, alignment, and frequency operations.

# ---------- SLIDE: Creating a Sequence Object ----------

# A single DNA sequence
dna_seq <- DNAString("ATGCGTACGTAGCTAG")
print(dna_seq)

# A single RNA sequence
rna_seq <- RNAString("AUGCGUACGUAGCUAG")
print(rna_seq)

# Coerce between classes with as()
rna_converted <- as(rna_seq, "RNAString")
print(rna_converted)

# A single amino acid sequence
aa_seq <- AAString("MVLSPADKTNVKAAW")
print(aa_seq)

# A collection of DNA sequences
dna_set <- DNAStringSet(c("ATGCGTACGACAGTAGCTAG",
                          "GATTACAAACATAAA",
                          "TTACATGACCCTTTACATG"))
print(dna_set)

# TASK 1:
# Create a DNAStringSet containing three named sequences of your choice
# (names like geneA, geneB, geneC). Print the object and its names().


################################################################################
# SECTION 2: BASIC SEQUENCE OPERATIONS
################################################################################

# ---------- SLIDE: Basic Sequence Operations ----------

# width()             — length of each sequence
# subseq()            — extract a subsequence
# reverseComplement() — reverse complement a DNA sequence

seq_length <- width(dna_set)
print(seq_length)

# Extract positions 4..8 from each sequence
subseq_dna <- subseq(dna_set, start = 4, end = 8)
print(subseq_dna)

# ---------- SLIDE: Reverse Complement ----------

rev_comp_dna <- reverseComplement(dna_set)
print(dna_set)
print(rev_comp_dna)

# ---------- SLIDE: Translation ----------

# translate() converts nucleotide sequences into protein sequences using the
# standard genetic code (or a custom one via the genetic.code argument).
#
# Signature:
#   translate(x, genetic.code = GENETIC_CODE,
#             no.init.codon = FALSE, if.fuzzy.codon = "error")

print(GENETIC_CODE)

translated_seq <- translate(dna_set)
print(translated_seq)

# TASK 2:
# Take dna_set and:
#   (a) compute its reverse complement,
#   (b) translate the reverse complement to protein.
# Hint: translate(reverseComplement(dna_set)).


################################################################################
# SECTION 3: PATTERN MATCHING
################################################################################

# ---------- SLIDE: Pattern Matching — Functions ----------

# matchPattern  / countPattern   — single sequence
# vmatchPattern / vcountPattern  — multiple sequences (vectorized)
#
# Key arguments:
#   max.mismatch, min.mismatch, with.indels, fixed, algorithm

# ---------- SLIDE: Pattern Matching in a Single Sequence ----------

pattern <- "CGTACG"
print(dna_seq)

match_result <- matchPattern(pattern, dna_seq)
print(match_result)

pattern_count <- countPattern(pattern, dna_seq)
print(paste("Number of occurrences of", pattern, ":", pattern_count))

# ---------- SLIDE: Pattern Matching in a Set of Sequences ----------

pattern <- "ACAT"
print(dna_set)
match_result_set <- vmatchPattern(pattern, dna_set)
print(match_result_set)

# ---------- SLIDE: Matching with Mismatches ----------

# max.mismatch allows near-matches — useful for finding approximate motifs.
pattern <- "ACAT"
dna_seq_mm <- DNAString("ACATGACATGACATGACGT")
match_result_mm <- matchPattern(pattern, dna_seq_mm, max.mismatch = 1)
print(match_result_mm)

# TASK 3:
# In dna_set, search for the motif "GTAG" allowing up to 1 mismatch.
# How many matches are found per sequence?
# Hint: vcountPattern(pattern, dna_set, max.mismatch = 1).


################################################################################
# SECTION 4: DATA IMPORT AND EXPORT
################################################################################

# ---------- SLIDE: Reading FASTA Files ----------

# readDNAStringSet() reads DNA sequences; readAAStringSet() reads proteins.
# Both handle plain text and gzipped files transparently.

fasta_file_1 <- "data/seq_data/someORF.fa"
fasta_file_2 <- "data/seq_data/dm3_upstream2000.fa.gz"

dna_seq_fasta <- readDNAStringSet(fasta_file_1)
print(dna_seq_fasta)

dna_seq_fasta_gz <- readDNAStringSet(fasta_file_2)
print(dna_seq_fasta_gz)

# ---------- SLIDE: Reading Protein FASTA ----------

protein_fasta_file <- "data/seq_data/UP000000625.fasta.gz"
protein_seq_fasta <- readAAStringSet(protein_fasta_file)
print(protein_seq_fasta)

# Writing a FASTA file:
#   writeXStringSet(dna_seq_fasta, "output.fa")
#   writeXStringSet(dna_seq_fasta, "output.fa.gz", compress = TRUE)


################################################################################
# SECTION 5: MULTIPLE SEQUENCE ALIGNMENTS
################################################################################

# ---------- SLIDE: Importing Multiple Sequence Alignments ----------

# readDNAMultipleAlignment()  — DNA alignments
# readAAMultipleAlignment()   — protein alignments
# Supported formats: "fasta", "clustal", "phylip".

aln_file <- "data/seq_data/example_msa.fasta"
aln <- readAAMultipleAlignment(aln_file, format = "fasta")
print(aln)

# ---------- SLIDE: Consensus from an MSA ----------

# consensusMatrix()  — per-position letter counts
# consensusString()  — single consensus sequence

cons_mat <- consensusMatrix(aln)
print(cons_mat)

consensus_seq <- consensusString(cons_mat)
print(consensus_seq)

# TASK 4:
# Compute the consensus matrix of 'aln' with as.prob = TRUE to get the
# per-position letter *frequencies* (rather than counts).
# Hint: consensusMatrix(aln, as.prob = TRUE).


################################################################################
# SECTION 6: PAIRWISE ALIGNMENT
################################################################################

# ---------- SLIDE: Pairwise Alignment ----------

# pairwiseAlignment() aligns two sequences. The type argument controls the
# alignment mode: "global" (default), "local", or "overlap".

seq1 <- DNAString("ATGCTAGCTAG")
seq2 <- DNAString("ATGCGGCTAG")

alignment <- pairwiseAlignment(seq1, seq2)
print(alignment)

# Extract the aligned sequences and the score
aligned_seq1 <- pattern(alignment)
aligned_seq2 <- subject(alignment)
print(aligned_seq1)
print(aligned_seq2)
print(paste("Alignment score:", score(alignment)))

# TASK 5:
# Re-run pairwiseAlignment(seq1, seq2, type = "local") and compare the score
# to the global alignment above.

# ---------- SLIDE: Scoring — substitution matrices for proteins ----------

# The default score penalises every mismatch equally. For proteins that is
# biologically wrong: a conservative substitution (I↔L, K↔R) should cost less
# than a disruptive one (W↔D). Substitution matrices such as BLOSUM62 and
# PAM250 encode the *evolutionary* probability of each substitution based on
# large databases of real homologous protein blocks.

BLOSUM62 <- load_sub_matrix("BLOSUM62")
BLOSUM45 <- load_sub_matrix("BLOSUM45")

# Peek at BLOSUM62 — rows/cols are amino acids, entries are log-odds scores
BLOSUM62[1:6, 1:6]

# Two short fragments from alpha- and beta-globin (around the heme pocket)
p1 <- AAString("PADAVMGNPKVKAHGKKVLGA")
p2 <- AAString("PKVKAHGKKVLGAFSDGLAHL")

# Same pair, different scoring models:
a62 <- pairwiseAlignment(p1, p2,
                         substitutionMatrix = BLOSUM62,
                         gapOpening = 10, gapExtension = 4)
a45 <- pairwiseAlignment(p1, p2,
                         substitutionMatrix = BLOSUM45,
                         gapOpening = 10, gapExtension = 4)

cat("BLOSUM62 score:", score(a62), "\n")
cat("BLOSUM45 score:", score(a45), "\n")

# BLOSUM45 is calibrated for more distant homologs and assigns a higher score
# to the same pair — the choice of matrix is part of the biological question
# you are asking.


################################################################################
# SECTION 7: FINDING PALINDROMES
################################################################################

# ---------- SLIDE: findPalindromes — restriction sites and hairpins ----------

# A palindrome in DNA is a region that equals its reverse complement.
# Two biological flavours:
#   - restriction enzyme recognition sites are palindromic by design
#     (EcoRI: GAATTC, BamHI: GGATCC, HindIII: AAGCTT, ...)
#   - hairpin / cruciform structures: two palindromic arms with a small loop

# Example 1 — strict palindromes (no loop): restriction-site scan
plasmid <- DNAString("AAAGAATTCCCCGAATTCAAGGATCCCCC")
res_sites <- findPalindromes(plasmid,
                             min.armlength = 3,
                             max.looplength = 0,
                             max.mismatch = 0)
print(res_sites)   # EcoRI (GAATTC) ×2 and BamHI (GGATCC)

# Example 2 — hairpins: allow a small loop between the two arms
rna_like <- DNAString("GGGAAATTTCCCAAAGGGTTTCCC")
hairs <- findPalindromes(rna_like,
                         min.armlength = 3,
                         max.looplength = 5,
                         max.mismatch = 0)
print(hairs)

# Arm extraction helpers:
#   palindromeArmLength(), palindromeLeftArm(), palindromeRightArm()
palindromeArmLength(res_sites)


################################################################################
# SECTION 8: INDEXING DNAStringSet
################################################################################

# ---------- SLIDE: Indexing XStringSet ----------

# Single brackets [] return an XStringSet (subset, preserves the collection).
# Double brackets [[]] return an XString  (a single sequence, unwrapped).

# First sequence as a DNAStringSet (length 1)
first_seq <- dna_set[1]
print(first_seq)

# First sequence as a DNAString (unwrapped)
first_seq_2 <- dna_set[[1]]
print(first_seq_2)

# Multiple sequences at once
multiple_seqs <- dna_set[1:2]
print(multiple_seqs)


################################################################################
# SECTION 9: SEQUENCE ANALYSIS
################################################################################

# ---------- SLIDE: Nucleotide Frequencies ----------

# alphabetFrequency()       — counts for every letter of the alphabet
# letterFrequency()         — counts for specified letters (incl. "GC")
# oligonucleotideFrequency()— counts of k-mers of a given width

print(dna_seq)

nucl_freq <- alphabetFrequency(dna_seq, baseOnly = TRUE)
print(nucl_freq)

# ---------- SLIDE: Top 4-mers ----------

oligo_freq <- oligonucleotideFrequency(dna_seq, width = 4)
sorted_oligo_freq <- sort(oligo_freq, decreasing = TRUE)[1:10]
print(sorted_oligo_freq)

# TASK 6:
# For dna_set, compute the GC content of each sequence as a proportion.
# Hint: letterFrequency(dna_set, "GC", as.prob = TRUE).


################################################################################
# SECTION 10: BSGENOME — PREPACKAGED REFERENCE GENOMES
################################################################################

# ---------- SLIDE: Overview of BSgenome ----------

# - Provides full reference genomes as R/Bioconductor data packages
# - Stores chromosomes as DNAString objects on disk, loaded lazily
# - Works seamlessly with Biostrings (getSeq, letterFrequency, ...) and
#   GenomicRanges (subset by GRanges intervals)

# ---------- SLIDE: Installing BSgenome Packages ----------

# Install the infrastructure:
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("BSgenome")
#
# Install specific genome packages:
#   BiocManager::install(c(
#     "BSgenome.Hsapiens.UCSC.hg38",
#     "BSgenome.Drerio.UCSC.danRer11"
#   ))

# ---------- SLIDE: Discovering Available Genomes ----------

# available.genomes()  — all BSgenome packages in your Bioconductor release
# installed.genomes()  — BSgenome packages installed on this machine
#
# head(available.genomes(), 10)
# installed.genomes()

# ---------- SLIDE: Inspecting a Genome ----------

# The examples below assume BSgenome.Hsapiens.UCSC.hg38 is installed.
# If it is not, install it first or skip this section.
#
# library(BSgenome.Hsapiens.UCSC.hg38)
# genome_human <- BSgenome.Hsapiens.UCSC.hg38
#
# seqnames(genome_human) |> head(10)     # chromosome names
# seqlengths(genome_human)["chr1"]       # chromosome length

# ---------- SLIDE: Basic Genome Queries ----------

# Extract a 10 kb window from chr1:
#   chr1_10kb <- getSeq(genome_human, "chr1", 100000, 110000)
#   print(chr1_10kb)
#
# Access a whole chromosome:
#   chr1 <- genome_human$chr1
#   print(chr1)
#
# Inspect provenance:
#   metadata(genome_human)

# ---------- SLIDE: Extracting Subsequences with getSeq() ----------

# getSeq() also works on a plain DNAStringSet with names set to chromosome-like
# identifiers — handy for quick demos without a full BSgenome package.

g <- DNAStringSet(c(chr1 = "ACGTACGTACGT",
                    chr2 = "TGCACTGCA",
                    chr3 = "GATCGATC"))
print(g)

s1 <- getSeq(g, GRanges("chr1", IRanges(1, 4)))
print(s1)

s2 <- getSeq(g, GRanges("chr2", IRanges(start = 1:4, end = 4:6)))
print(s2)

# ---------- SLIDE: Benefits of BSgenome ----------

# - Reproducibility: versioned, release-tied data packages
# - Efficiency: on-disk 2bit storage, lazy chromosome loading
# - Discoverability: available.genomes() / installed.genomes() / getBSgenome()
# - Integration: compatible with Biostrings and GenomicRanges out of the box


################################################################################
# SECTION 11: TxDb — TRANSCRIPT DATABASES
################################################################################

# ---------- SLIDE: TxDb — Transcript Databases ----------

# TxDb packages store gene/transcript/exon/CDS annotations in an on-disk
# SQLite database for a given genome build — e.g.
#   TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#   TxDb.Dmelanogaster.UCSC.dm6.ensGene
# Browse the catalog: https://bioconductor.org/packages/release/data/annotation/

# ---------- SLIDE: Installing and Loading a TxDb ----------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#
# library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# ---------- SLIDE: Core TxDb Extraction Functions ----------

# genes(txdb)                         — GRanges of all genes
# transcripts(txdb)                   — GRanges of all transcripts (with strand)
# exonsBy(txdb, by = "gene")          — GRangesList of exons per gene
# intronsByTranscript(txdb, use.names = TRUE)  — introns per transcript
# cdsBy(txdb, by = "tx")              — CDS per transcript
# promoters(txdb, upstream = 2000, downstream = 200)  — promoter windows
#
# All return GRanges / GRangesList objects that can be passed directly to
# getSeq() from BSgenome or to overlap functions in GenomicRanges.


################################################################################
# SECTION 12: BIOCONDUCTOR — THE SHARED INFRASTRUCTURE
################################################################################

# ---------- SLIDE: What is Bioconductor? ----------

# Bioconductor is an open-source project providing R packages for the analysis
# of high-throughput genomic data. Every package we touched in this session —
# Biostrings, BSgenome, GenomicRanges — is part of Bioconductor and builds on
# the same foundation: S4 classes, slots, and generic methods.

# ---------- SLIDE: Everything is an S4 object ----------

# The objects we created (dna_set, aln, g, seq1, ...) are not plain lists or
# character vectors — they are S4 objects with defined classes.

class(dna_set)     # DNAStringSet
class(aln)         # AAMultipleAlignment
class(seq1)        # DNAString

# isS4() confirms they use the S4 object system:
isS4(dna_set)
isS4(aln)

# ---------- SLIDE: Slots — the internal storage ----------

# S4 objects store their data in named slots. Peek at them with slotNames().
# You rarely touch slots directly — use accessor functions instead.

slotNames(dna_set)      # internal fields of a DNAStringSet
slotNames(aln)          # internal fields of an AAMultipleAlignment

# ---------- SLIDE: Accessors, not $ ----------

# For S4 objects, `$` does NOT work the same way as for data frames.
# Use the accessor functions provided by the package.

# Correct — use accessors:
width(dna_set)
names(g)
alphabet(dna_seq)

# This would error or return NULL for an S4 object (uncomment to try):
# dna_set$width

# Why? Accessors enforce the class contract: they return what the class
# promises, even if the internal slot representation changes.

# ---------- SLIDE: Methods dispatch on class ----------

# A generic function (e.g., print, width, length) runs a *different* method
# depending on the class of its argument.

methods(class = "DNAStringSet")           # functions defined for DNAStringSet
methods(class = "MultipleAlignment")      # and for MultipleAlignment

# The same function name — different behaviour per class:
length(dna_set)               # number of sequences in the set
length(dna_seq)               # number of characters in the one sequence

# And `width` is defined for the set but not for a single DNAString:
width(dna_set)                # length of each sequence in the set

# ---------- SLIDE: Generic functions across packages ----------

# The real power: one verb works across many classes, even across packages.
# getSeq() is defined for plain DNAStringSet objects, BSgenome objects, and
# FaFile (indexed FASTA) — you call it the same way every time:

# On a DNAStringSet with named "chromosomes":
getSeq(g, GRanges("chr1", IRanges(1, 4)))

# On a BSgenome (if installed):
#   getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr1", 100000, 110000)

# Same verb, same arguments — different class, different implementation.
# This is what lets Bioconductor packages compose cleanly.

# ---------- SLIDE: Finding help and vignettes ----------

# ?DNAStringSet            — help page for a class
# methods("translate")     — every class translate() dispatches on
# showMethods("getSeq")    — signatures for an S4 generic
# vignette(package = "Biostrings")  — list vignettes for a package
# browseVignettes("Biostrings")     — open them in a browser
#
# Vignettes are worked-example tutorials — usually the fastest way to learn a
# new Bioconductor package after this course.


################################################################################
# FINAL EXERCISE
################################################################################

# Using the FASTA file data/seq_data/dm3_upstream2000.fa.gz (Drosophila
# promoter regions) and the tools from this session:
#
#   1. Load the sequences with readDNAStringSet().
#   2. Report how many sequences there are and their mean width.
#   3. Compute the GC content of each sequence (letterFrequency(..., "GC",
#      as.prob = TRUE)) and plot the distribution with hist().
#   4. Pick a short motif (e.g. "TATAAA") and count matches per sequence with
#      vcountPattern(); also allow 1 mismatch and compare the totals.
