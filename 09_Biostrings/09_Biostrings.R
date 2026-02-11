#!/usr/bin/env Rscript
# =============================================================================
# Title: Biostrings: Bioconductor Package for Sequence Analysis
# Description: This script demonstrates how to use the Biostrings package along
# with related Bioconductor packages to perform a variety of sequence operations.
# Operations include creating sequence objects, basic manipulations, translation,
# pattern matching, data import/export, multiple sequence alignments, pairwise
# alignment, accessing sequences, and basic sequence analysis.
# =============================================================================

# Load required libraries
library(Biostrings)       # For sequence handling and analysis
library(BSgenome)         # For genome sequence handling and getSeq()
library(GenomicRanges)    # For GRanges and IRanges objects used in getSeq()

# -----------------------------------------------------------------------------
# Section: Creating Sequence Objects
# -----------------------------------------------------------------------------
# Create and print various sequence objects (DNA, RNA, and protein).

# Create a DNA sequence
dna_seq <- DNAString("ATGCGTACGTAGCTAG")
print(dna_seq)

# Create an RNA sequence
rna_seq <- RNAString("AUGCGUACGUAGCUAG")
print(rna_seq)

# Convert the RNA sequence to an RNAString object (illustrative conversion)
rna_converted <- as(rna_seq, "RNAString")
print(rna_converted)

# Create an amino acid (protein) sequence
aa_seq <- AAString("MVLSPADKTNVKAAW")
print(aa_seq)

# Create a set of DNA sequences (DNAStringSet)
dna_set <- DNAStringSet(c("ATGCGTACGACAGTAGCTAG",
                          "GATTACAAACATAAA",
                          "TTACATGACCCTTTACATG"))
print(dna_set)

# -----------------------------------------------------------------------------
# Section: Basic Sequence Operations
# Description: Compute sequence lengths, extract subsequences and obtain the
# reverse complement of sequences.
# -----------------------------------------------------------------------------

# Get the length (width) of each sequence in the DNAStringSet
seq_length <- width(dna_set)
print(seq_length)

# Extract a subsequence (positions 4 to 8) from each sequence in the DNAStringSet
subseq_dna <- subseq(dna_set, start = 4, end = 8)
print(subseq_dna)

# Get the reverse complement of each DNA sequence in the set
rev_comp_dna <- reverseComplement(dna_set)
print("Original DNA sequences:")
print(dna_set)
print("Reverse complement sequences:")
print(rev_comp_dna)

# -----------------------------------------------------------------------------
# Section: Translation
# Description: Translate nucleotide sequences into protein sequences.
# The translate() function converts a nucleotide sequence (or set) to its
# corresponding amino acid sequence using the genetic code.
#
# The following commented code illustrates the translate() function parameters:
# translate(x, genetic.code = GENETIC_CODE, no.init.codon = FALSE,
#           if.fuzzy.codon = "error")
# -----------------------------------------------------------------------------

# Print the standard genetic code table
print(GENETIC_CODE)

# Translate the DNA sequences in the DNAStringSet to protein sequences
translated_seq <- translate(dna_set)
print("Translated protein sequences:")
print(translated_seq)

# -----------------------------------------------------------------------------
# Section: Pattern Matching
# Description: Search for specific motifs (patterns) in sequences.
# -----------------------------------------------------------------------------

# Define a short pattern (motif) to search within the DNA sequence
pattern <- "CGTACG"
print("Searching the following DNA sequence:")
print(dna_seq)

# Use matchPattern() to find the motif in the DNA sequence
match_result <- matchPattern(pattern, dna_seq)
print("Match result for pattern in single DNA sequence:")
print(match_result)

# Count the number of occurrences of the pattern in the DNA sequence
pattern_count <- countPattern(pattern, dna_seq)
print(paste("Number of occurrences of", pattern, ":", pattern_count))

# Now, search for a different motif in a set of DNA sequences using vmatchPattern
pattern <- "ACAT"
print("DNA sequence set:")
print(dna_set)
match_result_set <- vmatchPattern(pattern, dna_set)
print("Match result in DNA sequence set:")
print(match_result_set)

# Matching with mismatches:
# Define a DNA sequence with repeated segments and allow for up to 1 mismatch.
dna_seq_mismatch <- DNAString("ACATGACATGACATGACGT")
match_result_mismatch <- matchPattern(pattern, dna_seq_mismatch, max.mismatch = 1)
print("Match result with up to 1 mismatch:")
print(match_result_mismatch)

# -----------------------------------------------------------------------------
# Section: Data Import and Export
# Description: Read DNA and protein sequence data from FASTA files. This
# includes handling both plain text and gzipped files.
# -----------------------------------------------------------------------------

# Define file paths for FASTA files (adjust the paths as needed)
fasta_file_1 <- "data/seq_data/someORF.fa"
fasta_file_2 <- "data/seq_data/dm3_upstream2000.fa.gz"

# Read a FASTA file containing DNA sequences
dna_seq_fasta <- readDNAStringSet(fasta_file_1)
print("DNA sequences from FASTA file:")
print(dna_seq_fasta)

# Read a gzipped FASTA file containing DNA sequences
dna_seq_fasta_gz <- readDNAStringSet(fasta_file_2)
print("DNA sequences from gzipped FASTA file:")
print(dna_seq_fasta_gz)

# Read a protein sequence from a FASTA file
protein_fasta_file <- "data/seq_data/UP000000625.fasta.gz"
protein_seq_fasta <- readAAStringSet(protein_fasta_file)
print("Protein sequences from FASTA file:")
print(protein_seq_fasta)

# -----------------------------------------------------------------------------
# Section: Importing Multiple Sequence Alignments (MSA)
# Description: Load a multiple sequence alignment file in FASTA format.
# -----------------------------------------------------------------------------

# Define the file path for the MSA file
aln_file <- "data/seq_data/example_msa.fasta"
# Read the multiple sequence alignment; here we assume an amino acid alignment.
aln <- readAAMultipleAlignment(aln_file, format = "fasta")
print("Multiple sequence alignment:")
print(aln)

# -----------------------------------------------------------------------------
# Section: Analyzing Multiple Sequence Alignments
# Description: Compute the consensus matrix and generate a consensus sequence
# from an MSA.
# -----------------------------------------------------------------------------

# Compute the consensus matrix from the alignment
cons_mat <- consensusMatrix(aln)
print("Consensus matrix:")
print(cons_mat)

# Generate the consensus sequence from the consensus matrix
consensus_seq <- consensusString(cons_mat)
print("Consensus sequence from MSA:")
print(consensus_seq)

# -----------------------------------------------------------------------------
# Section: Pairwise Alignment
# Description: Compare two sequences using pairwise alignment.
# -----------------------------------------------------------------------------

# Create two DNA sequences for pairwise alignment
seq1 <- DNAString("ATGCTAGCTAG")
seq2 <- DNAString("ATGCGGCTAG")

# Perform pairwise alignment with default parameters
alignment <- pairwiseAlignment(seq1, seq2)
print("Pairwise alignment result:")
print(alignment)

# Retrieve and print aligned sequences and the alignment score
aligned_seq1 <- pattern(alignment)
aligned_seq2 <- subject(alignment)
print("Aligned sequence 1:")
print(aligned_seq1)
print("Aligned sequence 2:")
print(aligned_seq2)
print(paste("Alignment score:", score(alignment)))

# -----------------------------------------------------------------------------
# Section: Accessing Individual Sequences from DNAStringSet
# Description: Extract specific sequences from a DNAStringSet object using
# single ([) and double ([[) bracket indexing.
# -----------------------------------------------------------------------------

# Access the first sequence in the DNAStringSet (returns an XStringSet)
first_seq <- dna_set[1]
print("First sequence (using single brackets):")
print(first_seq)

# Access the first sequence using double brackets (returns an XString)
first_seq_2 <- dna_set[[1]]
print("First sequence (using double brackets):")
print(first_seq_2)

# Access multiple sequences (e.g., first two sequences) from the DNAStringSet
multiple_seqs <- dna_set[1:2]
print("Multiple sequences (first two):")
print(multiple_seqs)

# -----------------------------------------------------------------------------
# Section: Accessing Subsequences using getSeq
# Description: Demonstrate extracting subsequences using getSeq() and GRanges.
# -----------------------------------------------------------------------------

# Create a DNAStringSet with named sequences (e.g., simulating chromosomes)
g <- DNAStringSet(c(chr1 = "ACGTACGTACGT",
                    chr2 = "TGCACTGCA",
                    chr3 = "GATCGATC"))
print("Named DNA sequences:")
print(g)

# Extract a subsequence from 'chr1' spanning positions 1 to 4
s1 <- getSeq(g, GRanges("chr1", IRanges(1, 4)))
print("Subsequence from 'chr1' (positions 1 to 4):")
print(s1)

# Extract subsequences from 'chr2' using multiple IRanges
s2 <- getSeq(g, GRanges("chr2", IRanges(start = 1:4, end = 4:6)))
print("Subsequences from 'chr2' (multiple ranges):")
print(s2)

# -----------------------------------------------------------------------------
# Section: Sequence Analysis
# Description: Analyze sequence content by computing nucleotide and
# oligonucleotide frequencies.
# -----------------------------------------------------------------------------

# Print the original DNA sequence for analysis
print("DNA sequence for frequency analysis:")
print(dna_seq)

# Compute the frequency of each nucleotide (only A, C, G, T)
nucl_freq <- alphabetFrequency(dna_seq, baseOnly = TRUE)
print("Nucleotide frequency (A, C, G, T):")
print(nucl_freq)

# Compute the frequency of 4-mers (oligonucleotides of width 4)
oligo_freq <- oligonucleotideFrequency(dna_seq, width = 4)

# Sort the frequencies in descending order and retrieve the top 10 most frequent 4-mers
sorted_oligo_freq <- sort(oligo_freq, decreasing = TRUE)[1:10]
print("Top 10 most frequent 4-mers:")
print(sorted_oligo_freq)
