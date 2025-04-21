# Codon Usage Bias in Yeast (Saccharomyces cerevisiae)
# Not all synonymous codons (different triplets encoding the same amino acid) are used equally often – this phenomenon is known as codon usage bias
# Oganisms often prefer certain codons, which can correlate with tRNA abundance and affect translational efficienc
#
#
# Highly expressed genes tend to use optimal codons corresponding to abundant tRNAs.
# In this example, we examine codon usage in the yeast genome to see the bias in usage
# of synonymous codons.
# Approach:
# - Retrieve all protein-coding sequences (CDS) from the yeast genome using a TxDb annotation and BSgenome.
# - Count the frequency of each codon in the whole genome.
# - Identify bias by comparing frequencies of synonymous codons (e.g. which codon for the same amino acid is most used).
#   We will also check which stop codon is most common.



# Load yeast genome and annotation (Saccharomyces cerevisiae S288C, UCSC sacCer3)
if (!requireNamespace("BSgenome.Scerevisiae.UCSC.sacCer3", quietly=TRUE)) {
    BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
}

# part of Bioconductor repository are also corresponding TxDb objects - which are genome annotation databases
if (!requireNamespace("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", quietly=TRUE)) {
    BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
}
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(Biostrings)
genome <- Scerevisiae  # BSgenome object
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

# Extract all CDS (coding DNA sequences) by transcript
# The cdsBy function retrieves the coding sequences for each transcript (by='tx') in the TxDb object.
cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
# Get the spliced CDS sequences for each transcript
cds_seqs <- extractTranscriptSeqs(genome, cds_by_tx)

# Compute codon frequencies across all CDS
codon_freq <- trinucleotideFrequency(cds_seqs, step=3, simplify.as="collapsed")
codon_freq <- sort(codon_freq, decreasing=TRUE)
head(codon_freq)  # show the six most frequent codons

# Translate codon names to amino acids and sum counts per amino acid
aa_for_codon <- GENETIC_CODE[names(codon_freq)]        # map each codon to its amino acid (including "*" for stop)
# Note tapply() function is used to apply a function (sum) to subsets of data (codon_freq) defined by another vector (aa_for_codon).

# Calculat the usage of each amino acid
aa_totals <- tapply(codon_freq, aa_for_codon, sum)     # total counts per amino acid
aa_totals <- sort(aa_totals, decreasing=TRUE)
print(aa_totals)

# For example codons for Leucine (L) are "TTA" "TTG" "CTA" "CTT" "CTG" "CTC"
codons_for_L <- names(aa_for_codon[aa_for_codon == "L"])
codons_for_L

# codon frequencies for Leucine (L)
codon_freq[codons_for_L]

# codon frequencies for Stop codons
stop_codons <- names(aa_for_codon[aa_for_codon == "*"])
codon_freq[stop_codons]

# TASK - Calculate codon usage for all amino acids
codon_usage <- list() # as relative codon usage

for (AA in names(AMINO_ACID_CODE)){
    codons <- names(aa_for_codon[aa_for_codon == AA])
    codon_usage[[AA]] <- codon_freq[codons]
    # calculate observed vs expected codon usage
    # expected codon usage is the frequency of the codon divided by the number of synonymous codons
    codon_usage[[AA]] <- codon_usage[[AA]]/sum(codon_usage[[AA]])/(1/length(codons))
}


# TASK - are all start codons ATG?