
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


# If starting fresh, load the yeast genome and TxDb as in Use Case 1, then:
cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
cds_seqs <- extractTranscriptSeqs(genome, cds_by_tx)

# Translate DNA sequences to amino acid sequences
aa_seqs <- Biostrings::translate(cds_seqs)
# count sequences with stops (*" in the middle - they are not valid, occurance at the end is OK
not_valid <- grepl("\\*[A-Z]", as.character(aa_seqs))  # check for stop codons
table(not_valid)
aa_seqs <- aa_seqs[!not_valid]  # remove invalid sequences
print(length(aa_seqs))                # number of proteins (translated transcripts)


# Search for N-glycosylation motif N-X-[S/T] (exclude X = P)
# We'll use a regex with negative lookahead to exclude P after N.
aa_strings <- as.character(aa_seqs)  # convert AAStringSet to character vector
has_glyco_site <- grepl("N(?!P)[ST]", aa_strings, perl=TRUE)
glyco_proteins <- sum(has_glyco_site)
cat(sprintf("Proteins with N-linked glycosylation motif (N-X-[S/T]): %d (%.1f%%)\n",
            glyco_proteins, 100 * glyco_proteins / length(aa_strings)))
