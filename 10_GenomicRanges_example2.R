# Analyze splice site motifs at introns of C. elegans
# - Install C. elegans genome and annotation
# - Get exons by gene and introns by transcript
# - extract intron sequences, keep only introns longer than 40 bp
# - Calculate consensus matrix for donor and acceptor sites
# - Plot sequence logos for donor and acceptor sites

# Load C. elegans genome and annotation (UCSC ce11)
if (!requireNamespace("BSgenome.Celegans.UCSC.ce11", quietly=TRUE)) {
    BiocManager::install("BSgenome.Celegans.UCSC.ce11")
}
if (!requireNamespace("TxDb.Celegans.UCSC.ce11.refGene", quietly=TRUE)) {
    BiocManager::install("TxDb.Celegans.UCSC.ce11.refGene")
}
library(BSgenome.Celegans.UCSC.ce11)
library(TxDb.Celegans.UCSC.ce11.refGene)
genome <- Celegans
txdb <- TxDb.Celegans.UCSC.ce11.refGene

# Get exons by gene and introns by transcript
exons_by_gene <- exonsBy(txdb, by="gene")
introns_by_tx <- intronsByTranscript(txdb, use.names=TRUE)

# Basic stats: number of exons per gene
exon_counts <- elementNROWS(exons_by_gene)  # this is equivalent to sapply(exons_by_gene, length)
plot(table(exon_counts), main="Exon counts per gene", xlab="Number of exons", ylab="Frequency")
# Length distributions
exon_lengths <- width(unlist(exons_by_gene))
hist(exon_lengths, breaks=1000, main="Exon lengths", xlab="Length (bp)", ylab="Frequency")
intron_lengths <- width(unlist(introns_by_tx))
hist(intron_lengths, breaks=1000, main="Intron lengths", xlab="Length (bp)", ylab="Frequency")

# Examine splice site motifs for introns
intron_seq <- getSeq(genome, unlist(introns_by_tx))
# remove short introns  shortet that 20 bp
intron_seq <- intron_seq[width(intron_seq) > 40]
donor_sites  <- subseq(intron_seq, start=1, width=20)
acceptor_sites <- subseq(intron_seq, end=-1, width=20)

# calculate consensus matrix for donor and acceptor sites
donor_consensus <- consensusMatrix(donor_sites, as.prob=TRUE, baseOnly=TRUE)
acceptor_consensus <- consensusMatrix(acceptor_sites, as.prob=TRUE, baseOnly=TRUE)

# optional plot logo for donor and acceptor sites
# optional: install ggseqlogo package if not already installed
if (!requireNamespace("ggseqlogo", quietly=TRUE)) {
    install.packages("ggseqlogo")
}
library(ggseqlogo)

ggseqlogo(donor_consensus, seq_type="dna")
ggseqlogo(acceptor_consensus, seq_type="dna")

ggseqlogo(list(
    donor = donor_consensus,
    acceptor = acceptor_consensus
), seq_type="dna")