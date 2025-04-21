# Load Drosophila genome and annotation (UCSC dm6)
if (!requireNamespace("BSgenome.Dmelanogaster.UCSC.dm6", quietly=TRUE)) {
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
}
if (!requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene", quietly=TRUE)) {
    BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
}
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(GenomicRanges)
genome <- Dmelanogaster
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Define promoters as 500 bp upstream of transcript start sites
tx <- transcripts(txdb)                         # all transcripts with strand info
prom_regions <- promoters(tx, upstream=500, downstream=0)
prom_regions <- trim(prom_regions)             # trim to genome boundaries
prom_seqs <- getSeq(genome, prom_regions)

# Define background sequences: 500 bp sequences sampled from intergenic regions
set.seed(1)
# Sample random genomic positions for background (not overlapping known transcripts)
gene_regions <- genes(txdb)
intergenic <- gaps(gene_regions)               # intergenic regions as the gaps between genes
intergenic <- intergenic[strand(intergenic)=="*"]   # ignore strand for intergenic
# Sample 1000 random intergenic regions of width 500
bg_candidates <- intergenic[width(intergenic) >= 500]
bg_indices <- sample(length(bg_candidates), size=1000, replace=FALSE)
bg_regions <- resize(bg_candidates[bg_indices], width=500, fix="start")
bg_seqs <- getSeq(genome, bg_regions)

# Count TATAAA motif occurrences in promoters vs background
motif <- "TATAAA"
prom_hits <- vcountPattern(motif, prom_seqs)
bg_hits <- vcountPattern(motif, bg_seqs)
prop_prom_with_motif <- mean(prom_hits > 0)
prop_bg_with_motif <- mean(bg_hits > 0)

cat("Proportion of promoters with TATAAA motif:", prop_prom_with_motif)
cat("Proportion of background with TATAAA motif:", prop_bg_with_motif)