# ## Overview of BSgenome
#
# - **Purpose:**  Infrastructure for full reference genomes in R via Bioconductor
# - **Data classes:**  Stores genomes as `DNAString` objects in on‑disk packages
# - **Integration:**  Works seamlessly with Biostrings, GenomicRanges, and other Bioconductor tools

# ## Installation & Loading

# - **Install BSgenome infrastructure**  via BiocManager

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome")

# - **Install specific genome packages**

BiocManager::install(c(
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Drerio.UCSC.danRer11"
))

# - **Load packages into session**

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Drerio.UCSC.danRer11)

# ## Discovering Available Genomes
# - **List available BSgenome packages**  for your Bioconductor version:

available.genomes()

# - **List installed BSgenome packages**  on your system:

installed.genomes()

# ## Inspecting Human Genome (hg38)

# - **Attach and assign**

library(BSgenome.Hsapiens.UCSC.hg38)
genome_human <- BSgenome.Hsapiens.UCSC.hg38

# - **View sequence names**

seqnames(genome_human) |>
  head(10)

# - **Get chromosome lengths**

seqlengths(genome_human)["chr1"]

# ## Basic Genome Queries

# - **Extract subsequence**  (e.g., chr1 positions 100000–110000):

chr1_10kb <- getSeq(genome_human, "chr1", 100000, 110000)
print(chr1_10kb)

# - Accessing entire chromosome sequences:

chr1 <- genome_human$chr1
print(chr1)

# ## Basic Genome Queries

metadata(genome_human)

