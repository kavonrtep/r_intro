# R Package Installation Instructions

All packages required by this course, organized by installation method.

## Step 0: Build dependencies

Required for compiling R/Bioconductor packages:

```bash
sudo apt install \
  build-essential r-base-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
  libpng-dev libjpeg-dev libtiff-dev \
  zlib1g-dev libgit2-dev cmake pkg-config libcairo2-dev
```

## Step 1: CRAN packages

### Option A: apt (Debian bookworm)

Verified against Debian bookworm main repository. The `bookworm-cran40` repo only provides base R — these packages come from the standard Debian repos.

```bash
sudo apt install \
  r-cran-ggplot2 r-cran-tidyverse r-cran-readxl \
  r-cran-gridextra r-cran-plotly r-cran-pheatmap r-cran-ggrepel \
  r-cran-ggforce r-cran-ape r-cran-optparse r-cran-knitr \
  r-cran-matrixstats r-cran-htmlwidgets r-cran-tibble r-cran-stringr \
  r-cran-phangorn r-cran-heatmaply r-cran-palmerpenguins \
  r-cran-ggseqlogo
```

`writexl` is not in Debian repos — install from R console: `install.packages("writexl")`

### Option B: R console (any platform)

```r
install.packages(c(
  "ggplot2", "tidyverse", "readxl", "writexl",
  "gridExtra", "plotly", "pheatmap", "ggrepel",
  "ggforce", "ape", "optparse", "knitr",
  "matrixStats", "htmlwidgets", "tibble", "stringr",
  "phangorn", "heatmaply", "palmerpenguins",
  "ggseqlogo"
))
```

## Step 2: Bioconductor packages (BiocManager)

Run in R console:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(version = "3.22")
install.packages("remotes")
remotes::install_version("tidytree", version = "0.4.6", repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "Biostrings",
  "GenomicRanges",
  "DESeq2",
  "BSgenome",
  "rtracklayer",
  "AnnotationDbi",
  "AnnotationHub",
  "Gviz",
  "phyloseq",
  "biomformat",
  "ggbio",
  "treeio",
  "ggtree",
  "trackViewer",
  "InteractionSet"
))
```

## Step 3: Genome and annotation data packages

These are large data packages (several hundred MB each). Run in R console:

```r
BiocManager::install(c(
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Drerio.UCSC.danRer11",
  "BSgenome.Scerevisiae.UCSC.sacCer3",
  "BSgenome.Celegans.UCSC.ce11",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
  "TxDb.Celegans.UCSC.ce11.refGene",
  "org.Hs.eg.db",
  "org.Dm.eg.db"
))
```

## Per-session package requirements

| Session | CRAN (apt) | Bioconductor (BiocManager) |
|---------|------------|---------------------------|
| 01 Introduction | (base R only) | — |
| 02 Data Import/Export | ggplot2, readxl, tibble | — |
| 03 Control/Apply | (base R only) | — |
| 04 Base Graphics | (base R only) | — |
| 05 ggplot2 | ggplot2, gridExtra | — |
| 06 DESeq2 | ggplot2, plotly, matrixStats, pheatmap, heatmaply, htmlwidgets | DESeq2, AnnotationDbi, org.Dm.eg.db |
| 07 Tidyverse | tidyverse, stringr, pheatmap | — |
| 08 RMarkdown/CLI | tidyverse, ggplot2, knitr, plotly, htmlwidgets, palmerpenguins, optparse | — |
| 09 Biostrings | — | Biostrings, BSgenome, GenomicRanges, rtracklayer, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Drerio.UCSC.danRer11 |
| 10 GenomicRanges | ggplot2, ggseqlogo | GenomicRanges, Biostrings, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Scerevisiae.UCSC.sacCer3, BSgenome.Celegans.UCSC.ce11, TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, TxDb.Celegans.UCSC.ce11.refGene |
| 11 Genomic Visualization | ape, tidyverse, ggplot2 | Gviz, trackViewer, GenomicRanges, rtracklayer, AnnotationHub, AnnotationDbi, InteractionSet, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg38.knownGene, BSgenome.Hsapiens.UCSC.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene, biomformat, phyloseq, ggtree |
| 12 Phylogenetics | ape, phangorn | treeio, ggtree |
