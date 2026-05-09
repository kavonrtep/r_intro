# Setup

This is a small starter project. You need:

1. A coding agent (GitHub Copilot CLI or OpenCode) — see the full setup guide in `../SARS-CoV_analysis_project/setup.md`.
2. **VS Code** as an editor (recommended R + Quarto extensions).
3. **R packages** for tree handling, spatial plotting, and reporting.

## R packages

CRAN packages:

```r
install.packages(c(
  "ape", "treeio", "ggplot2", "dplyr", "tidyr", "readr",
  "tibble", "lubridate", "patchwork",
  "maps", "phytools"
))
```

Bioconductor (`ggtree`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
```

Quarto (CLI) for the report — install from <https://quarto.org/docs/get-started/>.

Verify:

```r
library(ape); library(treeio); library(ggtree)
library(ggplot2); library(maps); library(phytools); library(patchwork)
library(dplyr); library(tidyr); library(readr)
```

If all `library()` calls succeed without errors, you are ready.

## What you do **not** need

This starter does not use `ggtreeExtra`, `ggnewscale`, or `Biostrings`. Add them only if you do the optional extensions described in `project_description.md`.

## A note on `phytools`

`phytools` is a heavy dependency (it pulls in `mvtnorm`, `coda`, `clusterGeneration`, …) because it implements many comparative-methods algorithms beyond what we use here. Installation can take a few minutes the first time. If installation fails on a system without compilers, install via your OS package manager instead — on Ubuntu this is `r-cran-phytools` (older but adequate for `phylo.to.map`).
