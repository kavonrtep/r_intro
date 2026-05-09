# Setup

This is a small starter project. You need:

1. A coding agent (GitHub Copilot CLI or OpenCode) — see the full setup guide in `../SARS-CoV_analysis_project/setup.md`.
2. **VS Code** as an editor (recommended R + Quarto extensions).
3. **R packages** for tree handling, plotting, and reporting.

## R packages

CRAN packages:

```r
install.packages(c(
  "ape", "treeio", "ggplot2", "dplyr", "tidyr", "readr",
  "tibble", "lubridate"
))
```

Bioconductor (`ggtree` provides both `ggtree()` and `gheatmap()`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
```

Quarto (CLI) for the report — install from <https://quarto.org/docs/get-started/>.

Verify:

```r
library(ape); library(treeio); library(ggtree)
library(ggplot2); library(dplyr); library(tidyr); library(readr); library(lubridate)
```

If all `library()` calls succeed without errors, you are ready.

## What you do **not** need

This starter does not use `ggtreeExtra`, `ggnewscale`, `phytools`, or `Biostrings`. Add them only if you do the optional extensions described in `project_description.md`.
