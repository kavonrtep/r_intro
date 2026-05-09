# Assignment 7: Convergent Evolution of SARS-CoV-2 Spike Mutations

## Background

The SARS-CoV-2 spike protein mediates cell entry and is the main target of neutralising antibodies. A notable feature of the pandemic was that certain spike amino-acid substitutions — such as N501Y, E484K, and L452R — appeared independently in multiple lineages (Alpha, Beta, Delta, Omicron) without being inherited from a common ancestor. When the same mutation arises on separate branches of a phylogeny, it is evidence of convergent evolution: the mutation confers a fitness advantage strong enough to be picked up repeatedly under independent selective pressures.

The visual test for convergence is straightforward. If you align a presence/absence matrix of mutations to a phylogenetic tree — one row per tip, one column per mutation — a mutation that spread clonally appears as a single coherent block of coloured cells. A mutation that arose independently multiple times appears as several separate blocks scattered across different parts of the tree.

In this assignment you will build that figure for a set of well-known SARS-CoV-2 spike substitutions.

## Data

All files are in `data_SARS_CoV/`:

| File | Description |
|---|---|
| `ncov_subsample.nwk` | Newick phylogeny, 400 tips. Branches are in **years** (time-calibrated). |
| `ncov_metadata.tsv` | One row per tip. Key columns: `strain`, `date`, `decimal_date`, `Nextstrain_clade`, `pango_lineage`. |
| `spike_mutations_long.csv` | Long-format table: one row per `(strain, mutation)` pair. Only spike (`S:...`) substitutions; reversions already collapsed. |

The join key between all three files is `strain` (tree tip labels = `strain` in the metadata = `strain` in the mutation table).

## Output

Produce a **self-contained Quarto HTML document** containing all code, figures, and written answers. Use `embed-resources: true` in the YAML.

This assignment is expected to be completed with a coding agent — approach it the same way as the Proteome Analysis project (`../Proteome_analysis_project/`), where the agent handles implementation while you direct the analysis and interpret the results. You are responsible for the biological reasoning and the final interpretation.

## Useful packages

```r
install.packages(c("ape", "ggplot2", "dplyr", "tidyr", "readr", "tibble", "lubridate"))
BiocManager::install(c("treeio", "ggtree"))
```

Relevant functions to look up: `ape::read.tree()`, `ggtree::ggtree()`, `ggtree::gheatmap()`, `tidyr::pivot_wider()`, `dplyr::full_join()`, `treeio::as.treedata()`.

---

## Task 1: Load and explore the data

Load all three files and explore their structure. Report how many tips the tree has, what the metadata contains, how many unique spike mutations appear in the long table, and which mutations are most prevalent.

---

## Task 2: Build the mutation matrix

Select a set of well-known spike substitutions to track — consider using:

```
S:D614G, S:N501Y, S:E484K, S:E484A, S:L452R,
S:T478K, S:K417N, S:K417T, S:P681H, S:P681R, S:H655Y
```

From the long mutation table, build a presence/absence matrix suitable for use with `gheatmap`: rows are tip names, columns are mutations, values indicate whether each tip carries each mutation.

---

## Task 3: Plot the time-calibrated tree

Plot the phylogeny as a time-calibrated tree. Attach the metadata and colour tips by `Nextstrain_clade` so the heatmap rows can be read in clade context.

---

## Task 4: Add the mutation heatmap

Attach the mutation matrix to the tree plot using `gheatmap`. Make the figure legible — adjust offset, width, and colour scale as needed.

---

## Task 5: Interpret the results

Write a short summary of what the heatmap shows. For at least two of the mutations, describe whether the pattern is clonal (one block) or convergent (multiple separate blocks), name the lineages involved, and explain what the pattern means biologically.

Upload the rendered `.html` file as your submission.
