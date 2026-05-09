# Clade-waves starter — instructions for the coding agent

This is a single-question starter project. Solve **one** question well; do not expand scope.

## The question

> Did SARS-CoV-2 clades replace each other through time?

The deliverable is a Quarto report containing a time-calibrated phylogeny with tips colored by `Nextstrain_clade`, plus 200 words of interpretation. See `project_description.md` for full framing.

## Project goals

- Render a clean, time-axis tree with tip color encoding `Nextstrain_clade`.
- Produce a per-clade count table and (optionally) a clade-by-month abundance plot.
- Answer the single question in plain prose, citing the figure.

## Expected layout

- `data/` — already populated, do not modify
- `scripts/` — your R scripts (create as you go)
- `results/` — your CSV outputs (create as you go)
- `reports/` — your Quarto report and figures (create as you go)

If a directory is missing and is needed, create it with that exact name.

## Data

`data/ncov_subsample.nwk` and `data/ncov_metadata.tsv` are already in the project. Read them with `ape::read.tree()` and `readr::read_tsv()`. Do **not** download anything; everything you need is already on disk.

## R coding style

- Tidyverse for data wrangling, ggplot2 + ggtree for plotting, Quarto for the report.
- Short, clear, teaching-friendly code. One script per logical step.
- Use relative paths from the project root: `data/...`, `results/...`, `reports/...`.

## Pitfalls to avoid

- **`mrsd` requires branches in years.** Check `max(node.depth.edgelength(tree))` — should be a few years, *not* a few thousandths.
- **Join key invariant.** Always run `setdiff(tree$tip.label, meta$strain)` before binding metadata to the tree.
- **Don't try to do everything.** This is the clade-waves question only. Not the mutation heatmap, not the geographic map, not the ancestral state reconstruction. Those live in the full `SARS-CoV_analysis_project/` and the sibling `SARS-CoV_mutation_heatmap/`.

## Scope discipline

If you find yourself writing helper functions, generating multiple analysis pipelines, or planning extensions before the core figure works, stop. Get the time-tree with clade-colored tips rendering first. Iterate from there.
