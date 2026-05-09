# SARS-CoV-2 clade waves

A small, single-question starter project. The full version is in `SARS-CoV_analysis_project/`; this directory tackles **one** of its questions in isolation.

> **Did SARS-CoV-2 clades replace each other through time?**

You are given a time-calibrated phylogeny of 400 SARS-CoV-2 genomes sampled between January 2020 and April 2026, plus per-tip metadata including the Nextstrain clade and sampling date. Your job is to render a figure that makes the temporal succession of clades visible, and to interpret what it shows.

See [`project_description.md`](./project_description.md) for the full framing.

## Structure

- `data/` — pre-cached inputs (ready to use, no downloading needed)
  - `ncov_subsample.nwk` — Newick tree, branches in **years**, 400 tips
  - `ncov_metadata.tsv` — one row per tip
- `scripts/` — your R scripts (create as you go)
- `results/` — your CSV outputs (create as you go)
- `reports/` — your Quarto report and figures (create as you go)

## Tools and style

- **tidyverse** for data handling
- **ggplot2** + **ggtree** for plotting
- **Quarto** for the report

Keep code short, clear, and functional.

## Where to start

Read [`project_description.md`](./project_description.md) and [`example_prompt.md`](./example_prompt.md), then ask your coding agent to draft an analysis plan from the latter.
