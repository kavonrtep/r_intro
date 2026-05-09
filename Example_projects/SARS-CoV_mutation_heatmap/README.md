# SARS-CoV-2 spike mutation heatmap

A small, single-question starter project. The full version is in `SARS-CoV_analysis_project/`; this directory tackles **one** of its questions in isolation.

> **Did key spike mutations evolve convergently across independent SARS-CoV-2 lineages?**

You are given a time-calibrated phylogeny of 400 SARS-CoV-2 genomes, per-tip metadata, and a pre-built table of spike-protein amino-acid substitutions for every tip. Your job is to render a heatmap of presence/absence for ~10 well-known spike mutations beside the tree, and to interpret which mutations recur on independent branches (the visual signature of convergent evolution).

See [`project_description.md`](project_description.md) for the full framing.

## Structure

- `data/` — pre-cached inputs (ready to use, no downloading needed)
  - `ncov_subsample.nwk` — Newick tree, branches in **years**, 400 tips
  - `ncov_metadata.tsv` — one row per tip
  - `spike_mutations_long.csv` — long-format mutation table: one row per `(strain, mutation)` pair, restricted to spike (`S:...`) substitutions
- `scripts/` — your R scripts (create as you go)
- `results/` — your CSV outputs (create as you go)
- `reports/` — your Quarto report and figures (create as you go)

## Tools and style

- **tidyverse** for data handling
- **ggplot2** + **ggtree** for the tree, **`gheatmap`** (in ggtree) for the heatmap
- **Quarto** for the report

Keep code short, clear, and functional.

## Where to start

Read [`project_description.md`](project_description.md) and [`example_prompt.md`](example_prompt.md), then ask your coding agent to draft an analysis plan from the latter.
