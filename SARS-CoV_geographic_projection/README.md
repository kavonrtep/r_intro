# SARS-CoV-2 geographic projection

A small, single-question starter project. The full version is in `SARS-CoV_analysis_project/`; this directory tackles **one** of its questions in isolation.

> **How can we project a SARS-CoV-2 phylogeny onto a world map?**

You are given a time-calibrated phylogeny of 400 SARS-CoV-2 genomes, per-tip metadata including country and `latitude`/`longitude`, and a small per-country COVID-19 case-count snapshot. Your job is to render the same data in two complementary geographic encodings — a country-fill world map and a `phytools::phylo.to.map` projection that connects each tip to its sampling coordinates — and to interpret what they show (and what they hide).

See [`project_description.md`](./project_description.md) for the full framing.

## Structure

- `data/` — pre-cached inputs (ready to use, no downloading needed)
  - `ncov_subsample.nwk` — Newick tree, branches in **years**, 400 tips
  - `ncov_metadata.tsv` — one row per tip, includes `country`, `latitude`, `longitude`
  - `covid_cases_per_country.csv` — 233-row OWID snapshot, used in the optional sampling-bias extension
- `scripts/` — your R scripts (create as you go)
- `results/` — your CSV outputs (create as you go)
- `reports/` — your Quarto report and figures (create as you go)

## Tools and style

- **tidyverse** for data handling
- **ggplot2** + **maps** for the country-fill world map
- **phytools** for `phylo.to.map`
- **patchwork** for combining tree + map panels
- **Quarto** for the report

Keep code short, clear, and functional.

## Where to start

Read [`project_description.md`](./project_description.md) and [`example_prompt.md`](./example_prompt.md), then ask your coding agent to draft an analysis plan from the latter.
