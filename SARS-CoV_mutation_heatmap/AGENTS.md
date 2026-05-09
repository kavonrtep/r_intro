# Mutation-heatmap starter — instructions for the coding agent

This is a single-question starter project. Solve **one** question well; do not expand scope.

## The question

> Did key spike mutations evolve convergently across SARS-CoV-2 lineages?

The deliverable is a Quarto report containing a `gheatmap` of ~10 well-known spike mutations beside a time-calibrated tree, plus 200 words of interpretation. See `project_description.md` for full framing.

## Project goals

- Build a tip × mutation **wide** presence/absence matrix from the long-format CSV provided in `data/`.
- Render it as a `gheatmap` next to a time-calibrated tree with tips colored by `Nextstrain_clade`.
- Identify visually which mutations show the convergent-evolution pattern (multiple separated blocks on the tree).

## Expected layout

- `data/` — already populated, do not modify
- `scripts/` — your R scripts (create as you go)
- `results/` — your CSV outputs (create as you go)
- `reports/` — your Quarto report and figures (create as you go)

If a directory is missing and is needed, create it with that exact name.

## Data

Three files are already in `data/`:

- `ncov_subsample.nwk` — 400-tip Newick tree, branches in years.
- `ncov_metadata.tsv` — per-tip metadata.
- `spike_mutations_long.csv` — pre-built long-format mutation table (one row per `(strain, mutation)` pair, S substitutions only, reversions already collapsed).

Read with `ape::read.tree()`, `readr::read_tsv()`, `readr::read_csv()`. Do **not** download anything; everything you need is on disk.

## R coding style

- Tidyverse for wrangling, `ggtree` + `gheatmap` for the plot, Quarto for the report.
- Short, clear, teaching-friendly code. One script per logical step.
- Use relative paths from the project root: `data/...`, `results/...`, `reports/...`.

## Pitfalls to avoid (these are the ones that bite first-time `gheatmap` users)

- **`gheatmap` aligns by row NAME, not row order.** Never `arrange()` the matrix manually before passing it to `gheatmap`. Just make sure every tip in the tree appears as a row name, and every `key_mut` appears as a column.
- **Missing columns.** If a `key_mut` is absent from every tip in the sample, `pivot_wider` will silently drop the column. Force-create empty columns:

   ```r
   for (m in setdiff(key_muts, colnames(mut_wide))) mut_wide[[m]] <- "absent"
   mut_wide <- mut_wide[, key_muts]
   ```

- **Tips with no mutations** are also dropped by the pivot. They need to appear as all-`absent` rows in the matrix or they will not align with the tree.
- **`mrsd` requires branches in years.** Check `max(node.depth.edgelength(tree))` — should be a few years, *not* a few thousandths.
- **Join key invariant.** Always run `setdiff(tree$tip.label, meta$strain)` before binding metadata to the tree.

## Scope discipline

Do **not** also build the geographic projection, the ancestral state reconstruction, the multi-strip facet plot, or any of the other encodings from the full project. Those live in `../SARS-CoV_analysis_project/` and a sibling `../SARS-CoV_clade_waves/`. Keep this project focused on `gheatmap`.

If you finish the core figure quickly, the optional extensions in `project_description.md` are the right next step — not adding new encodings.
