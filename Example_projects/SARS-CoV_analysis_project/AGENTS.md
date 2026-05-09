# SARS-CoV-2 phylogeny project — instructions for the coding agent

This directory is an educational example project for using a coding agent to drive a phylogenetic-visualization workflow. Keep solutions simple, readable, and functional. Prefer straightforward scripts and explicit steps over abstraction-heavy designs.

## Project goals

- Render a SARS-CoV-2 time-calibrated tree together with rich Nextstrain-style metadata.
- Produce the four required figures (time-tree, branch-colored tree, mutation `gheatmap`, multi-strip facet plot) plus optional MSA, geographic, and integrated-final-figure variants.
- Write a 300-word critical-thinking essay on one figure's biological claim and one way that claim could be wrong.
- Favor teaching value, reproducibility, and maintainability over clever code.

## Expected project structure

Use and preserve this layout:

- `data/` - input data: pre-cached Newick tree, metadata TSV, optional RBD alignment FASTA
- `scripts/` - R scripts for tree loading, joining, and per-figure rendering
- `results/` - intermediate tabular outputs (CSV/TSV)
- `reports/` - rendered Quarto report and figure files

If a directory is missing and is needed by the workflow, create it with that exact name.

## Data acquisition

Data files are **pre-cached** in `data/`. Do not call Nextstrain, GISAID, or any other external API during the analysis. Acquisition details live in:

- `./nextstrain_howto.md` — source of truth for how the cached files were obtained, what columns they contain, and how to refresh them.

If `data/` is empty or inconsistent, stop and ask before refetching.

## R coding style

Write R code in a teaching-friendly style:

- Prefer the **tidyverse** for data import, joins, and reshaping
- Prefer **ggplot2** + **ggtree** family for plotting
- Prefer **Quarto** for the final report
- Keep code short and easy to follow
- Use clear variable names: `tree`, `meta`, `td` for the joined `treedata`, `p_base`, `p1`, `p2`, …
- Split work into small functions only when it improves readability
- Avoid unnecessary metaprogramming, deep nesting, and overly generic helper layers
- Do not introduce complexity just to make the code look advanced

## Implementation guidance

When adding or editing analysis code:

- Assume scripts are run from the project root
- Use relative paths such as `data/...`, `results/...`, and `reports/...`
- Save tabular outputs as CSV/TSV
- Save figures as PDF or PNG under `reports/figures/`; set `width` and `height` explicitly
- Prefer one script per logical step (load → join → time-tree → tip color → branch color → heatmap → facet → MSA → map → final)
- Add brief comments only where they help a learner understand the step

## Reporting guidance

When creating reports:

- Use a Quarto document under `reports/`
- Keep narrative concise and focused on biological or analytical interpretation
- Show code and outputs in a way that is easy for students to inspect
- The report should reproduce the four required figures, plus the critical-thinking essay

## Pitfalls to avoid (already burned in prior runs)

- **Tip-label mismatch.** Always run `setdiff(tree$tip.label, meta$strain)` before joining. Nextstrain strain names contain `/` and `|` — tree-builders sometimes mangle them.
- **Heatmap row order.** `gheatmap` aligns by row *name*, not row order. Don't `arrange()` the matrix.
- **Forgetting `new_scale_fill()` / `new_scale_color()`** between two layers that use the same aesthetic — produces a single merged scale that maps both variables to the same colors.
- **Branch length units.** `mrsd` requires branches in time. Verify with `max(node.depth.edgelength(tree))` — should be ≈3–5 for a multi-year tree, not ≈0.001.
- **`facet_widths()` last.** It operates on the assembled gtable, so it must be the final call in the chain.
- **Naive ASR.** `ape::ace(..., model = "ER")` is fast but ignores sampling bias and assumes time-reversibility. State this caveat in any branch-color figure.
- **Geography is misleading.** Country of sampling ≠ country of infection origin. Pair any map with sequences-per-country vs. cases-per-country before drawing conclusions.

## What to avoid

- Overly complex object-oriented designs
- Large frameworks or unnecessary package dependencies
- Hidden side effects
- Clever but hard-to-read one-liners
- Live API calls during analysis (acquisition is offline, see `nextstrain_howto.md`)
- Writing code that is optimized for scale instead of clarity

# Ignores

Ignore these files and directories: notes
