# Did SARS-CoV-2 clades replace each other through time?

## The single question

Public sequence repositories contain millions of SARS-CoV-2 genomes annotated with sampling date and lineage. The Nextstrain project assigns every genome to a **clade** — a coarse-grained label like `19A`, `20I (Alpha)`, `21K (Omicron)`, `25C (XFG)` — that summarizes its evolutionary position.

When the genomes are placed on a **time-calibrated phylogeny** and tips are colored by clade, a striking pattern usually appears: clades come in *waves*. A clade rises, dominates the sample for a few months, then is displaced by a new clade. The clade transitions are not gradual — they look like discrete sweeps.

Your task is to reproduce this picture from the data shipped in `data/`, and to interpret what it shows.

## Working hypothesis

> SARS-CoV-2 clades succeed each other in time-limited waves rather than co-existing for long periods. A figure that places tips on a calendar x-axis and colors them by `Nextstrain_clade` will reveal this pattern as horizontally banded coloring on the tree.

## Data

Two pre-cached files in `data/` — ready to use, no downloading required:

| File                  | What it is                                                            |
|-----------------------|-----------------------------------------------------------------------|
| `ncov_subsample.nwk`  | Newick phylogeny of 400 tips. **Branches are in years** (time-calibrated). |
| `ncov_metadata.tsv`   | Per-tip metadata. One row per tip, joined to the tree by `strain`.    |

Verify the join key invariant before doing anything else:

```r
setdiff(tree$tip.label, meta$strain)   # MUST be empty
```

The metadata columns relevant to this question are:

- `strain` — join key
- `date` — sampling date (`YYYY-MM-DD`)
- `decimal_date` — same date as a numeric, used as `mrsd` anchor for `ggtree`
- `Nextstrain_clade` — the categorical label you will color by
- `pango_lineage` — finer-grained Pango notation (optional)

Other columns (`region`, `host`, `length`, …) are present but not needed for this question.

## What to produce

A short Quarto report under `reports/` containing:

1. **Data summary.** Number of tips, date range, number of distinct clades, and the count of tips per clade. Save the count table as `results/clade_counts.csv`.
2. **Figure 1 — base time-tree.** A `ggtree` plot of the tree with `mrsd` set to the most recent sampling date and `as.Date = TRUE`, so the x-axis is calendar time. No tip color yet — this is the foundation.
3. **Figure 2 — tips colored by clade.** Same tree, with `geom_tippoint(aes(color = Nextstrain_clade))`. This is the figure that answers the question.
4. **Figure 3 — clade abundance over time** (optional but encouraged). Bin tips by month and clade; plot as a stacked bar or area chart. This makes the wave pattern visible without the tree.
5. **A 200-word interpretation paragraph** answering: *do clades replace each other in waves, or do they co-circulate?* Use the figure as evidence and name at least three specific clade transitions you can see.

## Pitfalls already burned (carry these forward)

- **Tip-label mismatch.** Always run `setdiff(tree$tip.label, meta$strain)` before joining. Should be empty.
- **Branch-length units.** `mrsd` only works if branches are in years. Verify with `max(node.depth.edgelength(tree))` — should be ≈ 6 for this dataset, *not* ≈ 0.001.
- **Too many clades for one palette.** This sample has ~15 distinct clades. `viridis::turbo` handles that range better than `Set1`. Use `guide_legend(ncol = 2)` to keep the legend readable.

## Why this is a good first project

Every figure you make builds on `ggtree(td, mrsd = ..., as.Date = TRUE)`. Once you have that one expression working, layering the metadata on top is just `+ geom_tippoint(aes(color = ...))`. The harder techniques — heatmaps, ancestral state reconstruction, geographic projection — come later in the full project.

## Going further

If you finish early and want to extend:

- Add a second tip layer encoding `region` with `ggnewscale::new_scale_color()`.
- Compute and report **clade duration** — first and last sampling date of each clade — as a second `results/` table.
- Try the same plot in `layout = "fan"` for comparison.

These are exactly the kind of additions you would make in the full `SARS-CoV_analysis_project/`.

## Reference

Hadfield et al. (2018). Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics* 34:4121.
