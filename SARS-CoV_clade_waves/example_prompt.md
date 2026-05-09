This is a single-question starter project. Read `project_description.md`, then plan and implement the analysis below.

## The question

Did SARS-CoV-2 clades replace each other through time, or do they co-circulate for long periods?

## What to do

1. **Load the data.** Read `data/ncov_subsample.nwk` (Newick, 400 tips, branches in years) and `data/ncov_metadata.tsv` (one row per tip).

2. **Verify the join key.** Run `setdiff(tree$tip.label, meta$strain)` — must be empty. Run `max(node.depth.edgelength(tree))` — must be a few years (not a few thousandths).

3. **Summarize the data.** Count tips per `Nextstrain_clade`, and report the date range. Save the count table to `results/clade_counts.csv`.

4. **Build the base time-tree.** Convert `date` to `decimal_date` with `lubridate::ymd() |> lubridate::decimal_date()`. Make a treedata: `td <- full_join(as_tibble(tree), meta, by = c("label" = "strain")) |> as.treedata()`. Plot with `ggtree(td, mrsd = max(meta$decimal_date), as.Date = TRUE) + theme_tree2()` and a `scale_x_date()` with reasonable breaks.

5. **Color tips by clade.** Add `geom_tippoint(aes(color = Nextstrain_clade))` and a `scale_color_viridis_d(option = "turbo")` palette. Use `guide_legend(ncol = 2)` to keep the legend readable.

6. **(Optional) Clade abundance over time.** Bin tips by month × clade, plot as a stacked area or stacked bar with `ggplot(meta, aes(x = floor_date(ymd(date), "month"), fill = Nextstrain_clade)) + geom_bar()`. Save to `reports/figures/clade_waves.pdf`.

7. **Write a 200-word interpretation.** Look at the figure(s). Do clades come in waves, or do they overlap for a long time? Name at least three specific clade transitions you can see (e.g. "20A → 20I (Alpha) around early 2021"). Cite the figure.

8. **Render a Quarto report** under `reports/` that runs steps 1–7 end-to-end and produces the final figures and prose.

## Constraints

- Do not download anything — `data/` is already populated.
- Do not expand scope. The mutation heatmap and the geographic map are *not* part of this project. They live in sibling starter projects.
- Keep code short and teaching-friendly. Tidyverse + ggplot2 + ggtree + Quarto. Nothing exotic.
- Save tabular outputs as CSV under `results/`. Save the report and figures under `reports/`.

## What "done" looks like

- A rendered HTML report under `reports/` containing a time-calibrated tree with clade-colored tips, a clade count table, and a 200-word interpretation paragraph.
- `results/clade_counts.csv` exists.
- `setdiff(tree$tip.label, meta$strain)` returns `character(0)` somewhere visible in the report (so a reader can confirm the join is clean).

First, write a brief plan (≤ 1 page) to `analysis_plan.md`. Then implement.
