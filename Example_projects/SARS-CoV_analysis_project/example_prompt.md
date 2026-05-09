This is a workspace for phylogenetic-visualization data analysis based on the @project_description.md document. Inspect the document and create a plan for the analysis.

The main steps of analysis will include:

1. **Load and inspect data.** Read `data/ncov_subsample.nwk` and `data/ncov_metadata.tsv`.
   - Report number of tips, time span (`min`/`max` `date`), branch-length units (`max(node.depth.edgelength(tree))`).
   - Summarize metadata: counts per `Nextstrain_clade`, `region`, `country`, `host`.
   - Verify the join-key invariant: `setdiff(tree$tip.label, meta$strain)` must be empty.
   - Make a histogram of sampling dates and a barplot of clade counts.

2. **Build the time-calibrated base tree.** Construct `td <- full_join(as_tibble(tree), meta, by = c("label" = "strain")) |> as.treedata()` and produce the foundation plot `p_base` with `ggtree(td, mrsd = ..., as.Date = TRUE) + theme_tree2()`.
   - Use `mrsd = max(meta$decimal_date, na.rm = TRUE)` after converting `date` with `lubridate::ymd()` and `lubridate::decimal_date()`.

3. **Tip aesthetics.** Render `p1` = `p_base + geom_tippoint(aes(color = Nextstrain_clade))` with a viridis turbo palette. Then `p2` adding a second tip layer for `region`, demonstrating `ggnewscale::new_scale_color()`.

4. **Branch coloring via ancestral state reconstruction.** Run `ape::ace()` with `model = "ER"` on `Nextstrain_clade`, propagate node states with `groupOTU`, and render `p3` with branch colors. Document explicitly that this is a naive ASR (Morel et al. 2021) and that the colors should not be trusted as epidemiology.

5. **Spike mutation heatmap.** Build a tip × mutation presence/absence matrix for at least six key spike mutations (e.g. N501Y, E484K, D614G, K417N, L452R, T478K, P681H, P681R). Use `gheatmap(p1, mut_matrix, ...)`. Save the matrix to `results/spike_mutations.csv`.

6. **Multi-strip facet plot.** Use `geom_facet()` to add panels for region (tile), age (point), and genome length (segment). Apply `facet_widths()` last.

7. **MSA alongside the tree (optional).** If `data/spike_RBD_aligned.fasta` is present, render `msaplot(p1, fasta = rbd, window = c(...))`.

8. **Geography.** Two figures:
    - Country-fill world map joined to dominant clade per country.
    - `phytools::phylo.to.map` projection.
   Before drawing conclusions, compute `results/sequences_per_country.csv` and discuss sampling bias.

9. **Subclade zoom.** Use `tree_subset()` to isolate the Omicron clade and render it on its own time axis.

10. **Final integrated figure.** Combine three encodings (time tree + tip clade color + spike `gheatmap`) into one publication-quality panel. Save to `reports/figures/ncov_figure.pdf`.

11. **Critical-thinking essay.** 300 words. Pick one figure, identify a biological conclusion, and explain one way that conclusion could be wrong because of sampling, the inference method, or the visual encoding. Required reading: Morel et al. 2021 *MBE* 38:1777.

First, write a big-picture plan with the necessary implementation details. Divide the analysis into smaller controllable phases. Each phase should produce intermediate CSV/TSV outputs under `results/` and a Quarto-rendered section in `reports/`.

Write the plan to `analysis_plan.md`.

All necessary R packages are listed in `setup.md` and assumed to be installed:
- ape, treeio, ggtree, ggtreeExtra, ggplot2, ggnewscale, dplyr, tidyr, readr, lubridate, patchwork, phytools, Biostrings, quarto.

Data acquisition is documented in `nextstrain_howto.md` and is **already done** — read from `data/`, do not call any external API during the analysis.
