# Did key spike mutations evolve convergently across SARS-CoV-2 lineages?

## The single question

The SARS-CoV-2 spike protein is the entry receptor and the dominant target of neutralizing antibodies. A small set of spike mutations — `D614G`, `N501Y`, `E484K/A/Q`, `L452R`, `T478K`, `P681H/R`, `K417N/T`, `H655Y` — has appeared independently in multiple variant lineages. This is the textbook fingerprint of **convergent evolution under selection**: when the same mutation arises on independent branches of the tree, the most parsimonious explanation is that the mutation gives a fitness benefit large enough to be picked up repeatedly.

Convergence is *visible*. Place all tips on a phylogeny, paint a coloured cell beside each tip for every mutation it carries, and the colour patterns across the tree make repeat appearances obvious. Mutations that arose once and propagated form a single coherent block of cells; mutations that arose multiple times form **multiple separated blocks** at different parts of the tree.

Your task is to reproduce this picture.

## Working hypothesis

> A heatmap of ≈10 well-known spike substitutions, rendered as adjacent presence/absence cells next to a time-calibrated SARS-CoV-2 tree, will show at least one mutation (almost certainly `D614G`, plausibly `N501Y` and `L452R` as well) appearing in topologically separated parts of the tree — visual evidence of convergent evolution.

## Data

Three pre-cached files in `data/` — ready to use, no downloading required:

| File                          | What it is                                                                       |
|-------------------------------|----------------------------------------------------------------------------------|
| `ncov_subsample.nwk`          | Newick phylogeny of 400 tips. Branches are in years (time-calibrated).           |
| `ncov_metadata.tsv`           | Per-tip metadata. Joined to the tree by `strain`.                                 |
| `spike_mutations_long.csv`    | **Long-format** mutation table: one row per `(strain, mutation)`, S substitutions only. |

The long-format mutation table is provided so you can skip the most error-prone part of the wrangling (extracting per-spike mutations from a comma-separated text column). You still need to **pivot it to wide** for `gheatmap`.

Verify the join key invariant before doing anything else:

```r
setdiff(tree$tip.label, meta$strain)             # MUST be empty
setdiff(spike_long$strain, tree$tip.label)       # OK to be non-empty
```

The metadata columns relevant to this question are `strain`, `date`, `decimal_date`, and `Nextstrain_clade` (used to *colour* tips so you can read the heatmap rows in clade context).

## What to produce

A short Quarto report under `reports/` containing:

1. **Data summary.** Number of tips, number of unique spike mutations in `spike_mutations_long.csv`, and the **count of strains carrying each mutation** sorted descending. Save as `results/mutation_counts.csv`.

2. **Mutation matrix.** Pick 10 well-known spike mutations to track:

    ```r
    key_muts <- c("S:D614G", "S:N501Y", "S:E484K", "S:E484A", "S:L452R",
                  "S:T478K", "S:K417N", "S:K417T", "S:P681H", "S:P681R",
                  "S:H655Y")
    ```

   Build a tip × mutation **wide** presence/absence matrix from `spike_mutations_long.csv`:

    ```r
    mut_wide <- spike_long |>
      filter(mutation %in% key_muts) |>
      mutate(present = "present") |>
      pivot_wider(names_from = mutation, values_from = present,
                  values_fill = "absent") |>
      tibble::column_to_rownames("strain")
    ```

   Make sure every `key_mut` is a column even if entirely absent in this sample, every tip in the tree is a row even if it carries none of the mutations, and the column order matches `key_muts`. Save as `results/spike_mutations_wide.csv`.

3. **Figure 1 — heatmap on the tree.** Build a base `ggtree` with tips colored by `Nextstrain_clade` (so you can read the heatmap rows in clade context), then add the mutation matrix with `ggtree::gheatmap`:

    ```r
    p1 <- ggtree(td, mrsd = mrsd_value, as.Date = TRUE) +
      geom_tippoint(aes(color = Nextstrain_clade), size = 1.2) +
      scale_color_viridis_d(option = "turbo")

    p_heat <- gheatmap(p1, mut_wide, offset = 50, width = 0.5,
                       colnames_angle = 45, font.size = 3) +
      scale_fill_manual(values = c(present = "firebrick", absent = "grey92"),
                        name = "Spike mutation")
    ```

4. **A 200-word interpretation.** Look at the figure. For at least **two** of your `key_muts`, state whether you see *one* coherent block or *multiple separated blocks* on the tree, and what that means biologically. Specifically: **is there visible convergent evolution in this dataset**, and which mutation(s) show it most clearly?

## Pitfalls already burned (carry these forward)

- **Heatmap row order.** `gheatmap` aligns matrix rows to tip y-coordinates **by row name**, not by row order. Do **not** `arrange()` the matrix manually. Just make sure every `tree$tip.label` value exists as a row name (use a left-join trick or row-bind missing strains as all-`absent` rows).
- **Missing columns.** If a `key_mut` happens to be absent from every tip in this sample, `pivot_wider` will not create a column for it. Force every key mutation to appear as a column with all-`absent` values:

   ```r
   for (m in setdiff(key_muts, colnames(mut_wide))) mut_wide[[m]] <- "absent"
   mut_wide <- mut_wide[, key_muts]
   ```

- **Branch-length units.** `mrsd` only works if branches are in years. Verify with `max(node.depth.edgelength(tree))` — should be ≈ 6 for this dataset.
- **Reversion handling.** `spike_mutations_long.csv` already has reversions filtered out (a tip is listed as carrying a mutation only if the substitution is the one observed at the tip, after collapsing back-mutations along its lineage). You do not need to handle reversions yourself.

## Why this is a good first project

The full SARS-CoV-2 project layers six different visual encodings on the same tree. This starter focuses on the *most informationally dense* encoding — `gheatmap` — which is also the one most students get wrong on the first try (because of the row-name alignment trap). Master `gheatmap` here and the other encodings come more easily later.

## Going further

- Replace `Nextstrain_clade` tip color with `pango_lineage`, which is finer-grained — does the convergence pattern look different at the lineage level?
- Add a *second* heatmap stacked next to the first, e.g. for ORF1ab mutations. The pattern is `gheatmap(...) + new_scale_fill() + gheatmap(...)`.
- Compute a numeric **convergence score** per mutation: e.g. count distinct clades that carry it, or compute the parsimony score on the tree. Cross-check it against your visual reading.

## References

- Yu G. (2022). *Data Integration, Manipulation and Visualization of Phylogenetic Trees.* <https://yulab-smu.top/treedata-book/> — chapter 7 covers `gheatmap`.
- Harvey, W. T. et al. (2021). SARS-CoV-2 variants, spike mutations and immune escape. *Nat. Rev. Microbiol.* 19:409.
