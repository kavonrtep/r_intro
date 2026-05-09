This is a single-question starter project. Read `project_description.md`, then plan and implement the analysis below.

## The question

Did key spike mutations evolve convergently across SARS-CoV-2 lineages? Specifically: do any of the 10 well-known spike substitutions appear on **multiple separated parts of the tree** in this 400-tip sample?

## What to do

1. **Load the data.** Read `data/ncov_subsample.nwk`, `data/ncov_metadata.tsv`, and `data/spike_mutations_long.csv`. The third file is in **long format**: one row per `(strain, mutation)` pair, restricted to spike (`S:...`) substitutions.

2. **Verify invariants.** `setdiff(tree$tip.label, meta$strain)` must be empty. `max(node.depth.edgelength(tree))` should be a few years.

3. **Inspect the long table.** Count how many strains carry each mutation (`spike_long |> count(mutation, sort = TRUE)`). Save the top 30 to `results/mutation_counts.csv`. Do `D614G` and the Omicron-defining mutations (`H655Y`, `K417N`, `N501Y`, …) dominate?

4. **Pick the key mutations.** Use:

    ```r
    key_muts <- c("S:D614G", "S:N501Y", "S:E484K", "S:E484A", "S:L452R",
                  "S:T478K", "S:K417N", "S:K417T", "S:P681H", "S:P681R",
                  "S:H655Y")
    ```

5. **Build the wide matrix.** `pivot_wider` from `spike_mutations_long.csv` to get one row per strain and one column per mutation, with `"present"` / `"absent"` values:

    ```r
    mut_wide <- spike_long |>
      filter(mutation %in% key_muts) |>
      mutate(present = "present") |>
      pivot_wider(names_from = mutation,
                  values_from = present,
                  values_fill = "absent") |>
      tibble::column_to_rownames("strain")

    # Ensure every key_mut is a column
    for (m in setdiff(key_muts, colnames(mut_wide))) mut_wide[[m]] <- "absent"
    mut_wide <- mut_wide[, key_muts]

    # Ensure every tip is a row (tips with no key mutations were dropped by the pivot)
    missing_tips <- setdiff(tree$tip.label, rownames(mut_wide))
    if (length(missing_tips)) {
      pad <- as.data.frame(matrix("absent", nrow = length(missing_tips),
                                  ncol = ncol(mut_wide),
                                  dimnames = list(missing_tips, colnames(mut_wide))))
      mut_wide <- rbind(mut_wide, pad)
    }
    ```

   Save as `results/spike_mutations_wide.csv` (with row names as a `strain` column for legibility).

6. **Build the base tree** with tips colored by clade so the heatmap rows are readable in clade context:

    ```r
    meta <- meta |> mutate(decimal_date = lubridate::decimal_date(lubridate::ymd(date)))
    td <- full_join(as_tibble(tree), meta, by = c("label" = "strain")) |>
      treeio::as.treedata()

    p1 <- ggtree(td, mrsd = max(meta$decimal_date), as.Date = TRUE,
                 size = 0.25, color = "grey40") +
      geom_tippoint(aes(color = Nextstrain_clade), size = 1.2) +
      scale_color_viridis_d(option = "turbo", name = "Clade",
                            guide = guide_legend(ncol = 2)) +
      theme_tree2()
    ```

7. **Add the heatmap:**

    ```r
    p_heat <- gheatmap(p1, mut_wide,
                       offset = 50, width = 0.5,
                       colnames_angle = 45, colnames_offset_y = -2,
                       font.size = 3) +
      scale_fill_manual(values = c(present = "firebrick", absent = "grey92"),
                        name = "Spike mutation")
    ```

   Save the figure to `reports/figures/spike_heatmap.pdf` at sensible width/height.

8. **Write a 200-word interpretation.** For at least two of your `key_muts`, state whether the cells form **one coherent block** (single origin, propagated by descent) or **multiple separated blocks** (convergent evolution). Specifically discuss `D614G` (expected: near-fixed across all post-2020 tips, single block at top of tree), `N501Y` (expected: separate blocks in Alpha and in Omicron descendants), and one mutation of your choice.

9. **Render a Quarto report** under `reports/` that runs steps 1–8 end-to-end.

## Constraints

- Do not download anything — `data/` is already populated.
- Do not expand scope. Time-tree-by-clade-only and geographic-only analyses live in sibling starter projects, not here.
- Tidyverse + ggplot2 + ggtree + Quarto. Nothing exotic.
- Save tabular outputs as CSV under `results/`. Save the report and figures under `reports/`.

## What "done" looks like

- A rendered HTML report under `reports/` containing the `gheatmap` figure and the 200-word interpretation.
- `results/mutation_counts.csv` and `results/spike_mutations_wide.csv` exist.
- The interpretation explicitly names at least one mutation as showing a single-block (clonal-descent) pattern and at least one as showing a multi-block (convergent) pattern.

First, write a brief plan (≤ 1 page) to `analysis_plan.md`. Then implement.
