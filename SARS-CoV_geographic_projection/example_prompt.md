This is a single-question starter project. Read `project_description.md`, then plan and implement the analysis below.

## The question

How can we project a SARS-CoV-2 phylogeny onto a world map, and what does the projection hide about sampling bias?

## What to do

1. **Load the data.** Read `data/ncov_subsample.nwk` (Newick, 400 tips, branches in years), `data/ncov_metadata.tsv` (one row per tip, has `country`, `latitude`, `longitude`, `Nextstrain_clade`), and `data/covid_cases_per_country.csv` (per-country case counts, 233 rows).

2. **Verify invariants.** `setdiff(tree$tip.label, meta$strain)` must be empty. `max(node.depth.edgelength(tree))` should be a few years.

3. **Summarize sampling distribution.** Count tips per country, sort descending, save the top of the table to `results/sequences_per_country.csv`. Report which 5 countries supply the most sequences in this 400-tip sample.

4. **Align country names** between the three data sources. The metadata uses `USA` and `Czech Republic`; OWID uses `United States` and `Czechia`. `ggplot2::map_data("world")` uses `USA` and `Czech Republic`. Plan one tidy step to normalize the OWID country names to match the metadata before any join. (Hint: `cases <- cases |> mutate(country = case_when(country == "United States" ~ "USA", country == "Czechia" ~ "Czech Republic", TRUE ~ country))`.)

5. **Figure 1 — country-fill choropleth.** Compute the dominant `Nextstrain_clade` per country (mode), join to `world <- ggplot2::map_data("world")` by country (`world$region` matches `meta$country`), and render:

    ```r
    p_map <- ggplot() +
      geom_polygon(data = world, aes(long, lat, group = group),
                   fill = "grey95", color = "grey60", linewidth = 0.1) +
      geom_polygon(data = world |>
                     inner_join(country_clade, by = c("region" = "country")),
                   aes(long, lat, group = group, fill = Nextstrain_clade),
                   color = "grey50", linewidth = 0.15) +
      scale_fill_viridis_d(option = "turbo",
                           guide = guide_legend(ncol = 2)) +
      coord_quickmap() + theme_void() +
      ggtitle("Dominant Nextstrain clade per country (in this 400-tip sample)")
    ```

   Pair it with a clade-colored `ggtree` of the same data using `patchwork`:

    ```r
    (p_tree | p_map) + plot_layout(widths = c(1, 1.3))
    ```

   Save to `reports/figures/tree_plus_map.pdf`.

6. **Figure 2 — `phytools::phylo.to.map`.** Build the coordinate matrix:

    ```r
    library(phytools)
    coords <- meta |>
      select(strain, latitude, longitude) |>
      filter(strain %in% tree$tip.label,
             !is.na(latitude), !is.na(longitude)) |>
      tibble::column_to_rownames("strain") |> as.matrix()

    common  <- intersect(rownames(coords), tree$tip.label)
    tree_g  <- ape::keep.tip(tree, common)
    coords  <- coords[tree_g$tip.label, c("latitude", "longitude")]

    obj <- phylo.to.map(tree_g, coords, plot = FALSE)
    ```

   Render and save to `reports/figures/phylo_to_map.pdf` (wide aspect ratio, this figure is hard to read otherwise):

    ```r
    pdf("reports/figures/phylo_to_map.pdf", width = 14, height = 7)
    plot(obj, direction = "rightwards", colors = "darkblue",
         cex.points = 0.4, ftype = "off", lwd = 0.25)
    dev.off()
    ```

7. **Write a 200-word interpretation.** Pick one specific apparent geographic conclusion the figures invite (e.g. "Asia is dominated by clade X", "no sampling occurred in Africa", …). Use `results/sequences_per_country.csv` to describe one specific way that conclusion could be wrong because of sampling bias. Cite numbers — *which countries supply the most sequences*, and how does that distort the figure?

8. **(Optional extension) Bias-correction figure.** Join your `sequences_per_country.csv` to `data/covid_cases_per_country.csv`, compute `sequences_per_million_cases`, and plot the top-15 countries as a bar chart. Save to `results/sampling_intensity.csv` and `reports/figures/sampling_intensity.pdf`. Discuss in the interpretation.

9. **Render a Quarto report** under `reports/` that runs steps 1–8 end-to-end.

## Constraints

- Do not download anything — `data/` is already populated.
- Do not expand scope. Heatmaps, ancestral state reconstruction, multi-strip facets, and clade-waves analyses live in sibling starter projects, not here.
- Tidyverse + ggplot2 + maps + phytools + patchwork + Quarto. Nothing exotic.
- Save tabular outputs as CSV under `results/`. Save the report and figures under `reports/`.

## What "done" looks like

- A rendered HTML report under `reports/` containing both figures and the 200-word sampling-bias interpretation.
- `results/sequences_per_country.csv` exists.
- The interpretation cites at least one specific country count from `results/sequences_per_country.csv` to support its claim about bias.

First, write a brief plan (≤ 1 page) to `analysis_plan.md`. Then implement.
