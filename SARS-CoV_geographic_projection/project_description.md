# How can we project a SARS-CoV-2 phylogeny onto a world map?

## The single question

Every tip in the SARS-CoV-2 sample is annotated with a sampling **country** and a pair of `(latitude, longitude)` coordinates. A phylogeny is one geometry; the globe is another. There are several principled ways to bind one to the other:

1. **Country-fill choropleth.** Compute one summary statistic per country (e.g. dominant clade, number of sequences, or earliest sampling date), then fill each country polygon by that value. The map is the picture; the tree never appears.
2. **Tree-to-map projection.** Draw the tree on one side and a world map on the other. Connect each tip with a line to the geographic point at its sampling coordinates. Both geometries are visible side-by-side. `phytools::phylo.to.map` is the canonical implementation.
3. **Spatial-tree composite.** Place the tree on the map directly (each tip plotted at its sampling coordinates, edges drawn as great-circle arcs). Visually striking; computationally heavier; outside the scope of this starter.

Your task is to produce **(1)** and **(2)** for the same dataset, and to interpret what each makes visible — and, equally important, what each obscures.

## Working hypothesis

> A naive geographic visualization of a SARS-CoV-2 phylogeny will appear to show clear geographic structure (e.g. "this clade is European"), but the apparent structure is largely an artifact of *sampling effort*, not of biology. Two countries with identical case burdens but very different sequencing capacity will look very different on the map. A critical reader of any phylogeographic figure must check the sampling distribution before reading off geographic conclusions.

## Data

Three pre-cached files in `data/` — ready to use, no downloading required:

| File                              | What it is                                                                                |
|-----------------------------------|-------------------------------------------------------------------------------------------|
| `ncov_subsample.nwk`              | Newick phylogeny of 400 tips. Branches are in years (time-calibrated).                    |
| `ncov_metadata.tsv`               | Per-tip metadata, joined to the tree by `strain`. Includes `country`, `latitude`, `longitude`. |
| `covid_cases_per_country.csv`     | Per-country COVID-19 case counts (233 countries, snapshot from Our World in Data).        |

Verify the join key invariant before doing anything else:

```r
setdiff(tree$tip.label, meta$strain)   # MUST be empty
```

The metadata columns relevant to this question:

- `strain` — join key
- `country` — sampling country (most relevant geographic resolution)
- `division` — sub-national region (state/province), present but not always populated
- `latitude`, `longitude` — coordinates of `country` (centroid-ish, supplied by Nextstrain)
- `Nextstrain_clade` — used as a categorical fill in the country-level map

The case-count file uses **slightly different country names** for two cases: `United States` (vs. metadata's `USA`) and `Czechia` (vs. metadata's `Czech Republic`). Plan a small `mutate(country = case_when(...))` step to align the names before joining; this is a normal real-life nuisance.

## What to produce

A short Quarto report under `reports/` containing:

1. **Data summary.** Number of tips, number of distinct countries, distribution of tips per country (top 10). Save the per-country tip count as `results/sequences_per_country.csv`.

2. **Figure 1 — country-fill world map.** Compute the dominant Nextstrain clade per country (mode), join to the world polygons from `ggplot2::map_data("world")`, fill each country polygon by clade. Pair it with the tree (clade-colored tips) using `patchwork`:

    ```r
    (p_tree | p_map) + plot_layout(widths = c(1, 1.3))
    ```

3. **Figure 2 — `phytools::phylo.to.map`.** Project the tree onto the world map by connecting each tip to its sampling coordinates. Save with sensible width/height — these figures benefit from being wide:

    ```r
    library(phytools)
    obj <- phylo.to.map(tree, coords, plot = FALSE)
    plot(obj, direction = "rightwards", colors = "darkblue",
         cex.points = 0.4, ftype = "off", lwd = 0.25)
    ```

4. **A 200-word interpretation.** Look at both figures. Pick **one specific apparent geographic claim** that one of them invites a reader to make ("clade X is concentrated in Europe", "sampling is uniform across the globe", etc.) and explain *one* way that claim could be wrong because of **sampling bias**. Use `results/sequences_per_country.csv` as evidence.

5. **(Optional extension) The bias-correction figure.** Join `covid_cases_per_country.csv` to your sequences-per-country table, compute `sequences_per_million_cases = sequences / total_cases * 1e6`, and plot it as a bar chart or as a second world-map fill. Countries that look "important" on the original map may look much less so once normalized — and vice versa.

## Pitfalls already burned (carry these forward)

- **`map_data("world")` uses different country names than Nextstrain.** Notable mismatches you will hit: `USA` → `United States` (US) → various; `UK` → `Great Britain` and `Ireland` (separated). When in doubt, `unique(world$region)` lists every name `map_data` knows about.
- **`phytools::phylo.to.map` requires coordinates as a 2-column matrix with `tree$tip.label` as row names**, in the order `c("latitude", "longitude")`. Subsetting tips and reordering the matrix is essential — pass it the wrong way and the map is silently wrong.
- **Tips with `NA` coordinates must be dropped.** `phylo.to.map` errors on missing values. Use `ape::keep.tip(tree, common)` where `common <- intersect(tips_with_coords, tree$tip.label)`.
- **Country of *sampling* ≠ country of *infection origin*.** Travel-related cases are sequenced in the destination country; sampling clusters appear in countries with strong genomic surveillance, not necessarily in countries where the variant is most prevalent.

## Why this is a good first project

The two encodings here introduce ggplot2-on-real-spatial-data, the `phytools` ecosystem (which extends `ape` with comparative methods), and the data wrangling needed to align country-name conventions across three independent data sources. None of those are unique to SARS-CoV-2 — once you can do this for ncov, you can do it for any phylogeny with geographic metadata.

## Going further

- Reproduce the dominant-clade fill, but with **earliest sampling date** instead of clade — does it tell a different story about variant emergence?
- Use `ggrepel` to label countries with the most sequences directly on the map.
- Try the spatial-tree composite — plot tree edges as great-circle arcs over the map.

## References

- Revell, L. J. (2012). phytools. *Methods Ecol. Evol.* 3:217.
- Yu G. (2022). *Data Integration, Manipulation and Visualization of Phylogenetic Trees.* <https://yulab-smu.top/treedata-book/>.
- Mathieu et al. (2020). *Coronavirus Pandemic (COVID-19).* OurWorldInData.org. The case-count snapshot in `data/` is sourced from this dataset.
