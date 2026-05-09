# Geographic-projection starter — instructions for the coding agent

This is a single-question starter project. Solve **one** question well; do not expand scope.

## The question

> How can we project a SARS-CoV-2 phylogeny onto a world map?

The deliverable is a Quarto report containing two complementary figures (a country-fill choropleth and a `phytools::phylo.to.map` projection), each paired with a brief written interpretation, plus 200 words on **sampling bias**. See `project_description.md` for full framing.

## Project goals

- Render a per-country dominant-clade choropleth using `ggplot2::map_data("world")`.
- Render a tree-to-map projection using `phytools::phylo.to.map`.
- Document one concrete way that sampling bias undermines an apparent geographic conclusion drawn from the figures, using the per-country sequence count saved to `results/`.

## Expected layout

- `data/` — already populated, do not modify
- `scripts/` — your R scripts (create as you go)
- `results/` — your CSV outputs (create as you go)
- `reports/` — your Quarto report and figures (create as you go)

If a directory is missing and is needed, create it with that exact name.

## Data

Three files are already in `data/`:

- `ncov_subsample.nwk` — 400-tip Newick tree, branches in years.
- `ncov_metadata.tsv` — per-tip metadata, includes `country`, `latitude`, `longitude`, `Nextstrain_clade`.
- `covid_cases_per_country.csv` — 233 countries, snapshot from Our World in Data, used in the optional sampling-bias extension.

Read with `ape::read.tree()`, `readr::read_tsv()`, `readr::read_csv()`. Do **not** download anything; everything you need is on disk.

## R coding style

- Tidyverse for wrangling, `ggplot2` + `maps` + `patchwork` for the choropleth, `phytools` for `phylo.to.map`, Quarto for the report.
- Short, clear, teaching-friendly code. One script per logical step.
- Use relative paths from the project root: `data/...`, `results/...`, `reports/...`.

## Pitfalls to avoid (these are the ones that bite first-time map+tree users)

- **Country-name mismatches across three sources.** Nextstrain metadata uses `USA` and `Czech Republic`; OWID uses `United States` and `Czechia`; `map_data("world")` uses `USA` (matches metadata) and `Czech Republic` (matches metadata). Plan one tidy `case_when()` to align names *before* any join. `unique(world$region)` is the authoritative list of names that `map_data` knows about.
- **`phylo.to.map` matrix shape.** It needs a 2-column matrix, columns `c("latitude", "longitude")` *in that order*, with `tree$tip.label` as row names, every value non-`NA`. Drop tips with `NA` coordinates with `ape::keep.tip(tree, common)`.
- **Country of sampling vs. country of origin.** Always state in the report that "country" means "country of sampling," not "country where the variant emerged." Travel cases land where they are sequenced.
- **Join key invariant.** Always run `setdiff(tree$tip.label, meta$strain)` before binding metadata to the tree.
- **`map_data("world")` returns POLYGONS at long/lat granularity** — this is a long data frame, not a spatial object. Pipe it into `ggplot() + geom_polygon(aes(long, lat, group = group, fill = ...))`.

## Scope discipline

Do **not** also build the spike mutation heatmap, the ancestral state reconstruction, or the multi-strip facet plot from the full project. Those live in `../SARS-CoV_analysis_project/` and the sibling starters. Keep this project focused on geographic projection.

If you finish the core figures quickly, the optional **sampling-bias extension** in `project_description.md` is the right next step — not adding new encodings.
