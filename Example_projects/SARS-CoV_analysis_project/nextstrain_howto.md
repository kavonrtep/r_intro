# Acquiring SARS-CoV-2 phylogeny data — project-scoped cheat sheet

This project ships with the data files already cached in `data/`. **You do not need to re-download anything to run the analysis.** This document exists so that:

1. an instructor can refresh the snapshot for a new term, and
2. a student who wants to swap in a different subsample (regional build, host-jump build, etc.) knows where to look.

Three files are expected in `data/`:

| File                          | Format                       | Source                                  |
|-------------------------------|------------------------------|-----------------------------------------|
| `ncov_subsample.nwk`          | Newick, branches in **years**| Nextstrain `ncov` build                 |
| `ncov_metadata.tsv`           | TSV, one row per tip         | Nextstrain `ncov` build                 |
| `spike_RBD_aligned.fasta`     | FASTA, aligned, optional     | Nextstrain alignment, RBD window        |

---

## 1. Where Nextstrain publishes its builds

Nextstrain hosts its public SARS-CoV-2 builds at <https://nextstrain.org/ncov>. The builds are partitioned by region, time window, and subsampling strategy. Examples:

| Build URL fragment                       | What it is                                            |
|------------------------------------------|-------------------------------------------------------|
| `ncov/gisaid/global/all-time`            | Full global build, all-time, GISAID-derived           |
| `ncov/gisaid/global/6m`                  | Global, last 6 months                                 |
| `ncov/open/europe/all-time`              | Europe, all-time, open-data (no GISAID auth needed)   |
| `ncov/open/north-america/6m`             | North America, last 6 months, open-data               |

The **open** builds use ENA/NCBI/CNCB data and have no access restriction. Prefer them for a teaching project — students can refresh the data without a GISAID login. They are smaller (≈10–30k tips) but for our purposes a regional 6-month build subsamples to a teaching-friendly size out of the box.

For this project we want **300–500 tips**. Pick a regional build (e.g. `open/europe/6m`) or a globally subsampled build.

---

## 2. Downloading the tree + metadata pair

Each Nextstrain build exposes its source files at predictable URLs under `data.nextstrain.org`. The pair we need:

```
https://data.nextstrain.org/files/ncov/open/<region>/<window>/<dataset>.json
https://data.nextstrain.org/files/ncov/open/<region>/<window>/metadata.tsv.zst
```

The `.json` is an **Auspice v2 dataset** — a single file containing the time tree, internal-node annotations, and per-tip attributes. It is not a Newick file directly.

### Two paths to a Newick + TSV pair

**Path A — extract from the Auspice JSON** (no Nextstrain pipeline needed):

```r
# install.packages("jsonlite"); install.packages("ape")
library(jsonlite); library(ape); library(dplyr); library(readr)

aus <- read_json("ncov_open_europe_6m.json")

# Auspice v2 puts the tree under aus$tree (a nested node tree).
# treeio's read.beast.newick is heavy-handed for this; the simplest
# extractor is the helper in the `auspice2treedata` snippet on the
# ggtree mailing list, or — pragmatically — the `nextstrain-cli` tool:
```

**Path B — use `nextstrain-cli` and the `augur` pipeline**:

```bash
# pip install nextstrain-cli  (or via mamba)
nextstrain remote download \
  https://nextstrain.org/ncov/open/europe/6m \
  ./ncov_open_europe_6m

# This pulls the Auspice JSON + the upstream metadata.tsv.zst alongside it.
zstd -d metadata.tsv.zst   # produces metadata.tsv
```

For the Newick: the cleanest route is to run `augur tree` + `augur refine` on the source alignment, but that is heavy. For teaching, the practical alternative is to **export the time-tree from the Auspice JSON via `auspice_tree_to_table.py`** (a small upstream helper) or via the R `auspicer` / `treeio::read.nextstrain.json()` wrappers if present in your version.

The instructor-cached file `data/ncov_subsample.nwk` was produced by one of these paths.

---

## 3. Subsampling to 300–500 tips

If the build is bigger than 500 tips, subsample with `ape`:

```r
library(ape)
tree <- read.tree("data/ncov_full.nwk")
set.seed(42)
keep <- sample(tree$tip.label, 400)
tree_sub <- keep.tip(tree, keep)
write.tree(tree_sub, "data/ncov_subsample.nwk")
```

Then prune the metadata accordingly:

```r
library(readr); library(dplyr)
meta <- read_tsv("data/metadata_full.tsv")
meta_sub <- meta %>% filter(strain %in% tree_sub$tip.label)
write_tsv(meta_sub, "data/ncov_metadata.tsv")
```

**Stratified subsampling** is better than random for teaching — sample N tips per Nextstrain clade so that all major variants (Alpha, Delta, BA.1, BA.2, BA.5, XBB, etc.) are represented:

```r
target_per_clade <- 30
keep <- meta %>%
  group_by(Nextstrain_clade) %>%
  slice_sample(n = target_per_clade) %>%
  pull(strain)
```

---

## 4. Required metadata columns

The analysis scripts expect `data/ncov_metadata.tsv` to have at least:

| Column                | Type        | Notes                                          |
|-----------------------|-------------|------------------------------------------------|
| `strain`              | character   | **Join key.** Must match `tree$tip.label`      |
| `date`                | YYYY-MM-DD  | Sampling date                                  |
| `region`              | character   | Continental region                             |
| `country`             | character   | Country                                        |
| `division`            | character   | Sub-national division (state/province)         |
| `Nextstrain_clade`    | character   | e.g. `21K (Omicron)`, `20I (Alpha)`            |
| `pango_lineage`       | character   | e.g. `BA.1`, `B.1.1.7`                         |
| `host`                | character   | Usually `Homo sapiens`; cross-host builds vary |
| `age`                 | numeric     | May be missing                                 |
| `sex`                 | character   | May be missing                                 |
| `length`              | integer     | Genome length in nt                            |
| `latitude`            | numeric     | For geographic projection                      |
| `longitude`           | numeric     | For geographic projection                      |
| `nextclade_mutations` | character   | Comma-separated, e.g. `S:N501Y,S:D614G,...`    |

If any of these are missing in a refreshed build, either pad the column with `NA` or skip the figure that needs it (and document the skip in the report).

---

## 5. The join-key invariant

Always run **before** joining tree to metadata:

```r
setdiff(tree$tip.label, meta$strain)   # MUST be empty
setdiff(meta$strain, tree$tip.label)   # tips not in tree (often non-empty, fine)
```

Nextstrain strain names look like `Australia/NSW-1234/2021` or `hCoV-19/USA/CA-CDC-...|2022-03-15`. Tree-building tools sometimes:

- Replace `|` with `_`
- Replace `/` with `_`
- Truncate at the first whitespace

If the first `setdiff` is non-empty, normalize **both** sides with the same regex before joining. Don't fix one side only.

---

## 6. Branch-length units

```r
max(ape::node.depth.edgelength(tree))
```

- ≈ 3–5 → branches are in **years** (time-calibrated). Use `mrsd` directly.
- ≈ 0.001–0.01 → branches are in **substitutions/site**. Time-calibrate first via `treeio::read.beast()` for BEAST output, or `ape::chronos()` as a quick approximation.

The cached `data/ncov_subsample.nwk` is in years.

---

## 7. The optional spike RBD alignment

`data/spike_RBD_aligned.fasta` is a multiple sequence alignment of the spike RBD, ≈ 230 nt window, with sequence headers matching `tree$tip.label` exactly.

To produce it from a full-genome alignment:

```r
library(Biostrings)
genome_aln <- readDNAStringSet("data/full_alignment.fasta")
rbd <- subseq(genome_aln, start = 22517, end = 23185)   # RBD nt range
rbd <- rbd[names(rbd) %in% tree$tip.label]
writeXStringSet(rbd, "data/spike_RBD_aligned.fasta")
```

The RBD nt range varies by reference. Check against Wuhan-Hu-1 (NC_045512.2) before refreshing.

---

## 8. The shipped extraction script

The cached files in `data/` were produced by `data/build_data.R`, which:

1. Downloads `https://data.nextstrain.org/ncov_open_global_6m.json` (≈7 MB after `gzip` decoding) into `data/_auspice_cache.json`.
2. Walks the Auspice v2 JSON, builds an `ape::phylo` object with branches in **years** (derived from `num_date` differences between parent and child).
3. Accumulates spike (`S`) amino-acid mutations along each root → tip path, handling reversions, into the `nextclade_mutations` column.
4. Pulls per-country `latitude`/`longitude` from `aus$meta$geo_resolutions`.
5. Subsamples to 400 tips stratified by `Nextstrain_clade` (~10 per clade by default).
6. Writes `data/ncov_subsample.nwk` and `data/ncov_metadata.tsv`.

Run it from the project root:

```bash
Rscript data/build_data.R
```

The 6.6 MB `_auspice_cache.json` is gitignore-friendly: re-running `build_data.R` will fetch it again if absent. The student-facing files are the 16 KB Newick and the 200 KB TSV.

## 9. Refresh checklist

When refreshing `data/` for a new term:

1. Pick a Nextstrain build — open-data, ≤500 tips after subsampling. Edit the `dl_url` line in `data/build_data.R` if you want a regional build (e.g. `ncov_open_europe_6m.json`).
2. Delete `data/_auspice_cache.json` to force re-download.
3. `Rscript data/build_data.R` — the script asserts the join-key invariant and the branch-unit sanity check (`max(node.depth.edgelength) > 0.1`).
4. Inspect the new outputs (clade distribution should still be diverse, ~10 per clade).
5. (Optional) Slice the spike RBD window from a reference alignment → `data/spike_RBD_aligned.fasta`.
6. Update the date stamp below.

Last refreshed: 2026-05-09.

---

## 9. References

- Nextstrain ncov build: <https://nextstrain.org/ncov>
- Nextstrain CLI: <https://docs.nextstrain.org/projects/cli/>
- Augur (the Nextstrain pipeline): <https://docs.nextstrain.org/projects/augur/>
- ggtree treedata book, chapter on Auspice import: <https://yulab-smu.top/treedata-book/>
