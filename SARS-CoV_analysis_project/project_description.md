# SARS-CoV-2 Phylogeny with Rich Metadata Visualization

## Background

Between December 2019 and the time of writing, SARS-CoV-2 became the most densely sampled organism in the history of biology. Public sequence repositories — GISAID, ENA, NCBI Virus — hold millions of complete genomes annotated with sampling date, geography, host, and (since 2021) Pango lineage and Nextstrain clade. The Nextstrain project distills this firehose into time-calibrated phylogenies updated several times per week.

A phylogenetic tree on this scale is fundamentally a **2-D coordinate system**: the y-axis indexes tips, the x-axis carries either substitutions per site or calendar time. Once a tree is drawn, **any metadata table keyed by tip label can be bound to that coordinate system**. The resulting figures make biology that is invisible in a spreadsheet — convergent evolution of receptor-binding mutations, geographic spread of lineages, host jumps — *visually obvious*.

The R ecosystem around `ggtree` provides six principled visual encodings for this binding:

| Encoding             | ggtree mechanism                              | Best for                                           |
|----------------------|-----------------------------------------------|----------------------------------------------------|
| Tip aesthetics       | `geom_tippoint(aes(...))`                     | 1–2 categorical/continuous variables               |
| Branch aesthetics    | `aes(color = ...)` in `ggtree()`              | Inferred ancestral state along the tree            |
| Adjacent strips      | `geom_facet()` / `ggtreeExtra::geom_fruit()`  | Many variables, mixed types                        |
| Heatmap matrix       | `gheatmap()`                                  | Wide categorical matrices (mutations, presence)    |
| MSA alignment        | `msaplot()`                                   | Linking sequence sites to tree topology            |
| Geographic projection| `phytools::phylo.to.map`, ggplot2 maps        | Spatial spread                                     |

This project walks through all six on the same SARS-CoV-2 dataset so that the student sees the same biology rendered six different ways and learns to *choose* the encoding that matches the question being asked.

## Problem statement

Given a time-calibrated SARS-CoV-2 phylogeny of ≈300–500 tips and a Nextstrain-style metadata table (date, region, country, clade, Pango lineage, key spike mutations, host, age, sex, latitude/longitude), produce a series of figures that:

1. place each tip in calendar time;
2. encode clade or lineage on tips and on branches;
3. show the **spike mutation profile** as a heatmap aligned to the tree;
4. show heterogeneous metadata (categorical + continuous) as adjacent strip panels;
5. show a multiple-sequence-alignment of the spike receptor-binding domain (RBD) beside the tree;
6. project tips onto a world map.

Each figure is paired with a written critical reflection on what biological conclusion the figure invites the viewer to draw, and on **one** way that conclusion could be wrong because of sampling, the inference method, or the visual encoding.

## Working hypothesis

> The Nextstrain-style annotations of a SARS-CoV-2 phylogeny — clade, lineage, key spike mutations, region — when rendered as visual layers on a time-calibrated tree, are sufficient to *see* the major epidemiological events of 2020–2023: the emergence of D614G, the rise of Alpha, the convergent appearance of N501Y in Alpha and Omicron, and the geographic stratification of sampling effort. The same figures will simultaneously make several **artifacts of sampling** visually obvious to a critical reader.

## Data sources

Two pre-cached files distributed with the project, plus one optional alignment:

- `data/ncov_subsample.nwk` — Nextstrain-derived **time-calibrated** Newick tree (branches in years), ≈300–500 tips, regional or globally subsampled.
- `data/ncov_metadata.tsv` — tab-separated metadata, one row per tip, with at least: `strain`, `date`, `region`, `country`, `division`, `Nextstrain_clade`, `pango_lineage`, `host`, `age`, `sex`, `length`, `latitude`, `longitude`, `nextclade_mutations`.
- `data/spike_RBD_aligned.fasta` — multiple sequence alignment of the spike RBD (≈230 nt window), tip labels matching `tree$tip.label`. Optional; used in Step 7 only.

Acquisition, format details, and the join-key invariant `setdiff(tree$tip.label, meta$strain)` are documented in [`nextstrain_howto.md`](./nextstrain_howto.md).

## Research questions

1. **Time-axis tree.** When tips are placed on a calendar x-axis, what is the temporal envelope of each Nextstrain clade? Which clades co-circulate, and which displace one another?
2. **Tip vs. branch encoding.** A clade color on tips says "this is where the lineage ended up." A clade color on branches says "this is the inferred evolutionary history." Where do the two encodings disagree, and why?
3. **Mutation matrix.** Which spike mutations recur on multiple branches of the tree (convergent evolution)? N501Y, E484K, and L452R are textbook examples — does the figure reproduce them?
4. **Multivariate strip plot.** When region, host age, and genome length are bound to the same tree, do the strips reveal correlations between metadata columns that are independent of the tree itself?
5. **MSA window.** When the RBD alignment is rendered beside the tree, do mutation columns line up with branches that bear the corresponding `S:...` mutation in the heatmap? If not, why?
6. **Geography.** Country-level fill on a world map is seductive. After computing **sequences-per-country** vs. **cases-per-country**, which countries are over- or under-represented in the sample, and how does that bias every other figure?

## Scope

A clean comparative analysis on a teaching dataset requires deliberate restraint:

- **Tip count**: 300–500. Big enough to be biologically interesting, small enough to render in seconds.
- **Time span**: roughly 2020-01 through the most recent sampling date in the file (the *mrsd*).
- **Branch units**: time (years). The `mrsd` mechanism in `ggtree` only works when branches are in time. Verify with `max(node.depth.edgelength(tree))` — should be ≈3–5 for a multi-year SARS-CoV-2 tree. If the value is around 0.001 the tree is in substitutions/site and must be clock-calibrated first.
- **Layout**: rectangular for ≤500 tips; switch to a fan layout via `ggtreeExtra::geom_fruit()` when the tip count grows beyond that.
- **No live API calls during analysis.** Data acquisition happens once, ahead of time, and is documented in `nextstrain_howto.md`. Analysis scripts read from `data/`.

## Outputs

The project produces:

1. **A reproducible pipeline** — R scripts (or a Quarto notebook) that, starting from the three files in `data/`, generate every figure end-to-end without manual intervention.
2. **A Quarto report** containing:
    - the four required figures: time-tree with tip-clade color, branch-colored tree (ancestral state reconstruction), `gheatmap` with ≥6 spike mutations, and a multi-strip facet plot;
    - one optional integrated final figure suitable for publication;
    - a 300-word **critical-thinking essay**: pick *one* of the figures, identify a biological conclusion it invites, and explain *one* way that conclusion could be wrong because of sampling bias, the inference method, or the visual encoding.
3. **Intermediate data products** as CSV/TSV under `results/` — clade × date counts, the spike-mutation presence/absence matrix, per-country sampling counts.
4. **A README** in the project root describing how to reproduce the analysis.

## Optional extensions

- **Interactive HTML.** `plotly::ggplotly(p1)` makes tip metadata hoverable.
- **Bayesian time-tree.** Run a small BEAST analysis on 50 sequences and load via `treeio::read.beast()`. `geom_range("height_0.95_HPD")` adds HPD violins on internal nodes.
- **Tip silhouettes.** For cross-host coronaviruses (bat, pangolin), `rphylopic::geom_phylopic()` swaps tip labels for organism icons.
- **Subclade zoom.** `tree_subset(tree, node = MRCA(tree, omicron_tips), levels_back = 0)` extracts a clade, and `viewClade(p, node)` does the same interactively without recomputing the layout.
- **Transferability note.** The same pipeline transfers without change to other densely sampled organisms — *Pisum sativum* population genomics, *Mycobacterium tuberculosis* outbreak trees, etc. Swap `Nextstrain_clade` for ecotype or strain, and `nextclade_mutations` for marker SNPs.

## Pre-warned pitfalls

- **Tip-label mismatch** after `dplyr` cleanup — always run `setdiff(tree$tip.label, meta$strain)` before joining.
- **Heatmap row order** — `gheatmap` aligns by **row name**, not row order. Don't sort the matrix manually.
- **Forgetting `new_scale_fill()` / `new_scale_color()`** — produces a single merged scale that maps two variables to the same colors.
- **Branch length units** — `mrsd` only works if branches are in time. Always check `max(node.depth.edgelength(tree))`.
- **`facet_widths()` must come last** — it operates on the assembled gtable.
- **Naive ancestral state reconstruction.** `ape::ace(..., model = "ER")` ignores sampling bias and assumes time-reversibility. SARS-CoV-2 sampling is wildly geographically biased (UK and USA dominate). For real epidemiology, use **structured coalescent** methods (BEAST MASCOT). Read Morel et al. 2021 *MBE* before trusting any branch-color plot.
- **Country of sampling ≠ country of infection origin.** Geographic visualizations are easily misleading. Compute **sequences-per-country relative to reported cases-per-country** before drawing geographic conclusions.

## Key references

- Yu G. (2022). *Data Integration, Manipulation and Visualization of Phylogenetic Trees.* <https://yulab-smu.top/treedata-book/> — chapters 7 (`msaplot`, `gheatmap`) and 10 (`geom_facet`) are the canonical reference.
- Xu S. et al. (2022). ggtreeExtra: compact visualization of richly annotated phylogenetic data. *MBE* 39:msab166.
- Hadfield et al. (2018). Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics* 34:4121.
- Morel et al. (2021). Phylogenetic Analysis of SARS-CoV-2 Data Is Difficult. *MBE* 38:1777. **Required reading before the critical-thinking essay.**
