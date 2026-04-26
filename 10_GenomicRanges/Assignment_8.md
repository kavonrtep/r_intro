# Assignment 8: Annotating ChIP-seq Peaks in the Zebrafish Genome

## Background: Where does a transcription factor bind?

In a ChIP-seq experiment we sequence DNA fragments that were pulled down by an
antibody against a specific transcription factor (TF). After mapping and peak
calling we end up with a list of **genomic intervals** — the **peaks** — that
mark candidate binding sites. The biology only becomes interesting once we
**annotate** those peaks: which gene is each peak closest to? Does the peak
sit in a promoter, inside the gene body, or in an intergenic desert? Where
exactly relative to the **transcription start site (TSS)** does the TF
prefer to bind?

In this assignment you will analyse a (fictional) set of **~2000 ChIP-seq
peaks** distributed across **all 25 chromosomes** of the zebrafish genome
(*Danio rerio*, danRer11). The peaks were generated to mimic the binding
pattern of a typical promoter-binding TF — most peaks should cluster a few
hundred bp upstream of TSSs, with a smaller number sitting far away in
intergenic regions. Your job is to recover that pattern from the data.

You will use:

- the full **zebrafish RefSeq gene annotation** (UCSC RefSeq, all standard
  chromosomes), and
- the **danRer11 reference genome** from `BSgenome.Drerio.UCSC.danRer11`.

---

## Overview

In this assignment you will:

- Import the gene annotation (GFF3) and ChIP-seq peaks (BED) into GRanges.
- Build a **promoter** GRanges using `flank()` / `promoters()`.
- Classify each peak as **promoter / genic / intergenic**.
- Find the **nearest gene** to every peak and the **distance** to it.
- Compute the **signed peak-to-TSS offset** (strand-aware) and plot its
  distribution — this is the headline figure of the assignment.
- Use `subsetByOverlaps()` to keep only peaks near genes.
- Visualise individual peaks with the `plotGenomicTracks()` helper from
  class.
- Extract peak sequences from `BSgenome` and compute their GC content.

Total expected time: **45–60 min**.

---

## Setup

Make sure these packages are installed:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges",
  "rtracklayer",
  "Biostrings",
  "BSgenome.Drerio.UCSC.danRer11",  # zebrafish genome (~450 MB) — already used in lesson 9
  "Gviz"
))
```

Load libraries at the top of your script:

```r
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Drerio.UCSC.danRer11)
library(Gviz)
```

Paste in the `plotGenomicTracks()` helper from the lesson — you'll use it
in Tasks 2 and 7:

```r
plotGenomicTracks <- function(tracks_list,
                              colors     = NULL,
                              chromosome = NULL,
                              from       = NULL,
                              to         = NULL) {
  n <- length(tracks_list)
  if (is.null(names(tracks_list)) || any(names(tracks_list) == "")) {
    stop("tracks_list must be a named list (names used as track titles).")
  }
  if (!is.null(colors) && length(colors) != n) {
    stop("If colors is provided, its length must equal number of tracks.")
  }
  if (is.null(colors)) colors <- rep("gray", n)

  annotation_tracks <- mapply(function(gr, col, nm) {
    args <- list(range = gr, name = nm, fill = col,
                 col = "black", stacking = "dense")
    if ("id" %in% names(mcols(gr))) {
      ids <- mcols(gr)$id
      if (!all(is.na(ids)) && length(ids) == length(gr)) {
        args$showFeatureId     <- TRUE
        args$featureAnnotation <- "id"
        args$fontcolor.feature <- "black"
      }
    }
    do.call(AnnotationTrack, args)
  },
  gr = tracks_list, col = colors, nm = names(tracks_list),
  SIMPLIFY = FALSE)

  axis_track <- GenomeAxisTrack()
  all_tracks <- c(list(axis_track), annotation_tracks)
  plotTracks(all_tracks,
             chromosome       = chromosome,
             from             = from,
             to               = to,
             background.panel = "#FFFEDB",
             col.axis         = "black",
             cex.title        = 0.8,
             extend.left      = 5,
             extend.right     = 5)
  invisible(annotation_tracks)
}
```

Call it with a *named* list of GRanges — names become track titles. Set
`chromosome`, `from`, and `to` to control the plotted window.

The data files for this assignment live in the same folder:

| File | Content |
|---|---|
| `data/zebrafish_genes.gff3` | All RefSeq genes on chr1–25 (~13,000 genes) |
| `data/zebrafish_peaks.bed`  | ~2,000 fictional ChIP-seq peaks across the genome |

---

## Task 1: Load the Data and Inspect It

Read the gene annotation and the peak file into `GRanges` objects.

- Load both files into GRanges and report how many genes and peaks you have.
- How many peaks per chromosome? Is the distribution roughly proportional to
  chromosome size?
- What metadata columns does each object have?

<details>
<summary>💡 Hint</summary>

Look at `rtracklayer::import()`. To inspect a GRanges, you may need
`length()`, `seqnames()`, `table()`, `mcols()`.

</details>

---

## Task 2: Build a Promoter GRanges

Define a promoter as the region **2 kb upstream of the TSS** of each gene
(strand-aware).

- Build a GRanges of promoters from the gene annotation.
- How many promoter ranges do you get? What is their median width?
- Pick any gene (e.g. the first one on chr1) and plot the gene together with
  its promoter using `plotGenomicTracks()`.

<details>
<summary>💡 Hint</summary>

There is a one-line helper in `GenomicRanges` that makes a strand-aware
promoter GRanges directly — search the docs of `?promoters` (or look at
`flank()` if you prefer to do it manually). For the visualization, you'll
need to set `chromosome`, `from`, and `to` in `plotGenomicTracks()`.

</details>

---

## Task 3: Classify Each Peak as Promoter / Genic / Intergenic

Assign one of three labels to every peak:

- **`promoter`** — overlaps a 2 kb-upstream promoter,
- **`genic`**    — does not overlap a promoter but overlaps a gene body,
- **`intergenic`** — overlaps neither.

Add the label as a new metadata column `peaks$category` and report the count
in each category. Roughly what fraction of peaks land in promoters?

<details>
<summary>💡 Hint</summary>

Have a look at `overlapsAny()`. You will then need a way to combine two
logical vectors into a three-level label — base R `ifelse()` chained twice
does the job. To count categories, use `table()`.

The order of the test matters when a peak overlaps both a promoter and a
gene body: which label should win?

</details>

---

## Task 4: Annotate Each Peak with its Nearest Gene

For every peak, find the **nearest gene** and the **distance** to it.

- Add two new columns to `peaks`: `nearest_gene` (the gene Name) and
  `dist_to_gene` (integer, in bp; 0 = peak inside a gene).
- Print the **5 closest** peak–gene pairs (by distance, but skip peaks that
  literally overlap a gene — those have distance 0) and the **5 most
  distant**.

<details>
<summary>💡 Hint</summary>

`nearest()` returns indices, `distance()` returns bp distances (0 = overlap).
To assign metadata, use `mcols(peaks)$col <- ...` or `peaks$col <- ...`.
Use `order()` on `peaks$dist_to_gene` for the sorted slice, and a logical
filter to drop the overlapping peaks before sorting "closest non-overlapping".

</details>

---

## Task 5: Peak Position Relative to the TSS — the Headline Figure

For every peak, compute its **signed offset** from the TSS of its nearest
gene:

- `offset = peak_center - TSS` for genes on the **+** strand,
- `offset = TSS - peak_center` for genes on the **−** strand.

So **negative offsets are upstream** of the TSS and **positive offsets are
downstream** (into the gene body), regardless of which strand the gene is
on.

- Add the column `peaks$tss_offset` (integer bp).
- **Plot a histogram** of `tss_offset`, restricted to **±5 kb** so the
  promoter-proximal peak is visible. Add a vertical dashed line at 0 (the
  TSS).
- Save the histogram to `work_dir/task5_tss_offset_histogram.png`.
- **Interpret**: at roughly what offset does the histogram peak? What does
  that suggest about where this fictional TF prefers to bind?

<details>
<summary>💡 Hint</summary>

The TSS of a gene is its `start()` if it is on the `+` strand and its
`end()` if it is on the `-` strand — `ifelse()` on `strand(...)` is the
cleanest way to express this. The peak center is `(start + end) / 2`.

For the histogram use `hist()` with `xlim = c(-5000, 5000)` and a sensible
`breaks` (e.g. 100). Use `abline(v = 0, lty = 2)` to mark the TSS.

</details>

---

## Task 6: Subset to "Gene-Associated" Peaks

Many downstream analyses (e.g., motif discovery) only care about peaks that
sit close to a gene. Build a single GRanges containing **only the peaks that
fall in a promoter OR inside a gene body**.

- How many peaks remain? How does this match what you found in Task 3?

<details>
<summary>💡 Hint</summary>

`subsetByOverlaps()` works against any reference set. Two GRanges objects
can be combined with `c()`.

</details>

---

## Task 7: Visualise One Promoter Peak in Context

Pick one peak that was classified as `promoter` and draw a 3-track plot
showing the gene, the 2 kb promoter, and the peak itself in the surrounding
30 kb window.

- Use the `plotGenomicTracks()` helper.
- Set `from`/`to` to give 15 kb of padding on each side so you can see the
  gene structure.
- Save the figure to `work_dir/task7_promoter_peak_in_context.png`.

<details>
<summary>💡 Hint</summary>

Subset the peaks by `peaks$category == "promoter"` and pick one. Use
`nearest()` to look up its associated gene. The `from` / `to` arguments of
`plotGenomicTracks()` accept absolute coordinates. Don't forget to set
`chromosome` to whatever chromosome that peak is on.

</details>

---

## Task 8: Extract Peak Sequences and Compute GC Content

Pull the DNA sequence under each peak from the zebrafish genome and compute
its GC content.

- Add the result as `peaks$gc_content`.
- **Compare**: is mean GC content **higher** in promoter peaks than in
  intergenic peaks? Summarise per category.

<details>
<summary>💡 Hint</summary>

`getSeq()` from `BSgenome` extracts a `DNAStringSet`. `letterFrequency()`
with `as.prob = TRUE` gives you the proportion of each letter. To summarise
a numeric column by a grouping factor, use `tapply()` or `aggregate()`.

Real promoter regions in vertebrate genomes are typically GC-rich (CpG
islands), so on **real** ChIP-seq data you would expect promoter peaks to
have higher GC. In this **fictional** dataset the peaks were placed without
biasing for CpG islands, so the three categories may end up similar — that
is OK. Report what you observe and explain whether the expected biological
pattern shows up.

</details>

---

## Task 9 (Bonus): Find Peaks That Are Both Promoter-Proximal and Highly Scored

The BED file ships with a `score` column (peak signal — higher = stronger
binding). Identify the "best" candidate peaks: peaks that

- overlap a promoter, **and**
- have `score` above the median peak score.

How many such peaks are there? Print the **top 20** by score, showing
`name`, `nearest_gene`, `tss_offset`, and `score`.

<details>
<summary>💡 Hint</summary>

You already have `peaks$category` from Task 3, `peaks$nearest_gene` from
Task 4, and `peaks$tss_offset` from Task 5 — combine them with a logical
condition on `peaks$score`. `median()` and `order()` are the only new
functions you need.

In a real ChIP-seq study, this list would be the input for downstream motif
analysis or pathway enrichment.

</details>

---

## What to submit

Submit a single R script named
- **`lastname_firstname_genomicranges_assignment.R`** that:
- Figures as PNGs
