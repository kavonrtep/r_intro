# Assignment 7: Differential Gene Expression Analysis with DESeq2

## Background: The `airway` Dataset

In this assignment you will analyse the **`airway`** dataset — the standard teaching dataset for RNA-seq differential expression in R/Bioconductor.

The data come from the following publication:

> Himes et al. (2014) *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.* **PLoS ONE** 9(6): e99625. doi:10.1371/journal.pone.0099625

The study investigated the transcriptomic response of primary human **airway smooth muscle (ASM)** cells to **dexamethasone** — a synthetic glucocorticoid (GC) widely used to treat asthma and other inflammatory diseases. ASM cells were treated with 1 µM dexamethasone (or left untreated as control) for 18 hours, and RNA-Seq was performed on **four independent ASM cell lines** (each from a different donor). This gives **8 samples** in total: 4 treated (`trt`) and 4 untreated (`untrt`).

The study identified **316 differentially expressed genes**, including the previously uncharacterised gene *CRISPLD2*, which was found to have SNPs associated with two asthma pharmacogenetic phenotypes. Several well-known GC-responsive genes — *DUSP1*, *KLF15*, *PER1*, *TSC22D3* — were also confirmed.

The `airway` Bioconductor package contains these data as a `RangedSummarizedExperiment` object with **63,677 genes × 8 samples**.

---

## Overview

In this assignment you will:

- Load the `airway` dataset and extract the count matrix and sample metadata
- Build a `DESeqDataSet` with a design that accounts for both cell line and treatment
- Run the DESeq2 differential expression pipeline
- Explore the data with a PCA plot
- Visualise results with an MA plot and a volcano plot
- Create an annotated heatmap of the top differentially expressed genes
- Annotate results with human gene symbols using `org.Hs.eg.db`
- Export the annotated results table

---

## Setup

Install the required packages if needed:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("airway", "DESeq2", "org.Hs.eg.db", "AnnotationDbi"))
install.packages(c("ggplot2", "pheatmap"))
```

Load libraries at the top of your script:

```r
library(airway)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
```

---

## Task 1: Load and Inspect the Data

Load the `airway` dataset and extract the count matrix and sample annotation.

- Use `data(airway)` to load the `SummarizedExperiment` object.
- Extract the count matrix with `assay()` and the sample metadata with `colData()`.
- Inspect the dimensions of the count matrix and print the sample metadata table.
- Verify that the column names of the count matrix match the row names of the metadata.

<details>
<summary>💡 Hint</summary>

```r
# Load the dataset
data(airway)

# Extract count matrix (genes × samples) and sample metadata
countData <- assay(airway)
colData   <- as.data.frame(colData(airway))

# Inspect
dim(countData)
head(countData[, 1:4])
print(colData)

# Verify alignment
all(rownames(colData) == colnames(countData))
```

The `colData` has columns `cell` (cell line ID) and `dex` (treatment: `"trt"` or `"untrt"`).

</details>

---

## Task 2: Create the DESeqDataSet and Filter Low-Count Genes

Create a `DESeqDataSet` and filter out genes with very low counts.

- Build a `DESeqDataSet` using `DESeqDataSetFromMatrix()`.
- Use **`design = ~ cell + dex`** — this accounts for variation between the four cell lines (donors) before testing for the treatment effect.
- Set `"untrt"` as the reference level for `dex`.
- Remove genes whose total count across all samples is less than 10.
- Report how many genes remain after filtering.

<details>
<summary>💡 Hint</summary>

```r
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = colData,
                              design    = ~ cell + dex)

# Set reference level
dds$dex <- relevel(dds$dex, ref = "untrt")

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")
```

Why `~ cell + dex`? Including `cell` as a blocking factor removes donor-to-donor variation, leaving a cleaner estimate of the treatment effect — similar to how `~ type + condition` was used in the pasilla analysis in class.

</details>

---

## Task 3: Run the DESeq2 Pipeline and Extract Results

Run the full DESeq2 analysis and extract results at two FDR thresholds.

- Run `DESeq()` on the filtered dataset.
- Extract results at the default threshold (FDR < 0.1) and print a summary.
- Extract results at a stricter threshold (FDR < 0.05) and print a summary.
- How many genes are significantly differentially expressed at FDR < 0.05?

<details>
<summary>💡 Hint</summary>

```r
# Run DESeq2 pipeline
dds <- DESeq(dds)

# Default results (FDR < 0.1)
res    <- results(dds)
summary(res)

# Stricter threshold
res05  <- results(dds, alpha = 0.05)
summary(res05)

# Count significant genes
sum(res05$padj < 0.05, na.rm = TRUE)
```

`results()` tests the last variable in the design by default — here `dex`, which is what we want.

</details>

---

## Task 4: Variance-Stabilizing Transformation and PCA

Apply VST and explore sample relationships with a PCA plot.

- Apply `vst(dds, blind = FALSE)` to obtain stabilised expression values.
- Use `plotPCA()` with `intgroup = c("cell", "dex")` to get PCA coordinates.
- Create a `ggplot2` scatter plot of PC1 vs PC2.
- Color points by `dex` (treatment) and use different shapes for `cell` (cell line).
- Add variance-explained percentages to the axis labels.
- **Interpret:** Which factor — treatment or cell line — drives the most separation along PC1?

<details>
<summary>💡 Hint</summary>

```r
# VST transformation
vsd <- vst(dds, blind = FALSE)

# PCA data
pcaData   <- plotPCA(vsd, intgroup = c("cell", "dex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot
ggplot(pcaData, aes(PC1, PC2, color = dex, shape = cell)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of airway samples") +
  theme_minimal()
```

</details>

---

## Task 5: MA Plot and Volcano Plot

Visualise the differential expression results.

**MA plot:**
- Use `plotMA(res05, ylim = c(-5, 5))` to show fold change versus mean expression.
- Blue points are significant genes. What pattern do you observe?

**Volcano plot:**
- Convert `res05` to a data frame and add a gene ID column.
- Create a `ggplot2` scatter plot of `log2FoldChange` (x-axis) vs `-log10(padj)` (y-axis).
- Add a color column: `"up"` if log2FC > 1 and padj < 0.05, `"down"` if log2FC < -1 and padj < 0.05, `"ns"` otherwise.
- Use `scale_color_manual()` with meaningful colors (e.g., red/blue/grey).
- Add vertical lines at ±1 and a horizontal line at -log10(0.05).

<details>
<summary>💡 Hint</summary>

```r
# MA plot
plotMA(res05, main = "airway: DEX vs untreated (FDR < 0.05)", ylim = c(-5, 5))

# Volcano plot
volc <- as.data.frame(res05)
volc$gene <- rownames(volc)
volc$group <- "ns"
volc$group[volc$log2FoldChange >  1 & !is.na(volc$padj) & volc$padj < 0.05] <- "up"
volc$group[volc$log2FoldChange < -1 & !is.na(volc$padj) & volc$padj < 0.05] <- "down"

ggplot(volc, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c(up = "red", down = "steelblue", ns = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggtitle("Volcano Plot: DEX vs untreated") +
  xlab("Log2 Fold Change") + ylab("-Log10 Adjusted P-value") +
  theme_minimal()
```

</details>

---

## Task 6: Heatmap of Top Differentially Expressed Genes

Create an annotated heatmap showing the top 30 differentially expressed genes.

- Select the top 30 genes by absolute log2 fold change among those with padj < 0.05.
- Extract their VST-normalised expression values from the `vsd` object using `assay()`.
- Create sample annotation for the heatmap columns showing `dex` and `cell`.
- Use `pheatmap()` with `scale = "row"` and a blue-white-red colour palette.
- Set `show_rownames = TRUE` so gene IDs are visible.

<details>
<summary>💡 Hint</summary>

```r
# Select top 30 DE genes by |log2FC|
sig_genes <- subset(as.data.frame(res05), padj < 0.05 & !is.na(padj))
top30     <- head(sig_genes[order(-abs(sig_genes$log2FoldChange)), ], 30)

# Extract VST values for these genes
mat <- assay(vsd)[rownames(top30), ]

# Sample annotation
anno_col <- data.frame(treatment = colData$dex,
                       cell_line = colData$cell,
                       row.names = rownames(colData))

# Heatmap
pheatmap(mat,
         scale           = "row",
         annotation_col  = anno_col,
         show_rownames   = TRUE,
         cluster_cols    = TRUE,
         cluster_rows    = TRUE,
         color           = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
         main            = "Top 30 DE genes: DEX vs untreated")
```

</details>

---

## Task 7: Gene Annotation

Add human gene symbols and descriptions to the results using `org.Hs.eg.db`.

- The row names of the results are Ensembl gene IDs (e.g., `ENSG00000000003`).
- Use `mapIds()` to retrieve gene symbols (`SYMBOL`) and gene names (`GENENAME`) from `org.Hs.eg.db`.
- Add these as new columns to the results data frame.
- Sort by adjusted p-value and export to `work_dir/airway_deseq2_results.csv`.
- Print the top 10 genes (by padj) with their symbol and fold change.

<details>
<summary>💡 Hint</summary>

```r
# Convert results to data frame
res_df       <- as.data.frame(res05)
res_df$ensembl <- rownames(res_df)

# Map Ensembl IDs to gene symbols and names
# keytype = "ENSEMBL" because row names are Ensembl IDs
res_df$symbol   <- mapIds(org.Hs.eg.db,
                           keys      = res_df$ensembl,
                           column    = "SYMBOL",
                           keytype   = "ENSEMBL",
                           multiVals = "first")
res_df$genename <- mapIds(org.Hs.eg.db,
                           keys      = res_df$ensembl,
                           column    = "GENENAME",
                           keytype   = "ENSEMBL",
                           multiVals = "first")

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# Export
dir.create("work_dir", showWarnings = FALSE)
write.csv(res_df, "work_dir/airway_deseq2_results.csv", row.names = FALSE)

# Preview
print(head(res_df[, c("symbol", "log2FoldChange", "padj", "genename")], 10))
```

</details>

---

## Task 8 (Bonus): Check Known GC-Responsive Genes

The Himes et al. paper validated several known glucocorticoid-responsive genes by qRT-PCR: **DUSP1**, **FKBP5**, **KLF15**, **TSC22D3**, and **PER1**. Check whether these genes are significantly differentially expressed in your results.

- Subset your annotated results to rows where `symbol` is one of the five genes above.
- Print their log2 fold change, p-value, and adjusted p-value.
- Are all five significant at FDR < 0.05? Are they all upregulated by dexamethasone?

<details>
<summary>💡 Hint</summary>

```r
known_genes <- c("DUSP1", "FKBP5", "KLF15", "TSC22D3", "PER1")
gc_results  <- res_df[res_df$symbol %in% known_genes, ]
print(gc_results[, c("symbol", "log2FoldChange", "pvalue", "padj")])
```

Cross-reference these results with Table 1 in the Himes et al. paper.

</details>
