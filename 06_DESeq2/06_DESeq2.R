############################################################
# Introductory RNA-Seq Analysis with DESeq2 on the pasilla Dataset
#
# This lesson demonstrates:
# - Data loading and quality filtering
# - Differential expression analysis with DESeq2
# - Exploratory data analysis using PCA, MDS, and pairs plots
# - Interactive visualizations with plotly and heatmaply
#
# STUDENT TASKS are provided as comments within the code.
#
# Refer to the DESeq2 vignette for additional details:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
############################################################

# -------------------------------
# 1. Setup and Load Libraries
# -------------------------------
# Load required libraries. Each package serves a specific purpose:
# - DESeq2: Main package for RNA-Seq differential expression analysis.
# - ggplot2: Used for creating elegant and customizable static plots.
# - plotly: Provides interactive graphics for a more engaging data exploration.
# - matrixStats: Efficient functions for computing row statistics.
library(DESeq2)       # Differential expression analysis
library(ggplot2)      # Static plotting
library(plotly)       # Interactive plotting
library(matrixStats)  # For computing row variances

# -------------------------------
# 2. Data Loading and Preparation
# -------------------------------
# Load the count matrix and sample annotation data.
# - countData: Contains raw gene expression counts (rows: genes, columns: samples).
# - colData: Contains sample metadata (e.g., condition, sequencing type).
# Update file paths as needed.
countData <- read.delim("data/pasilla_gene_counts.tsv", row.names = 1)
colData <- read.csv("data/pasilla_sample_annotation.csv", row.names = 1)

# Background on the pasilla dataset:
# This dataset comes from an RNA-Seq experiment (Brooks et al. 2010) that investigates
# the impact of siRNA knockdown of pasilla—a gene implicated in mRNA splicing regulation.
# It includes multiple biological replicates for both the knockdown and the untreated control.

# Ensure the 'condition' variable is a factor (categorical variable).
colData$condition <- as.factor(colData$condition)

# Quick inspection of the data:
head(countData)
head(colData)

# IMPORTANT: The order of columns in countData must match the rows in colData.
# DESeq2 does not automatically match these; if they are misaligned, the analysis will be incorrect.
all(rownames(colData) == colnames(countData))  # Should return TRUE
# If the order is not correct, reorder colData to match countData.
colData <- colData[colnames(countData), ]
all(rownames(colData) == colnames(countData))  # Confirm the alignment

# -------------------------------
# 3. Create DESeqDataSet and Filter Low Count Genes
# -------------------------------
# Create the DESeqDataSet object, which combines the count data and the experimental design.
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# 'dds' is now a specialized object for DESeq2 analyses.

# Add gene names as metadata. This extra information is useful for later annotation steps.
metadata <- data.frame(gene = rownames(countData))
mcols(dds) <- DataFrame(mcols(dds), metadata)
mcols(dds)  # View the metadata to ensure gene names are attached

# Filter out genes with low counts across all samples to reduce noise.
# Low count genes can lead to unreliable statistical estimates.
# Here, we remove genes with a total count of less than 3.
# TASK: Plot a histogram of the row sums to examine the count distribution.
keep <- rowSums(counts(dds)) >= 3

# TASK: Use table() to determine how many genes have total counts less than 3.
# table(rowSums(counts(dds)) < 3)

# Retain only the genes that meet the threshold.
dds <- dds[keep, ]
cat("Number of genes after filtering:", nrow(dds), "\n")

# Set the reference level for 'condition'. This is crucial for proper interpretation of contrasts.
print(dds$condition)
# If the default (alphabetical) order is not desired, set the levels manually.
dds$condition <- factor(dds$condition, levels = c("untreated", "treated"))
# Alternatively, relevel so that 'untreated' is the reference.
dds$condition <- relevel(dds$condition, ref = "untreated")

# NOTE: For technical replicates (multiple sequencing runs of the same library),
# DESeq2 offers the collapseReplicates() function. Do not use this function to merge biological replicates.

# -------------------------------
# 4. Run the DESeq2 Analysis
# -------------------------------
# Run the full DESeq2 pipeline:
# 1. Normalization (estimating size factors)
# 2. Dispersion estimation for each gene
# 3. Fitting a negative binomial model for each gene
# 4. Hypothesis testing for differential expression
dds <- DESeq(dds)

# Extract results for differential expression.
# By default, results() uses an FDR threshold of 0.1.
res <- results(dds)
summary(res)

# You can adjust the FDR threshold. Here, we set it to 0.05.
res05 <- results(dds, alpha = 0.05)
summary(res05)

# -------------------------------
# 5. Variance Stabilizing Transformation (VST)
# -------------------------------
# RNA-Seq count data can span several orders of magnitude with variance increasing with the mean.
# VST transforms the data such that the variance becomes roughly constant across different mean values.
# This makes the data more suitable for visualizations like PCA, MDS, and heatmaps.
vsd <- vst(dds, blind = FALSE)
# 'vsd' now holds the variance-stabilized data.

# -------------------------------
# 6. Exploratory Data Analysis (EDA)
# -------------------------------
# 6a. Principal Component Analysis (PCA)
# -------------------------------
# PCA reduces the dimensionality of the data while preserving variance.
# We extract the PCA data and include the 'type' (sequencing type) for additional context.
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pcaData$type <- colData$type
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create a static PCA plot using ggplot2.
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of pasilla Samples") +
  theme_minimal()
print(p_pca)

# Generate an interactive version of the PCA plot using plotly.
p_pca_interactive <- ggplotly(p_pca)
print(p_pca_interactive)

# 6b. Multidimensional Scaling (MDS)
# -------------------------------
# MDS is used to visualize the similarity or dissimilarity of data points.
# Here, we compute the distance matrix based on the transformed data.
dists <- dist(t(assay(vsd)))
print(dists)  # Inspect the computed distances
mds <- cmdscale(dists)
mdsData <- data.frame(MDS1 = mds[, 1],
                      MDS2 = mds[, 2],
                      condition = colData$condition,
                      type = colData$type)

# Create a static MDS plot with ggplot2.
p_mds <- ggplot(mdsData, aes(MDS1, MDS2, color = condition, shape = type)) +
  geom_point(size = 3) +
  ggtitle("MDS Plot of pasilla Samples") +
  theme_minimal()
print(p_mds)

# Generate an interactive MDS plot using plotly.
p_mds_interactive <- ggplotly(p_mds)
print(p_mds_interactive)

# 6c. Pairs Plot
# -------------------------------
# A pairs plot shows scatterplots for every pair of samples.
# Here, we use the full transformed data. You might choose to focus on the most variable genes.
vsd_df <- as.data.frame(assay(vsd))
head(vsd_df)

# Create a pairs plot with transparency to handle overplotting.
# The color "#00000010" sets a low opacity (10% opacity) for the points.
pairs(vsd_df, col = "#00000010")


# -------------------------------
# 7. Visualization of Differential Expression Results
# -------------------------------
# 7a. MA Plot (Base R Plot)
# -------------------------------
# The MA plot displays the relationship between the log2 fold changes (M) and the mean of normalized counts (A).
# Significant genes (adjusted p-value below threshold) are typically highlighted.
plotMA(res, main = "DESeq2 MA-Plot", ylim = c(-2, 2))
plotMA(res, main = "DESeq2 MA-Plot", ylim = c(-5, 5))
plotMA(res05, main = "DESeq2 MA-Plot (FDR < 0.05)", ylim = c(-5, 5))

# 7b. Volcano Plot using ggplot2
# -------------------------------
# A volcano plot combines fold-change and statistical significance to highlight interesting genes.
volcanoData <- as.data.frame(res)
volcanoData$gene <- rownames(volcanoData)

p_volcano <- ggplot(volcanoData, aes(x = log2FoldChange, y = -log10(padj), text = gene)) +
  geom_point(alpha = 0.4) +
  ggtitle("Volcano Plot") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-Value") +
  theme_minimal()
print(p_volcano)

# Make the volcano plot interactive with plotly, so students can hover over points to see gene IDs.
p_volcano_interactive <- ggplotly(p_volcano, tooltip = c("text", "x", "y"))
print(p_volcano_interactive)

# -------------------------------
# 8. Accounting for Technical Variation
# -------------------------------
# PCA and MDS indicated that sequencing type (single-read vs. paired-end) affects the data.
# To correct for this technical variation, include 'type' as a covariate in the design.

colData$type <- as.factor(gsub("-.+","", colData$type))
print(colData)

dds2 <- DESeqDataSetFromMatrix(countData = countData,
                               colData = colData,
                               design = ~ type + condition)
# Ensure the reference level for condition is correctly set.
dds2$condition <- relevel(dds2$condition, ref = "untreated")
# Run the DESeq2 pipeline on the new dataset.
dds2 <- DESeq(dds2)

# Extract results comparing untreated and treated, accounting for sequencing type.
res2 <- results(dds2, contrast = c("condition", "untreated", "treated"))
summary(res2)

# Additionally, extract results to see the effect of sequencing type.
res2_type <- results(dds2, contrast = c("type", "paired", "single"))
summary(res2_type)

# Extract results with a more stringent FDR threshold (alpha = 0.05).
res2_05 <- results(dds2, alpha = 0.05, contrast = c("condition", "untreated", "treated"))
summary(res2_05)

# STUDENT TASK: Compare the number of differentially expressed genes (FDR < 0.05)
# from this analysis (accounting for 'type') with the previous analysis.


# Export the results sorted by adjusted p-value.
res2_05_reordered <- res2_05[order(res2_05$padj), ]
write.csv(res2_05_reordered, "work_dir/deseq2_results.csv", row.names = TRUE)

# -------------------------------
# 9. Additional Visualizations: Heatmaps
# -------------------------------
# Generate normalized and transformed data for heatmap visualization.
vsd2 <- vst(dds2, blind = FALSE)
ntd2 <- normTransform(dds2)  # Simple normalization for heatmap visualization

# Create a heatmap of the top 20 expressed genes using pheatmap.
library(pheatmap)
select <- order(rowMeans(counts(dds2, normalized = TRUE)), decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds2)[, c("condition", "type")])
pheatmap(assay(ntd2)[select, ], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df)

# Create a heatmap using hierarchical clustering to visualize differentially expressed genes.
select <- res2_05$padj < 0.01
select[is.na(select)] <- FALSE  # Replace NA values in padj with FALSE
df <- as.data.frame(colData(dds2)[, c("condition", "type")])
pheatmap(assay(ntd2)[select, ], cluster_rows = TRUE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df, scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100))
# TASK: Export this heatmap as a PNG file (e.g., using png(), pheatmap(), and dev.off()).

# Create an interactive heatmap using heatmaply.

# NOTE: it may be necessary to install the 'heatmaply' and 'htmlwidgets' packages.
# use command:
# install.packages(c("heatmaply", "htmlwidgets"))




library(heatmaply)
library(htmlwidgets)
p <- heatmaply(assay(ntd2)[select, ], scale = "row", k_row = 3, k_col = 3)
# Save the interactive heatmap as an HTML file for sharing.
saveWidget(p, "work_dir/heatmaply.html")

# -------------------------------
# 10. Gene Annotation
# -------------------------------
# Retrieve gene annotations (e.g., gene symbols and descriptions) from FlyBase.
# This uses the org.Dm.eg.db package from Bioconductor.
BiocManager::install("org.Dm.eg.db")  # Install if needed; this may take a while.
library(org.Dm.eg.db)
library(AnnotationDbi)
fly_ids <- rownames(countData)
# Check available annotation columns using columns(org.Dm.eg.db)
fly_genes <- select(org.Dm.eg.db, keys = fly_ids, columns = c("SYMBOL", "GENENAME"),
                    keytype = "FLYBASE")

# Merge the annotation data with the differential expression results.
res2_05_annotated <- merge(as.data.frame(res2_05_reordered),
                           fly_genes, by.x = "row.names",
                           by.y = "FLYBASE", all.x = TRUE)
# Sort the annotated results by adjusted p-value.
res2_05_annotated <- res2_05_annotated[order(res2_05_annotated$padj), ]

# Export the annotated results to a CSV file.
write.csv(res2_05_annotated, "work_dir/deseq
2_results_annotated.csv", row.names = FALSE)

# TASK: create new data frame which will include res2_05_annotated and  also expression
# values for each gene in each sample
# use assay(vsd2) to get normalized and transformed data, then merge it with res2_05_annotated