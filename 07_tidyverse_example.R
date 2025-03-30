# ---------------------------------------------------------------------
# Hands-On Data Cleaning & Manipulation with R Tidyverse
# Molecular Biology Example using Golub et al. Leukemia Expression Data
#
# This script demonstrates:
#   - Reading in a real-world gene expression dataset.
#   - Cleaning data by removing extraneous columns.
#   - Filtering rows using string detection.
#   - Creating new columns with mutate (e.g., flagging control probes).
#   - Computing row-wise summaries (average expression).
#   - Reshaping data between wide and long formats.
#   - Simple string manipulations.
#   - Creating basic ggplot2 visualizations.
#   - Generating heatmaps with the pheatmap package.
#
# The dataset is assumed to be in the "data/aml" folder with the following files:
#   - data_set_ALL_AML_train.csv        (Training dataset)
#   - annotation.csv                    (Sample annotation mapping patients to cancer
#   types)
#
# Students are expected to use, test, and modify this code as a starting point.
# ---------------------------------------------------------------------

# Load the core tidyverse packages.
# The tidyverse provides a suite of tools for data manipulation (dplyr, tidyr), string
# operations (stringr),
# and visualization (ggplot2), which are essential for reproducible data analysis.
library(tidyverse)

# ===============================
# 1. Load and Inspect Training Data
# ===============================
# Read the training dataset from a CSV file.
# Ensure the file path is correct relative to your working directory.
train_file <- "data/aml/data_set_ALL_AML_train.csv"
expr <- read_csv(train_file)
# Tip: You can also open this CSV file in Excel or LibreOffice Calc to inspect its
# structure.

# Inspect the dataset:
# - dim(): shows the number of rows (genes) and columns (samples plus metadata).
# - colnames(): lists the names of the columns (only the first 10 are printed for
# brevity).
# - head(): prints the first few rows of selected columns to give an overview of the data.
print(dim(expr))
print(colnames(expr)[1:10])
print(head(expr[, 1:5]))
# Note: In this dataset, numeric columns represent expression values for samples.
#       Additional columns with "call" indicate detection status ("P", "M", "A") per
#       microarray software.

# ===============================
# 2. Remove "call" Columns
# ===============================
# The dataset includes paired columns: one for expression values and one for detection
# calls.
# For this exercise, we want to focus on the expression values.
# We remove any column whose name contains the string "call" using select() with
# -contains().
expr_clean <- expr %>% select(-contains("call"))
print(colnames(expr_clean)[1:10])
print(dim(expr_clean))
# Reasoning: Removing the call columns simplifies downstream analysis, focusing only on
# quantitative expression.

# ===============================
# 3. Filter Genes by Keyword
# ===============================
# We often want to study specific sets of genes. Here, we extract genes that mention
# "kinase" in their description.
# Backticks (`) are used around column names with spaces.
kinase_genes <- expr_clean %>%
  filter(str_detect(`Gene Description`, "kinase"))
print(paste("Number of kinase-related genes:", nrow(kinase_genes)))
print(head(kinase_genes$`Gene Description`, 10))
# Optional: You can try filtering for "ribosomal" genes similarly to explore different
# gene families.
ribosomal_genes <- expr_clean %>%
  filter(str_detect(`Gene Description`, "ribosomal"))
print(paste("Number of ribosomal-related genes:", nrow(ribosomal_genes)))
# This step reinforces the use of string detection (str_detect) for subsetting data
# based on textual patterns.

# ===============================
# 4. Flag Control Probes
# ===============================
# In microarray data, control probes (often used for quality control) have IDs starting
# with "AFFX".
# We add a new Boolean column "IsControl" that flags these rows.
expr_clean <- expr_clean %>%
  mutate(IsControl = str_detect(`Gene Accession Number`, "^AFFX"))
# Inspect the new column along with gene names and accession numbers.
print(head(expr_clean %>% select(`Gene Description`, `Gene Accession Number`, IsControl)))
print(table(expr_clean$IsControl))
# Reasoning: Identifying control probes helps in filtering or analyzing them separately
# from biological genes.

# ===============================
# 5. Compute Average Expression per Gene
# ===============================
# Here we calculate the average expression value across all sample columns for each gene.
# We exclude non-numeric columns ("Gene Description", "Gene Accession Number", and
# "IsControl") from the calculation.
expr_clean <- expr_clean %>%
  mutate(AvgExpr = rowMeans(select(., -c(`Gene Description`, `Gene Accession Number`,
                                         IsControl)), na.rm = TRUE))
print(head(expr_clean %>% select(`Gene Description`, AvgExpr)))
# We then summarize average expression by control status (IsControl) using group_by()
# and summarize().
control_summary <- expr_clean %>%
  group_by(IsControl) %>%
  summarize(mean_expr = mean(AvgExpr), sd_expr = sd(AvgExpr))
print(control_summary)
# Reasoning: Computing row-wise summaries helps in understanding the overall expression
# distribution and comparing groups.

# ===============================
# 6. Reshape Data: Wide to Long Format
# ===============================
# Many analyses (especially plotting) benefit from a long-format data structure,
# where each row corresponds to a single measurement.
# pivot_longer() is used to collapse all sample columns into two: "Sample" and
# "Expression".
expr_long <- expr_clean %>%
  pivot_longer(
    cols = -c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr),
    names_to = "Sample",
    values_to = "Expression"
  )
print(head(expr_long))
print(dim(expr_long))  # Expected: (number of genes) * (number of samples)
# A long format is ideal for ggplot2, which expects one observation per row.

# ===============================
# 7. Join "Call" Information Back
# ===============================
# If you wish to include the detection calls in your analysis, you can reshape and join
# them back.
# First, we select the call-related columns from the original data.
calls_long <- expr %>%
  select(`Gene Accession Number`, contains("call")) %>%
  pivot_longer(
    cols = -`Gene Accession Number`,
    names_to = "SampleCall",
    values_to = "Call"
  )
print(head(calls_long))
# Next, we simplify the sample names by removing the "call" prefix.
calls_long <- calls_long %>%
  mutate(Sample = str_remove(SampleCall, "^call[.]+")) %>%
  select(-SampleCall)
print(calls_long)
# The sample numbers might be encoded as "4", "6", "8", etc. We renumber them
# sequentially.
f <- function(x) { (x - 2) / 2 }
calls_long$Sample <- calls_long$Sample %>%
  as.numeric() %>%
  f() %>%
  as.character()
print(calls_long)
# Finally, we join the call information with the long-format expression data using
# left_join().
expr_long <- expr_long %>% left_join(calls_long, by = c("Gene Accession Number",
                                                        "Sample"))
print(head(expr_long))
print(table(expr_long$Call, useNA = "ifany"))
# Reasoning: This optional step allows you to later assess how detection calls relate
# to expression values.

# ===============================
# 8. Visualize the Data with ggplot2
# ===============================
# Plot 1: Boxplot of Expression Distribution by Call Category.
# This plot shows how expression values vary based on whether a call was "P", "M", or "A".
# A log10 transformation is applied to the y-axis to better visualize a wide range of
# expression values.
p1 <- ggplot(expr_long, aes(x = Call, y = Expression)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(trans = 'log10') +
  labs(title = "Expression Distribution by Call Category",
       x = "Call (Detection Flag)", y = "Expression (log10)") +
  theme_minimal()
print(p1)


# Plot 2: Bar Plot of Top 10 Kinase Genes by Average Expression.
# We extract the top 10 kinase genes (based on the previously computed AvgExpr),
# reorder them for clarity, and then plot a horizontal bar chart.
top10_kinases <- expr_clean %>%
  filter(str_detect(`Gene Description`, "kinase")) %>%
  arrange(desc(AvgExpr)) %>%
  slice_head(n = 10)
p2 <- ggplot(top10_kinases, aes(x = reorder(`Gene Description`, AvgExpr), y = AvgExpr)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Kinase Genes by Average Expression",
       x = "Gene (Description)", y = "Average Expression") +
  theme_minimal()
print(p2)
# Reasoning: This plot highlights the most highly expressed kinase genes, which might
# be of biological interest.

# ===============================
# 9. Incorporate Annotation Data for Additional Visualization
# ===============================
# Load the annotation file that maps each patient (sample) to a cancer type (e.g., ALL
# or AML).
annotation <- read_csv("data/aml/annotation.csv")
print(head(annotation))
# Rename the 'patient' column to 'Sample' to match the sample identifiers in expr_long.
annotation <- annotation %>% rename(Sample = patient)
# Convert Sample to character to ensure consistency for joining.
annotation$Sample <- annotation$Sample %>% as.character()

# Join the annotation data with the long-format expression data.
expr_long <- expr_long %>% left_join(annotation, by = "Sample")
print(head(expr_long))
# Reasoning: Adding cancer type information enables stratified analyses and more
# informative visualizations.

# ------------------------------
# Subset and Facet Plot for Kinase Genes
# ------------------------------
# We now focus on a subset of the data: kinase genes.
# First, filter expr_long to include only rows with "kinase" in the gene description.
kinase_expr_long <- expr_long %>%
  filter(str_detect(`Gene Description`, "kinase"))

# Next, select the top 5 kinase genes by average expression from expr_clean.
# This selection focuses our visualization on the most prominent kinases.
top5_kinases <- expr_clean %>%
  filter(str_detect(`Gene Description`, "kinase")) %>%
  arrange(desc(AvgExpr)) %>%
  slice_head(n = 5) %>%
  pull(`Gene Description`)

# Further filter the long-format data to include only these top 5 genes.
kinase_expr_long_subset <- kinase_expr_long %>%
  filter(`Gene Description` %in% top5_kinases)

# Create a faceted boxplot with jittered points:
# - x-axis: cancer type (from the annotation)
# - y-axis: Expression (log10-transformed)
# - Facets: separate panel for each gene (top 5 kinases)
# - Color fill differentiates cancer types
p4 <- ggplot(kinase_expr_long_subset, aes(x = cancer, y = Expression, fill = cancer)) +
  geom_boxplot(outlier.shape = NA) +        # Boxplot summarizes the distribution per
  # cancer type
  geom_jitter(width = 0.2, alpha = 0.6) +     # Jitter adds individual data points for
  # detailed insight
  facet_wrap(~`Gene Description`) +          # Facet by gene to compare distributions
  # side-by-side
  labs(title = "Expression of Top 5 Kinase Genes by Cancer Type",
       x = "Cancer Type", y = "Expression (log10)") +
  theme_minimal() +
  theme(legend.position = "none")             # Legend is redundant here because x-axis
# labels already denote cancer types
print(p4)
# Reasoning: Faceting and color-coding enable visual comparisons across multiple gene
# subsets and conditions.

# ===============================
# 10. Generate Heatmaps using pheatmap
# ===============================
# Load the pheatmap package (install if not available: install.packages("pheatmap"))
library(pheatmap)

# ---- Heatmap Example: Top 100 Genes by Variability ----
# Here, we calculate variance and standard deviation for each gene to identify genes
# with high variability.
# Variability can reveal genes that are differentially expressed across samples.
expr_clean <- expr_clean %>%
  mutate(Variance = apply(select(., -c(`Gene Description`, `Gene Accession Number`,
                                       IsControl, AvgExpr)), 1, var),
         SD = apply(select(., -c(`Gene Description`, `Gene Accession Number`,
                                 IsControl, AvgExpr)), 1, sd),
         CoeffVar = SD / AvgExpr)
# Select the top 100 genes with the highest coefficient of variation.
top100_var_genes <- expr_clean %>%
  arrange(desc(CoeffVar)) %>%
  slice_head(n = 100)
# Extract only the expression data (numeric columns) for these genes.
top100_var_expr <- top100_var_genes %>%
  select(-c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr, Variance,
            SD, CoeffVar))
# Generate a heatmap. We scale rows to standardize each gene's expression.
pheatmap(top100_var_expr, scale = "row", show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 100 Genes with Highest Variance")
# Reasoning: A heatmap of variable genes may reveal clusters of samples or genes with
# similar expression patterns.

# ---- Heatmap with Annotation: Top 100 Genes by Variability ----
# Now, we add cancer type annotations to the heatmap columns.
pheatmap(top100_var_expr, scale = "row", show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 100 Genes with Highest Variance by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer))


# ---- Heatmap Example: Top 100 Genes by Average Expression ----
# We select the top 100 genes based on average expression.
top100AvgExpr_annotated <- expr_clean %>%
  arrange(desc(AvgExpr)) %>%
  slice_head(n = 100)
# Extract the expression data.
top100AvgExpr <- top100AvgExpr_annotated %>%
  select(-c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr, Variance,
            SD, CoeffVar))
# Generate a heatmap with correlation-based clustering.
pheatmap(top100AvgExpr, scale = "row", show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D", clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Top 100 Genes with Highest Average Expression by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer))
# Clustering based on correlation can help identify groups of genes with similar
# expression profiles.

# ---- Heatmap with Gene Names ----
# To display gene names on the heatmap, first extract the gene descriptions.
gene_description <- top100AvgExpr_annotated$`Gene Description`
# Use make.unique() if necessary to ensure unique labels.
# Then, pass these labels via the 'labels_row' argument in pheatmap.
pheatmap(top100AvgExpr, scale = "row", labels_row = gene_description, show_rownames =
  TRUE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D", clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Top 100 Genes with Highest Average Expression by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer))


# ---- Heatmap with Custom Color Palette ----
# The following heatmap uses a blue-white-red palette to emphasize low-to-high
# expression changes.
pheatmap(top100AvgExpr, scale = "row", labels_row = gene_description, show_rownames =
  TRUE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D", clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Top 100 Genes with Highest Average Expression by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer),
         color = colorRampPalette(c("blue", "white", "red"))(100)
)


# ---- Heatmap Using All Genes (Excluding Control Probes) ----
# Sometimes you may wish to visualize the full expression dataset excluding control
# probes.
expr_clean_no_control <- expr_clean %>%
  filter(IsControl == FALSE) %>%  # Exclude control probes
  select(-IsControl)
# Remove metadata columns so that only expression data remains.
expr_clean_no_control <- expr_clean_no_control %>%
  select(-c(`Gene Description`, `Gene Accession Number`, AvgExpr, Variance, SD, CoeffVar))
# Generate a heatmap of all genes, adding cancer type annotations.
p <- pheatmap(expr_clean_no_control, scale = "row", show_rownames = FALSE,
              cluster_rows = TRUE, cluster_cols = TRUE,
              clustering_method = "ward.D", clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              main = "All Genes by Cancer Type",
              annotation_col = data.frame(Cancer = expr_long$cancer),
              color = colorRampPalette(c("blue", "white", "red"))(100)
)
# Save the heatmap to a file for future reference or publication.
ggsave("work_dir/heatmap_all_genes.png", p, width = 10, height = 10, units = "in", dpi
  = 300)
