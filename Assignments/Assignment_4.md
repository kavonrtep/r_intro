## Assignment: Hands-On Data Cleaning & Manipulation with R Tidyverse 


### Overview 

In this assignment, you will work with a gene expression dataset from Golub et al.'s Leukemia Expression Data. The data is stored in the `data/aml` folder and includes:
 
- **data_set_ALL_AML_train.csv:**  A CSV file with gene expression values and additional metadata.
- **annotation.csv:**  A CSV file mapping patients (samples) to cancer types.


Throughout this assignment, you will practice:
 
- Loading data using the tidyverse (`readr` and `dplyr`)
- Cleaning data by removing extraneous columns
- Filtering rows based on keyword detection with `stringr`
- Creating new columns with `mutate`
- Calculating row-wise summaries (average expression)
- Reshaping data between wide and long formats using `tidyr`
- Merging datasets using `left_join`
- Visualizing data with `ggplot2`
- Generating heatmaps using the `pheatmap` package

---
### Task 1: Load and Inspect the Training Data 

 
**Load the Dataset:** 
 
  - **Step:**  Load `data_set_ALL_AML_train.csv` into R using `read_csv()`.
  - **Instructions:** 
     - Read the data from the file path `"data/aml/data_set_ALL_AML_train.csv"`.
     - Inspect the dimensions with `dim()`, list the column names using `colnames()`, and preview the first few rows with `head()`.

<details>
<summary>💡 Hint</summary>


```R
# Load the training dataset
train_file <- "data/aml/data_set_ALL_AML_train.csv"
expr <- read_csv(train_file)

# Inspect the data structure
print(dim(expr))
print(colnames(expr)[1:10])
print(head(expr[, 1:5]))
```
</details>

---

### Task 2: Remove "call" Columns 

 
**Clean the Dataset:** 
  - **Step:**  Remove any columns that contain the string `"call"` because these represent detection status (not expression values).
  - **Instructions:** 
     - Use the `select()` function with `-contains("call")` from `dplyr` to remove these columns.
     - Verify the remaining columns and check the new dimensions of the dataset.
<details>
<summary>💡 Hint</summary>

```R
expr_clean <- expr %>% select(-contains("call"))
print(colnames(expr_clean)[1:10])
print(dim(expr_clean))
```
</details>

---

### Task 3: Filter Genes by Keyword 
 
**Extract Genes Based on Description:** 
  - **Step:**  Identify genes whose descriptions contain the keyword `"kinase"`.
  - **Instructions:** 
    - Use `filter()` along with `str_detect()` (from `stringr`) on the column `Gene Description`.
    - Print the number of kinase-related genes and view some of their descriptions.

<details>
<summary>💡 Hint</summary>

```R
kinase_genes <- expr_clean %>% filter(str_detect(`Gene Description`, "kinase"))
print(paste("Number of kinase-related genes:", nrow(kinase_genes)))
print(head(kinase_genes$`Gene Description`, 10))
```
</details>
 
**Optional – Filter Ribosomal Genes:** 
 
  - **Instructions:** 
    - Repeat the process for genes with `"ribosomal"` in their description and print the count.

<details>
<summary>💡 Hint</summary>

```R
ribosomal_genes <- expr_clean %>% filter(str_detect(`Gene Description`, "ribosomal"))
print(paste("Number of ribosomal-related genes:", nrow(ribosomal_genes)))
```
</details>

---
### Task 4: Flag Control Probes 
 
**Identify Control Probes:** 
 
  - **Step:**  Create a new Boolean column called `IsControl` that flags control probes. Control probes are identified by a `Gene Accession Number` starting with `"AFFX"`.
  - **Instructions:** 
     - Use `mutate()` along with `str_detect()` to add the `IsControl` column.
     - Display the first few rows of columns: `Gene Description`, `Gene Accession Number`, and `IsControl`, and summarize the counts.

<details>
<summary>💡 Hint</summary>

```R
expr_clean <- expr_clean %>%
  mutate(IsControl = str_detect(`Gene Accession Number`, "^AFFX"))

# Inspect the new column and summary
print(head(expr_clean %>% select(`Gene Description`, `Gene Accession Number`, IsControl)))
print(table(expr_clean$IsControl))
```
</details>

---

### Task 5: Compute Average Expression per Gene 
**Calculate Mean Expression:** 
 
  - **Step:**  Compute the average expression for each gene across all sample columns.
  - **Instructions:** 
 
    - Use `mutate()` together with `rowMeans()`. Exclude non-numeric columns such as `Gene Description`, `Gene Accession Number`, and `IsControl`.
     - Print the first few rows to display `Gene Description` and the newly computed average expression (`AvgExpr`).

<details>
<summary>💡 Hint</summary>

```R
expr_clean <- expr_clean %>%
  mutate(AvgExpr = rowMeans(select(., -c(`Gene Description`, `Gene Accession Number`, IsControl)), na.rm = TRUE))

# View average expression for each gene
print(head(expr_clean %>% select(`Gene Description`, AvgExpr)))
```

</details>
 
**Summarize by Control Status:** 
 
  - **Instructions:** 
 
    - Group the data by `IsControl` and calculate the mean and standard deviation of `AvgExpr` using `group_by()` and `summarize()`.

<details>
<summary>💡 Hint</summary>


```R
control_summary <- expr_clean %>%
  group_by(IsControl) %>%
  summarize(mean_expr = mean(AvgExpr), sd_expr = sd(AvgExpr))
print(control_summary)
```

</details>

---

### Task 6: Reshape Data from Wide to Long Format 

**Transform the Dataset:** 
 
  - **Step:**  Convert the dataset from wide format (one row per gene) to long format (one row per gene-sample measurement).
   - **Instructions:** 
     - Use `pivot_longer()` from `tidyr` to collapse all sample columns into two columns: `Sample` and `Expression`.
     - Exclude the metadata columns (`Gene Description`, `Gene Accession Number`, `IsControl`, and `AvgExpr`).
     - Check the first few rows and dimensions of the reshaped data.

<details>
<summary>💡 Hint</summary>


```R
expr_long <- expr_clean %>%
  pivot_longer(
    cols = -c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr),
    names_to = "Sample",
    values_to = "Expression"
  )
print(head(expr_long))
print(dim(expr_long))
```

</details>

---

### Task 7: Join "Call" Information Back 

 
2. **Reintegrate Detection Calls:** 
 
  - **Step:**  The original dataset includes detection calls in columns whose names contain `"call"`. Join these back into the long-format data.
  - **Instructions:** 
 
    - Select the call-related columns from the original `expr` data.
    - Use `pivot_longer()` to reshape the call data into a long format with columns `SampleCall` and `Call`.
    - Clean the sample names (remove the `"call"` prefix and adjust if needed).
    - Use `left_join()` to merge the call data with `expr_long` based on `Gene Accession Number` and `Sample`.
     - Inspect the resulting dataset and summarize the `Call` frequencies.

<details>
<summary>💡 Hint</summary>


```R
# Reshape call columns from original dataset
calls_long <- expr %>%
  select(`Gene Accession Number`, contains("call")) %>%
  pivot_longer(
    cols = -`Gene Accession Number`,
    names_to = "SampleCall",
    values_to = "Call"
  )

# Simplify the sample names by removing the prefix
calls_long <- calls_long %>%
  mutate(Sample = str_remove(SampleCall, "^call[.]+")) %>%
  select(-SampleCall)

# Adjust sample numbering if needed
f <- function(x) { (x - 2) / 2 }
calls_long$Sample <- calls_long$Sample %>% as.numeric() %>% f() %>% as.character()

# Join with long-format expression data
expr_long <- expr_long %>% left_join(calls_long, by = c("Gene Accession Number", "Sample"))
print(head(expr_long))
print(table(expr_long$Call, useNA = "ifany"))
```

</details>

---



### Task 8: Visualize the Data with ggplot2 

 
2. **Boxplot of Expression by Call Category:** 
 
  - **Step:**  Create a boxplot to visualize how expression values vary by detection call category.
  - **Instructions:** 
    - Use `ggplot2` to plot a boxplot of `Expression` versus `Call`.
    - Apply a log10 transformation to the y-axis for better visualization.

<details>
<summary>💡 Hint</summary>


```R
p1 <- ggplot(expr_long, aes(x = Call, y = Expression)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(trans = 'log10') +
  labs(title = "Expression Distribution by Call Category",
       x = "Call (Detection Flag)", y = "Expression (log10)") +
  theme_minimal()
print(p1)
```

</details>
 
**Bar Plot of Top 10 Kinase Genes:** 
 
  - **Step:**  Create a horizontal bar plot for the top 10 kinase genes based on average expression.
  - **Instructions:** 
 
    - From the cleaned data, filter for kinase genes and sort them by `AvgExpr` in descending order.
    - Use `geom_col()` and `coord_flip()` to generate a horizontal bar plot.

<details>
<summary>💡 Hint</summary>


```R
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
```

</details>



---



### Task 9: Incorporate Annotation Data for Stratified Visualization 

 
**Merge Patient Annotation:** 
 
  - **Step:**  Load the annotation file and merge it with your long-format expression data.
  - **Instructions:** 
 
    - Read `annotation.csv` using `read_csv()`.
    - Rename the `patient` column to `Sample` and convert it to character.
    - Join the annotation data with `expr_long` using `left_join()`.
    - Verify the merged data by previewing the first few rows.

<details>
<summary>💡 Hint</summary>


```R
# Load annotation data
annotation <- read_csv("data/aml/annotation.csv")
annotation <- annotation %>% rename(Sample = patient)
annotation$Sample <- as.character(annotation$Sample)

# Merge annotation data with expression data
expr_long <- expr_long %>% left_join(annotation, by = "Sample")
print(head(expr_long))
```

</details>
 
**Faceted Plot for Kinase Genes:** 
 
  - **Step:**  Focus on kinase genes and visualize expression by cancer type.
 
  - **Instructions:** 
 
    - Filter `expr_long` to include only rows for kinase genes.
    - From `expr_clean`, identify the top 5 kinase genes by average expression.
    - Further filter the long-format data for these top 5 genes.
    - Create a faceted boxplot (using `facet_wrap()`) with jittered points to show the expression distribution by cancer type.

<details>
<summary>💡 Hint</summary>


```R
kinase_expr_long <- expr_long %>% filter(str_detect(`Gene Description`, "kinase"))

# Identify top 5 kinase genes
top5_kinases <- expr_clean %>%
  filter(str_detect(`Gene Description`, "kinase")) %>%
  arrange(desc(AvgExpr)) %>%
  slice_head(n = 5) %>%
  pull(`Gene Description`)

# Filter long data for these top 5 genes
kinase_expr_long_subset <- kinase_expr_long %>%
  filter(`Gene Description` %in% top5_kinases)

# Create faceted boxplot
p4 <- ggplot(kinase_expr_long_subset, aes(x = cancer, y = Expression, fill = cancer)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~`Gene Description`) +
  labs(title = "Expression of Top 5 Kinase Genes by Cancer Type",
       x = "Cancer Type", y = "Expression (log10)") +
  theme_minimal() +
  theme(legend.position = "none")
print(p4)
```

</details>



---



### Task 10: Generate Heatmaps Using pheatmap 

 
**Compute Variability Metrics:** 
 
  - **Step:**  Calculate variance, standard deviation, and coefficient of variation (CoeffVar) for each gene.
  - **Instructions:** 
    - Exclude non-numeric columns and compute the metrics using `apply()` and `mutate()`.

<details>
<summary>💡 Hint</summary>


```R
expr_clean <- expr_clean %>%
  mutate(Variance = apply(select(., -c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr)), 1, var),
         SD = apply(select(., -c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr)), 1, sd),
         CoeffVar = SD / AvgExpr)
```

</details>
 
4. **Heatmap of Top 100 Variable Genes:** 
 
  - **Step:**  Select the top 100 genes with the highest coefficient of variation and generate a heatmap.
  - **Instructions:** 
 
    - Filter the dataset based on `CoeffVar` (highest first) and slice the top 100 rows.
    - Extract only the numeric expression data.
    - Use `pheatmap()` to create a heatmap with row scaling.

<details>
<summary>💡 Hint</summary>


```R
top100_var_genes <- expr_clean %>%
  arrange(desc(CoeffVar)) %>%
  slice_head(n = 100)

top100_var_expr <- top100_var_genes %>%
  select(-c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr, Variance, SD, CoeffVar))

pheatmap(top100_var_expr, scale = "row", show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 100 Genes with Highest Variance")
```

</details>
 
6. **Additional Heatmap Variants:** 
 
  - **Heatmap with Annotation:** 
    - Add cancer type annotations to the heatmap.

<details>
<summary>💡 Hint</summary>


```R
pheatmap(top100_var_expr, scale = "row", show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 100 Genes with Highest Variance by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer))
```

</details>
 
  - **Heatmap for Top 100 Genes by Average Expression:** 
 
    - Select the top 100 genes based on `AvgExpr` and generate a heatmap with correlation-based clustering.

<details>
<summary>💡 Hint</summary>


```R
top100AvgExpr_annotated <- expr_clean %>%
  arrange(desc(AvgExpr)) %>%
  slice_head(n = 100)

top100AvgExpr <- top100AvgExpr_annotated %>%
  select(-c(`Gene Description`, `Gene Accession Number`, IsControl, AvgExpr, Variance, SD, CoeffVar))

pheatmap(top100AvgExpr, scale = "row", show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D", clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Top 100 Genes with Highest Average Expression by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer))
```

</details>
 
  - **Heatmap with Gene Names and Custom Colors:** 
 
    - Display gene names as row labels and use a blue-white-red palette.

<details>
<summary>💡 Hint</summary>


```R
gene_description <- top100AvgExpr_annotated$`Gene Description`
pheatmap(top100AvgExpr, scale = "row", labels_row = gene_description, show_rownames = TRUE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D", clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Top 100 Genes with Highest Average Expression by Cancer Type",
         annotation_col = data.frame(Cancer = expr_long$cancer),
         color = colorRampPalette(c("blue", "white", "red"))(100))
```

</details>
 
  - **Heatmap Using All Genes (Excluding Control Probes):** 
 
    - Filter out control probes and visualize all remaining genes with cancer type annotations.

<details>
<summary>💡 Hint</summary>


```R
expr_clean_no_control <- expr_clean %>%
  filter(IsControl == FALSE) %>%  # Exclude control probes
  select(-IsControl) %>%
  select(-c(`Gene Description`, `Gene Accession Number`, AvgExpr, Variance, SD, CoeffVar))

p <- pheatmap(expr_clean_no_control, scale = "row", show_rownames = FALSE,
              cluster_rows = TRUE, cluster_cols = TRUE,
              clustering_method = "ward.D", clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              main = "All Genes by Cancer Type",
              annotation_col = data.frame(Cancer = expr_long$cancer),
              color = colorRampPalette(c("blue", "white", "red"))(100))
# Save the heatmap to a file
ggsave("work_dir/heatmap_all_genes.png", p, width = 10, height = 10, units = "in", dpi = 300)
```

</details>

---





