## Follow-Up Assignment: Gene Expression Filtering and Normalization 

### Overview 

In this assignment, you will work with gene expression data provided in two files:

- **feature_counts.csv:**  A tab-delimited file containing gene expression counts with gene names as row names and sample names as column names.
- **feature_annotation.csv:**  A tab-delimited file containing feature annotations, including a `gender` column.

You will practice:

- Loading and inspecting data
- Summarizing annotation data
- Computing basic statistics (such as counting genes and identifying those with zero expression)
- Filtering genes based on expression criteria
- Using the apply family to calculate per-gene statistics
- Running a differential expression test using a provided function
- Creating a simple base R plot to visualize your results
---


### Task 1: Data Loading and Inspection 
 
1. **Load the Datasets:**  
  - **Step:**  Load the expression and annotation data into R using `read.delim()`.
  - **Instructions:**  
    - Load `feature_counts.csv` and `feature_annotation.csv` from the `data/` directory. 
    - Use `read.delim() function.
    - Inspect the first few rows of both datasets with `head()`.
<details>
<summary>💡 Hint</summary>

```R
feature_counts <- read.delim("data/feature_counts.csv", header = TRUE, sep = "\t")
feature_annotation <- read.delim("data/feature_annotation.csv", header = TRUE, sep = "\t")

# Inspect the data
head(feature_counts)
head(feature_annotation)
```
</details>


---
### Task 2: Summarizing Feature Annotation 
 
1. **Summarize the Gender Information:**  
  - **Step:**  Use the `table()` function to summarize the `gender` column in `feature_annotation` (`table` function is used to create a frequency table of the variable).
  - **Instructions:** 
    - Create a summary table of genders and print the result.
<details>
<summary>💡 Hint</summary>

```R
gender_summary <- table(feature_annotation$gender)
print(gender_summary)
```
</details>
---


### Task 3: Counting Total Genes 
 
1. **Determine the Number of Genes:**  
  - **Step:**  Use the `nrow()` function on `feature_counts` to find out how many genes are in the dataset.
  - **Instructions:** 
    - Print the total number of genes.
<details>
<summary>💡 Hint</summary>

```R
total_genes <- nrow(feature_counts)
cat("Total number of genes:", total_genes, "\n")
```
</details>


---
### Task 4: Identifying Genes with Zero Expression 
Feature counts are often sparse, with many genes having zero expression in certain samples. In this task, you will identify genes with zero expression across all samples. Such genes may be filtered out in downstream analyses.
1. **Count Genes with Zero Expression Across All Samples:**  
  - **Step:**  Identify genes that have zero counts in every sample.
  - **Instructions:**  
    - Use a vectorized approach (e.g., `rowSums()`) to determine which genes have all zero values.
    - Use comparison operator "==" and `sum` function to count the number of genes with zero expression in all samples.
    - Print the number of such genes.
<details>
<summary>💡 Hint</summary>

```R
# Calculate how many genes have zero expression in all samples
zero_expression_genes <- sum(rowSums(feature_counts == 0) == ncol(feature_counts))
cat("Number of genes with zero expression in all samples:", zero_expression_genes, "\n")
```
</details>

---

### Task 5: Filtering Genes Based on Expression 
 
1. **Filter Genes with Positive Counts in at Least 10 Samples:**  
  - **Step:**  Create a filtered version of `feature_counts` that retains only the genes which have a positive count in at least 10 samples.
 
  - **Instructions:**  
    - Use a logical condition with `rowSums()` to filter the dataset.
    - Use either `subset()` or logical indexing to create the filtered dataset.
    - Determine and print the number of genes retained after filtering.
<details>
<summary>💡 Hint</summary>

```R
# Filter to keep genes with positive expression in at least 10 samples
filtered_counts <- feature_counts[rowSums(feature_counts > 0) >= 10, ]
kept_genes <- nrow(filtered_counts)
cat("Number of genes kept after filtration:", kept_genes, "\n")
```
</details>


---


### Task 6: Calculating Gene Expression Statistics 
 
1. **Compute Standard Deviation and Mean Expression:**  
  - **Step:**  Using the `apply()` function, calculate the standard deviation and mean for each gene (i.e., each row) across all samples in the filtered dataset. 
 
  - **Instructions:** 
    - Compute the per-gene mean and standard deviation. Stanard deviation is calculated using `sd()` function.
    - Display the first few results for each statistic.

<details>
<summary>💡 Hint</summary>

```R
# Calculate per-gene mean expression
gene_means <- apply(filtered_counts, 1, mean)
# Calculate per-gene standard deviation
gene_sds <- apply(filtered_counts, 1, sd)

# Display the first few statistics
cat("First few gene means:\n")
print(head(gene_means))
cat("First few gene standard deviations:\n")
print(head(gene_sds))
```
</details>

2. **Calculate standard deviation using `for` loop**  
  - **Step:**  Calculate the standard deviation for each gene using a `for` loop.
 
  - **Instructions:** 
    - Create an empty vector to store the standard deviations.
    - Use a `for` loop to iterate over each row of the filtered dataset.
    - Calculate the standard deviation for each gene and store the result in the vector.
    - Display the first few standard deviations.
<details>
<summary>💡 Hint</summary>

```R
# Calculate standard deviation for each gene using a for loop
gene_sds_loop <- numeric(nrow(filtered_counts))
for (i in 1:nrow(filtered_counts)) {
  gene_sds_loop[i] <- sd(filtered_counts[i, ])
}
head(gene_sds_loop)
```
</details>
Which approach do you find more efficient for this task, `apply` or `for` loop? Why?

---


### Task 7: Differential Expression Analysis 
 
1. **Run the Differential Expression Test:**  
  - **Step:**  Use the provided `test_differential_expression` function to perform a t-test comparing expression between genders.
 
  - **Instructions:**  
    - Run the function on `feature_counts` (or on the filtered dataset if preferred) along with the `gender` vector from `feature_annotation`.

    - Inspect the first few rows of the resulting data frame.
<details>
<summary>💡 Hint</summary>

```R
test_differential_expression <- function(expr_data, groups) {
  expr_data <- as.matrix(expr_data)
  group_means <- apply(expr_data, 1, function(x) {
    tapply(x, groups, mean)
  })
  # t test for each gene
  p_value <- apply(expr_data, 1, function(x) {
    res <- t.test(x ~ groups)
    return(res$p.value)
  })
  results <- data.frame(t(group_means), p_value)
  return(results)
}

# Run the differential expression test
DE_results <- test_differential_expression(feature_counts, feature_annotation$gender)
head(DE_results)
```
</details>
 
2. **Identify Differentially Expressed Genes:**  
  - **Step:**  Determine how many genes are differentially expressed with at least a two-fold difference between genders.
 
  - **Instructions:**  
    - For each gene in `DE_results`, compute the fold change as the ratio between the higher and lower group means.
    - Count and print the number of genes with a fold change of at least 2.
<details>
<summary>💡 Hint</summary>

```R
# Calculate fold change for each gene (using the maximum and minimum between the two groups)
fold_changes <- apply(DE_results[, c("female", "male")], 1, function(x) {
  max(x) / min(x)
})

# Count genes with a fold change of at least 2
diff_expr_genes <- sum(fold_changes >= 2, na.rm = TRUE)
cat("Number of differentially expressed genes (>=2-fold difference):", diff_expr_genes, "\n")
```
</details>


---

### Task 8: Basic Plotting 
 
1. **Visualize Expression Distribution:**  
  - **Step:**  Create a histogram of the mean expression values (computed in Task 6) for the filtered genes using base R graphics.
 
  - **Instructions:** 
    - Plot a histogram with appropriate axis labels and a title.

<details>

<summary>💡 Hint</summary>

```R
# Plot a histogram of gene means
hist(gene_means,
     main = "Distribution of Mean Gene Expression",
     xlab = "Mean Expression",
     ylab = "Frequency",
     col = "lightblue",
     border = "black")
```

</details>

---


## Submission 
 
- **Script File:**  Save your completed assignment in an R script file (e.g., `gene_expression_assignment.R`).
- **Error-Free:**  Ensure your script runs without errors and that each task produces the expected output.
- **Submission:**  Send your R script file as an attachment via email.


