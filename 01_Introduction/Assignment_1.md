## Follow-Up Assignment: Bacterial Growth Curve Analysis

### Background

Your lab is studying the role of the **dnaA** gene in *Escherichia coli* replication initiation. The DnaA protein binds to the chromosomal origin of replication (*oriC*) and is essential for starting a new round of DNA replication. A colleague has constructed a *dnaA* temperature-sensitive mutant and measured growth at the restrictive temperature to assess the impact on cell division.

The experiment compared the **wild-type** (WT) strain with the **dnaA mutant** (dnaA) by measuring optical density at 600 nm (OD600) every 30 minutes over 10 hours in a plate reader. Each strain was measured in **three biological replicates**.

The data is provided in the file `data/growth_curves.tsv`.

---

### Task 1: Load and Inspect the Data

1. **Load the dataset:**
   - **Step:** Read the file `data/growth_curves.tsv` into a data frame.
   - **Instructions:**
     - Use an appropriate function for tab-delimited files.
     - Inspect the data with `head()`, `str()`, and `summary()`.

<details>
<summary>💡 Hint</summary>

```R
growth <- read.delim("data/growth_curves.tsv")
```
The function `read.delim()` reads tab-separated files by default.
</details>

---

### Task 2: Calculate Mean OD600 for Each Strain

1. **Compute the mean across replicates:**
   - **Step:** For each time point, calculate the mean OD600 across the three WT replicates and separately across the three dnaA replicates.
   - **Instructions:**
     - Create two new columns in your data frame: `WT_mean` and `dnaA_mean`.
     - Print the first few rows to verify.

<details>
<summary>💡 Hint</summary>

You can calculate a row-wise mean of specific columns using arithmetic on vectors:
```R
growth$WT_mean <- (growth$WT_rep1 + growth$WT_rep2 + growth$WT_rep3) / 3
```
Alternatively, look at the `rowMeans()` function — it can operate on a subset of columns selected with `growth[, c("col1", "col2", "col3")]`.
</details>

---

### Task 3: Plot Growth Curves

1. **Create a plot comparing both strains:**
   - **Step:** Plot time (x-axis) vs mean OD600 (y-axis) for both WT and dnaA on the same graph.
   - **Instructions:**
     - Use `plot()` for the first curve and `lines()` to add the second curve.
     - Use different colors for each strain.
     - Add axis labels (`"Time (min)"`, `"OD600"`), a title, and a legend.

<details>
<summary>💡 Hint</summary>

```R
plot(growth$time_min, growth$WT_mean, type = "b", col = "...",
     xlab = "...", ylab = "...", main = "...")
points(growth$time_min, growth$dnaA_mean, col = "...", type = "b")
legend("topleft", legend = c("WT", "dnaA"), col = c("...", "..."), lty = 1)
```
Use `type = "b"` for a line + points. You may need to set `ylim = c(0, 2)` in the `plot()` call so both curves fit.
</details>

---

### Task 4: Identify the Growth Phases

1. **Find the maximum OD600 for each strain:**
   - **Step:** Use `max()` on the mean columns to determine the maximum OD600 reached by each strain.
   - **Instructions:**
     - Print both values and comment on the difference.

2. **Estimate when exponential growth begins:**
   - **Step:** Find the first time point where the mean OD600 exceeds 0.1 for each strain.
   - **Instructions:**
     - Use logical indexing on the data frame and `min()` on the time column of the matching rows.
     - The difference between these values reflects the difference in lag phase duration.

<details>
<summary>💡 Hint</summary>

To find the first time point exceeding a threshold:
```R
wt_start <- min(growth$time_min[growth$WT_mean > 0.1])
```
Apply the same approach for `dnaA_mean`.
</details>

---

### Task 5: Plot All Replicates

1. **Visualize replicate variability:**
   - **Step:** Create a plot showing all six individual replicate curves (3 WT + 3 dnaA), each as a separate line.
   - **Instructions:**
     - Use `plot()` for the first replicate, then `lines()` for each additional replicate.
     - Use one color for all WT replicates (e.g., blue) and another for all dnaA replicates (e.g., red).
     - Add a legend indicating strain identity.

<details>
<summary>💡 Hint</summary>

Start with:
```R
plot(growth$time_min, growth$WT_rep1, type = "l", col = "blue",
     ylim = c(0, 2), xlab = "Time (min)", ylab = "OD600",
     main = "Growth Curves - All Replicates")
lines(growth$time_min, growth$WT_rep2, col = "blue")
```
Continue adding `lines()` for each remaining replicate. Setting `ylim` in the initial `plot()` ensures all curves fit.
</details>

---
