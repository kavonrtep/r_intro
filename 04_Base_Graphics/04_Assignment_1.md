## Assignment: Visualizing Chick Growth Data

### Overview

In this assignment you will work with the built-in **ChickWeight** dataset, which records the weight of chicks measured repeatedly over time while they were fed one of four diets. You will practice loading and exploring data, subsetting data frames, calculating basic statistics, and creating publication-style plots using base R graphics.

---

### Task 1: Data Exploration

Load and inspect the dataset:
- Use `data(ChickWeight)` to make the dataset available, then display the first few rows with `head()`.
- Use `str()` and `summary()` to understand its structure.
- How many chicks are in the dataset? How many diet groups?

<details>
<summary>💡 Hint</summary>

Use `length(unique(...))` on the `Chick` and `Diet` columns to count unique values.

</details>

---

### Task 2: Subsetting and Basic Statistics

1. Create a subset containing only chicks on **Diet 1**. Save it as `diet1`.
2. From `diet1`, extract only the measurements at the **last time point (day 21)**. Save it as `diet1_day21`.
3. Calculate the mean and standard deviation of chick weight at day 21 for Diet 1. Print the results.

<details>
<summary>💡 Hint</summary>

Use logical indexing `data[data$column == value, ]` or the `subset()` function to filter rows. Use `mean()` and `sd()` on the `weight` column of the filtered data.

</details>

---

### Task 3: Growth Curves — One Diet per Panel

Create a **2×2 multiplot** showing weight over time for each diet separately. Each panel should:
- Show individual chick trajectories as lines (`type = "l"`)
- Have a title indicating the diet number
- Have labelled axes

To draw multiple lines on one panel, use `plot()` for the first chick and `lines()` for the rest. Start with an empty plot (`plot(NA, xlim = ..., ylim = ...)`) so all trajectories share the same axis range.

<details>
<summary>💡 Hint</summary>

Use `par(mfrow = c(2, 2))` to set up the layout. Loop over diets with `for (d in 1:4)`, subset the data, then loop over `unique()` chick IDs within each diet. Use `range()` to set consistent axis limits across panels. Reset the layout with `par(mfrow = c(1, 1))` when done.

</details>

---

### Task 4: All Diets in One Plot — Color by Group

Create a **single scatter plot** of weight versus time using the full dataset, with points colored by diet. Add a legend.

- Define a vector of four colors, one per diet.
- Use the `Diet` column to index into the color vector to assign a color to each point.
- Add a legend with `legend()`.

<details>
<summary>💡 Hint</summary>

If `diet_colors` is a vector of four colors, then `diet_colors[ChickWeight$Diet]` gives each row its corresponding color. Pass this vector to the `col` argument of `plot()`.

</details>

---

### Task 5: Final Weight Distribution — Boxplot

Create a **boxplot** comparing final weight (day 21) across all four diets. Add:
- Axis labels and a title
- A horizontal dashed reference line at the overall mean final weight

<details>
<summary>💡 Hint</summary>

Subset to `Time == 21` first. Use `boxplot()` with a formula (`weight ~ Diet`) and the `col` argument for colors. Add the mean line with `abline(h = ...)`.

</details>

---

### Task 6: Export to PDF

Export the boxplot from Task 5 and the multiplot from Task 3 as **two separate PDF files**. Name them `final_weight_boxplot.pdf` and `growth_curves.pdf`.

<details>
<summary>💡 Hint</summary>

Wrap your plotting code between `pdf("filename.pdf")` and `dev.off()`. For the multiplot, set `par(mfrow = c(2, 2))` inside the PDF device before plotting.

</details>

---

## Submission

- Save your completed script as `your_name_chick_growth_assignment.R` and share it with me. (via email, or moodle)
