#!/usr/bin/env Rscript
# =====================================================
# Introduction to Base Graphics in R - Practice Script
# =====================================================
# This script demonstrates key examples from the presentation.
# Follow the commented instructions to explore and experiment
# with R’s base graphics functions.

# -------------------------------------
# Example 1: Scatter Plot
# -------------------------------------
set.seed(123)
x <- rnorm(50) + 1:50
y <- rnorm(50) + 1:50

# Create a basic scatter plot.
plot(x, y,
     main = "Scatter Plot",
     xlab = "X-axis",
     ylab = "Y-axis",
     pch = 19,       # solid circle
     col = "blue")

# TASK 1:
# Modify the scatter plot by:
# - Changing the point shape (try pch = 15, 17, etc.)
# - Changing the color (try "red", "green", etc.)
# Use the plot() function parameters or add new points with points().


# -------------------------------------
# Example 2: Line Plot
# -------------------------------------
x_line <- 0:10
y_line <- x_line^2

# Create a simple line plot.
plot(x_line, y_line,
     type = "l",   # line type
     main = "Line Plot",
     xlab = "X-axis",
     ylab = "Y-axis",
     col = "red")

# TASK 2:
# Create a new line plot using a different line type.
# For example, try using a dashed line (lty = 2) and add points with a different color.
# Hint: use the points() function after plot().

# -------------------------------------
# Example 3: Exploring Plot Options
# -------------------------------------
# Demonstrate different plot types by setting up a 2x2 canvas.
plot(x_line, y_line, type = "l", main = "type = 'l'")
plot(x_line, y_line, type = "p", main = "type = 'p'")
plot(x_line, y_line, type = "b", main = "type = 'b'")
plot(x_line, type = "h", main = "type = 'h'")

# TASK 3:
# Experiment with additional type options.
# For example, try using "s" (stair steps) or "o" (overplotted points and lines) in a new plot.


# -------------------------------------
# Example 4: Adding More Data with ToothGrowth Dataset
# -------------------------------------
# Load the ToothGrowth dataset.
data(ToothGrowth)
# Display the first few rows.
print(head(ToothGrowth))

# Create a simple plot: Dose vs. Tooth length.
plot(ToothGrowth$dose, ToothGrowth$len,
     main = "Tooth Growth",
     xlab = "Dose (mg/day)",
     ylab = "Tooth length",
     pch = 19,
     col = "blue")

# TASK 4:
# Create a boxplot of tooth length grouped by supplement type.
# Hint: Use the boxplot() function with the formula: len ~ supp

# -------------------------------------
# Example 5: Plotting by Supplement Type
# -------------------------------------
# Separate data by supplement type.
tg_vc <- ToothGrowth[ToothGrowth$supp == "VC", ]
tg_oj <- ToothGrowth[ToothGrowth$supp == "OJ", ]

# Plot data for VC.
plot(tg_vc$dose, tg_vc$len,
     main = "Tooth Growth by Supplement Type",
     xlab = "Dose (mg/day)",
     ylab = "Tooth length",
     pch = 19,
     col = "blue")
# Add points for OJ.
points(tg_oj$dose, tg_oj$len,
       pch = 19,
       col = "red")
# Add a legend.
legend("bottomright", legend = c("VC", "OJ"), col = c("blue", "red"), pch = 19)

# TASK 5:
# Modify the above plot by adding a regression line (or a smoothing line)
# for one of the supplement types.
# Hint: use lm() to fit a simple model and abline() to add the line.

# -------------------------------------
# Example 6: Using plot(), points(), and lines() Functions
# -------------------------------------
# Create a new plot with multiple graphical elements.
par(mfrow = c(1, 3))
xx <- 1:10
yy <- c(1, 3, 2, 5, 4, 7, 6, 9, 8, 10)

# Single plot using plot().
plot(xx, yy, type = "b", main = "Single Plot")

# Plot with additional points.
plot(xx, yy, type = "b", main = "Plot + Points")
points(xx, yy/2, type = "p", col = "red")

# Plot with additional points and lines.
plot(xx, yy, type = "b", main = "Plot + Points + Lines")
points(xx, yy/2, type = "l", col = "red")
lines(xx, yy*2, col = "blue")
# Reset layout.
par(mfrow = c(1, 1))

# TASK 6:
# On the third plot, try adding:
# - New points with different shapes/colors.
# - Annotations using the text() function to label some points.
# - A custom title using title() or mtext().
# Experiment with layering different elements.




################################################################################
# R SCRIPT: Analysis and Plotting of ChIP-seq Bedgraph Data for chr1
#
# INSTRUCTIONS FOR STUDENTS:
#
# Your task is to create a script that:
#   1. Loads a bedgraph file (without a header) into a data frame.
#      This bedgraph file contains ChIP-seq data with columns:
#      "chr_name", "start", "end", "score".
#      Score is a numeric value and represent enrichment obtained from ChIP-seq experiment.
#   2. Assigns the correct column names.
#   3. Subsets the data for rows corresponding to chromosome "chr1".
#   4. Creates a color vector where most points are "black" but changes to "red"
#      for rows with a score above a specified threshold (e.g., score > 5).
#   5. Plots the data with base R graphics using a scatter/histogram-like plot.
#   6. Adds a legend to explain the color coding.
#
# HINTS:
#
# 1. Loading the Data:
#    - Use the read.table() function.
#    - Since the file is tab-delimited and has no header, set header = FALSE
#      and sep = "\t".
#    - Example: my_data <- read.table("your_file.bedgraph", header = FALSE, sep = "\t")
#
# 2. Defining Column Names:
#    - Use the colnames() function to set names to: "chr_name", "start", "end", "score".
#    - Example: colnames(my_data) <- c("chr_name", "start", "end", "score")
#
# 3. Subsetting the Data:
#    - Extract only the rows where the chromosome name equals "chr1".
#    - You can use subset() or indexing (e.g., my_data[my_data$chr_name == "chr1", ])
#
# 4. Creating a Color Vector:
#    - Create a vector of colors (same length as the number of rows in your subset).
#    - Use rep() to fill the vector with "black" initially.
#    - Then, update the vector so that for any row where score > 5, the color becomes "red".
#    - Hint: Use conditional indexing, e.g., colors[subset_data$score > 5] <- "red"
#
# 5. Plotting the Data:
#    - Use the plot() function to create a scatter plot.
#    - Set the x-axis to the 'start' column and y-axis to the 'score' column.
#    - Use type = "h" for a plot with vertical lines (histogram-like look).
#    - Add labels using xlab, ylab, and main. You can include newline characters (e.g., "\n")
#      in the main title if needed.
#
# 6. Adding a Legend:
#    - Use the legend() function.
#    - Place it in the "topright" of the plot.
#    - The legend should describe the two categories:
#         "Score <= 5" (black) and "Score > 5" (red).
#
# 7. Experiment:
#    - After you complete the above, try modifying the threshold (e.g., change 5 to another value)
#      or adjust the plot parameters to better visualize your data.
#
# NOTES:
#   - Ensure your bedgraph file (e.g., "chip_seq.bedgraph") is saved in your working directory.
#   - Check your data by using functions such as head(), tail(), and str() before plotting.
#
# Write your code below by following these steps interactively.
################################################################################

# Step 1: Load the bedgraph data into a data frame.
# Hint: Use read.table() with header = FALSE and sep = "\t".
# path to the bedgraph file is "data/chip_seq.bedgraph"


# Step 2: Define the column names.
# Hint: Use colnames() to assign names: "chr_name", "start", "end", "score".


# Step 3: Inspect your data to ensure it's loaded correctly.
# Hint: Use head(), tail(), and str() to check the data structure.

# Step 4: Subset the data to include only rows for chromosome "chr1".
# Hint: Use subset() or indexing [ ] to filter the data.


# Step 5: Plot the data.
# Hint: Use plot() to create a scatter plot of 'start' vs. 'score'.
#       - test different values for pch, col, and type ('l', 'p', 'h' parameters.
#       - Label your axes with xlab and ylab.


# Step 6: Create a color vector based on the 'score' values.
# Hint: use ifelse() function to assign colors based on a condition.
#       - For example, colors <- ifelse(subset_data$score > 5, "red", "black")
#       - Use the colors vector in the plot() function to color the points.

# Once you have written your code, run it interactively to see your plot and make adjustments.



# =====================================================
# Introduction to Histograms in R - Practice Script
# =====================================================
# This script demonstrates key examples for working with histograms.
# Follow the commented instructions to explore and experiment with
# the hist() function and its parameters in base R.

# -------------------------------------
# Example 1: Simple Histogram
# -------------------------------------
# Using the built-in mtcars dataset, create a simple histogram of the 'mpg' variable.
hist(mtcars$mpg,
     main = "Histogram of Miles Per Gallon",
     xlab = "Miles Per Gallon",
     col = "lightblue")

# TASK 1:
# Modify the above histogram by:
# - Changing the title and axis labels.
# - Changing the color of the bars (try "pink", "yellow", etc.).
# - Adding additional parameters such as 'xlim' to adjust the x-axis limits.
# Experiment with these options and observe how the histogram appearance changes.


# -------------------------------------
# Example 2: Adjusting Histogram Breaks
# -------------------------------------
# Create a histogram of the 'mpg' variable with a specific number of breaks.
hist(mtcars$mpg,
     breaks = 40,  # Try modifying this value (e.g., 5 or 15) to see the effect on binning.
     main = "Histogram of MPG with 10 Breaks",
     xlab = "Miles Per Gallon",
     col = "lightgreen")

# TASK 2:
# Modify the above histogram to use a custom vector of breakpoints instead of a single number.
# Hints:
#   - Create a break vector using seq(), for example:
#         breaks_vec <- seq(min(mtcars$mpg), max(mtcars$mpg), length.out = 12)
#   - Pass this vector to the breaks argument.
# Experiment with different numbers of breakpoints and compare the resulting histograms.


# -------------------------------------
# Example 3: Overlaying Histograms for Two Datasets
# -------------------------------------
# Generate two sample datasets.
set.seed(123)
data1 <- rnorm(100, mean = 0, sd = 1)
data2 <- rnorm(100, mean = 1, sd = 1)

# Define common breaks for both histograms. Brakes must include the range of both datasets.
common_breaks <- seq(min(c(data1, data2)),
                     max(c(data1, data2)),
                     length.out = 15)

# Create the histogram for data1.
# We set 'xlim' and 'ylim' so that both histograms can be properly overlaid.
h1 <- hist(data1,
           breaks = common_breaks,
           main = "Overlayed Histograms",
           xlab = "Value",
           col = "#FFAAAAAA",  # semi-transparent
           border = "white",
           xlim = range(common_breaks),
)

# Overlay the histogram for data2.
hist(data2,
     breaks = common_breaks,
     col = "#AAFFAAAA",  # semi-transparent
     border = "white",
     add = TRUE)

# Add a legend to differentiate the two datasets.
legend("topright", legend = c("Data 1", "Data 2"),
       fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))

# TASK 3:
# Modify the overlayed histograms by:
# - Adjusting the transparency of the colors.
# - Changing the common break values (try a different number of breakpoints).
# - Experiment by altering the mean or standard deviation of data1 or data2, then overlay histograms.
# - Explore other parameters such as 'freq' (for frequency) vs. 'prob' (for density) by setting freq = FALSE.
# Run your modified code and observe how the histogram visualizations change.



################################################################################
# Histogram of Gene Locations - Practice Script
################################################################################
# This script guides you through:
#   1. Importing and preparing gene annotation data from a BED file.
#   2. Exploring the data by inspecting the number of genes per chromosome.
#   3. Creating a histogram of gene start positions for a single chromosome.
#   4. Creating a multiplot layout of histograms for chromosomes chr1 through chr6.
#
# STUDENT TASKS:
#   - Follow the instructions provided in the comments.
#   - Modify histogram parameters (breaks, colors, axis limits) to explore their effects.
#   - Experiment with different layout options using par(mfrow).
#
# NOTES:
#   - Make sure the "genes.bed" file is in the "data" folder within your working directory.
#   - This exercise assumes familiarity with data import, subsetting, and base R plotting.
################################################################################

# -----------------------------------------------------
# Step 1: Import and Prepare the Data
# -----------------------------------------------------
# Import the BED file. BED files are tab-delimited and typically have no header.
# use read.table function, data are located in "data/genes.bed", there is no header in the file
genes <- read.table("data/genes.bed", header = FALSE)

# Inspect the first few rows to understand the structure.

# By default, the column names are V1, V2, ..., V6.
# Assign meaningful column names based on the standard BED format.
colnames(genes) <- c("chrom", "start", "end", "name", "score", "strand")
head(genes)

# Alternatively, you could define the header while reading the file:
# genes <- read.table("data/genes.bed", header = FALSE,
#                     col.names = c("chrom", "start", "end", "name", "score", "strand"))

# -----------------------------------------------------
# Step 2: Explore the Data
# -----------------------------------------------------
# Check the number of genes in the file.
n_genes <- nrow(genes)
print(n_genes)

# Determine the number of genes per chromosome.
n_genes_chr <- table(genes$chrom)
print(n_genes_chr)

# Create a simple plot to visualize the number of genes per chromosome.
plot(n_genes_chr,
     main = "Number of Genes per Chromosome",
     xlab = "Chromosome",
     ylab = "Number of Genes")
# try to use barplot instead of plot

# Calculate gene lengths and add them as a new column.
genes$length <- genes$end - genes$start

# TASK: Create histogram of gene lengths



# -----------------------------------------------------
# Step 3: Histogram for a Single Chromosome
# -----------------------------------------------------
# TASK: Create a histogram of gene start positions for a single chromosome.
# Hints:
#   - Use subsetting to extract genes from a specific chromosome (e.g., "chr12").
#   - Use the hist() function to plot the distribution of the 'start' positions.
#   - Customize the histogram by adjusting the number of breaks, colors, and labels.
#
# Example Code:
chr12_genes <- genes[genes$chrom == "chr12", ]

# TASK: use hist function to plot the distribution of the 'start' positions


# TASK:
# - Experiment by changing the number of breaks (try 100 or 300).
# - Modify the color and add xlim or ylim parameters to better visualize the data.
# - Consider adding additional labels or gridlines using abline() or grid().

# -----------------------------------------------------
# Step 4: Multiplot of Histograms for Chromosomes chr1 - chr6
# -----------------------------------------------------
# TASK: Create a multiplot layout to display histograms of gene start positions for chromosomes chr1 to chr6.
#
# Hints:
#   - Use par(mfrow = c(2, 3)) to set up a layout with 2 rows and 3 columns.
#   - Loop through a vector of chromosome names (e.g., paste0("chr", 1:6)).
#   - Define a common breaks vector for all histograms to ensure consistency.
#
# Example Code:
# Set up a multiplot layout.
par(mfrow = c(2, 3))
# Define the chromosomes to plot.
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6")

# Define common breaks based on the overall range of gene start positions.
common_breaks <- seq(min(genes$start), max(genes$start), length.out = 50)

# Loop through each chromosome and create its histogram.
# Use a for loop to iterate over the chromosome names.
# select subset of genes for each chromosome using subset function or indexing
# use hist function to plot the distribution of the 'start' positions


# TASK:
# - Modify the layout dimensions (e.g., try par(mfrow = c(3, 2))).
# - Change the color or transparency in the histograms.
# - Experiment with different numbers of common breaks to see how they affect each histogram.
# - Optionally, add legends or annotations to your multiplot.
