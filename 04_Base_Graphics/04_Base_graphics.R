################################################################################
# R SCRIPT FOR SESSION 4: INTRODUCTION TO BASE GRAPHICS IN R
#
# This script covers:
# 1. Core plot types: scatter, line, histogram, boxplot, barplot, stripchart
# 2. Customizing plots: colors, shapes, sizes, axis limits, log scale
# 3. Layering: points(), lines(), legend(), annotations
# 4. Multi-plot layouts: par(mfrow), layout()
# 5. Exporting plots: png(), pdf()
#
# DATA SOURCES USED:
# - ToothGrowth: built-in dataset (guinea pig tooth growth experiment)
# - mtcars: built-in dataset (car performance data)
# - data/genes.bed: BED file with gene annotations
# - data/chip_seq.bedgraph: ChIP-seq bedgraph data
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################

# ---------- SLIDE: Base Graphics ----------

####################### 1. SCATTER AND LINE PLOTS ##############################

# ---------- SLIDE: Scatter Plot ----------

set.seed(123)
x <- rnorm(50) + 1:50
y <- rnorm(50) + 1:50

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

# ---------- SLIDE: Line Plot ----------

x_line <- 0:10
y_line <- x_line^2

plot(x_line, y_line,
     type = "l",
     main = "Line Plot",
     xlab = "X-axis",
     ylab = "Y-axis",
     col = "red")

# TASK 2:
# Create a new line plot using a different line type.
# Try a dashed line (lty = 2) and add points with a different color.
# Hint: use the points() function after plot().

# ---------- SLIDE: Plot types: `type` argument ----------

par(mfrow = c(2, 2))
plot(x_line, y_line, type = "l", main = "type = 'l'")
plot(x_line, y_line, type = "p", main = "type = 'p'")
plot(x_line, y_line, type = "b", main = "type = 'b'")
plot(x_line,         type = "h", main = "type = 'h'")
par(mfrow = c(1, 1))

# TASK 3:
# Experiment with additional type options.
# Try "s" (stair steps) or "o" (overplotted points and lines).

####################### 2. WORKING WITH DATASETS ##############################

# ---------- SLIDE: Adding data: ToothGrowth dataset ----------

data(ToothGrowth)
head(ToothGrowth)

# ---------- SLIDE: Adding data: simple plot ----------

plot(ToothGrowth$dose, ToothGrowth$len,
     main = "Tooth Growth",
     xlab = "Dose (mg/day)",
     ylab = "Tooth length",
     pch = 19,
     col = "blue")

# TASK 4:
# Create a boxplot of tooth length grouped by supplement type.
# Hint: Use boxplot() with the formula: len ~ supp

# ---------- SLIDE: Adding data: plotting by group ----------

tg_vc <- ToothGrowth[ToothGrowth$supp == "VC", ]
tg_oj <- ToothGrowth[ToothGrowth$supp == "OJ", ]

plot(tg_vc$dose, tg_vc$len,
     main = "Tooth Growth by Supplement Type",
     xlab = "Dose (mg/day)",
     ylab = "Tooth length",
     pch = 19,
     col = "blue")
points(tg_oj$dose, tg_oj$len, pch = 19, col = "red")
legend("bottomright", legend = c("VC", "OJ"),
       col = c("blue", "red"), pch = 19)

# TASK 5:
# Add a regression line for one of the supplement types.
# Hint: use lm() to fit a model and abline() to draw the line.

# ---------- SLIDE: `plot`, `points`, and `lines` ----------

par(mfrow = c(1, 3))
xx <- 1:10
yy <- c(1, 3, 2, 5, 4, 7, 6, 9, 8, 10)

plot(xx, yy, type = "b", main = "single plot")

plot(xx, yy, type = "b", main = "plot + points")
points(xx, yy / 2, type = "p", col = "red")

plot(xx, yy, type = "b", main = "plot + points + lines")
points(xx, yy / 2, type = "l", col = "red")
lines(xx, yy * 2, col = "blue")

par(mfrow = c(1, 1))

# TASK 6:
# On the third plot, add:
# - New points with different shapes/colors.
# - Annotations using text() to label some points.
# - A custom title using title() or mtext().

####################### 3. CATEGORICAL DATA ###################################

# ---------- SLIDE: Plotting categorical data: `boxplot` ----------

# len ~ supp is a formula: plot len on y-axis grouped by supp
boxplot(len ~ supp, data = ToothGrowth,
        main = "Tooth Growth by Supplement Type",
        xlab = "Supplement Type",
        ylab = "Tooth length")

# ---------- SLIDE: `boxplot` input formats ----------

# boxplot also accepts a list
example_list <- list(
  a = rnorm(100),
  b = rnorm(100, 3, 0.5),
  c = rnorm(100, -2, 5)
)
boxplot(example_list)

# ---------- SLIDE: Categorical data: `stripchart` ----------

stripchart(len ~ supp, data = ToothGrowth,
           main = "Tooth Growth by Supplement Type",
           xlab = "Supplement Type",
           ylab = "Tooth length",
           col = "blue",
           vertical = TRUE)

# ---------- SLIDE: `stripchart` with jitter ----------

stripchart(len ~ supp, data = ToothGrowth,
           main = "Tooth Growth by Supplement Type",
           xlab = "Supplement Type",
           ylab = "Tooth length",
           col = "blue",
           method = "jitter",
           vertical = TRUE)

# ---------- SLIDE: Bar plot vs. Box plot ----------

par(mfrow = c(1, 2))
boxplot(mpg ~ cyl, data = mtcars,
        main = "Boxplot of MPG by Cylinders",
        xlab = "Number of Cylinders",
        ylab = "Miles per Gallon",
        col = "lightblue")
barplot(table(mtcars$cyl),
        main = "Bar Chart of Cylinders",
        xlab = "Number of Cylinders",
        ylab = "Frequency",
        col = "lightgreen")
par(mfrow = c(1, 1))

# ---------- SLIDE: Bar plot example ----------

par(mfrow = c(1, 2))
barplot(table(mtcars$gear),
        main = "Bar Chart of Gears",
        xlab = "Number of Gears",
        ylab = "Frequency",
        col = "lightgreen")

counts <- table(mtcars$am, mtcars$gear)
barplot(counts,
        main = "Gears by Transmission",
        xlab = "Number of Gears",
        ylab = "Frequency",
        col = c("lightblue", "lightgreen"))
par(mfrow = c(1, 1))

####################### 4. HISTOGRAMS ##########################################

# ---------- SLIDE: Histogram ----------

x2 <- rnbinom(100, mu = 10, size = 10)
hist(x2, main = "Histogram of Negative Binomial Distribution")

# ---------- SLIDE: Histogram: controlling breaks ----------

par(mfrow = c(1, 2))
hist(x2, main = "breaks = 20", breaks = 20)
hist(x2, main = "breaks = 0:40", breaks = 0:40)
par(mfrow = c(1, 1))

# TASK 7:
# Experiment with breaks:
# - Use a custom vector from seq(), e.g.:
#     breaks_vec <- seq(min(x2), max(x2), length.out = 12)
# - Try freq = FALSE to show density instead of counts.
# - Overlay two histograms using add = TRUE.

####################### 5. CUSTOMIZING PLOTS ##################################

# ---------- SLIDE: Customizing plots: colors, shapes, size ----------

par(mfrow = c(2, 3))
plot(x, y, main = "default plot")
plot(x, y, col = "#FF0000",   main = "col sets point color")
plot(1:20, col = rainbow(20), main = "col as a color vector")
plot(1:20, pch = 1:20,        main = "pch sets point shape")
plot(x, y, cex = 3,           main = "cex sets point size")
plot(x, y, main = "lines and grid")
abline(0, 2)
grid()
par(mfrow = c(1, 1))

# ---------- SLIDE: Customizing plots: axis limits and log scale ----------

par(mfrow = c(2, 2))
plot(x, y, xlim = c(0, 30), ylim = c(0, 20), type = "b",
     main = "setting axis limits")
plot(x, y, log = "y", type = "b", main = "log scale on y-axis")
plot(x, y, log = "x", type = "b", main = "log scale on x-axis")
plot(x, y, log = "xy", type = "b", main = "log scale on both axes")
par(mfrow = c(1, 1))

# ---------- SLIDE: Adding legends and annotations ----------

plot(x, y, pch = 19, col = "blue")
points(x + 0.5, y + 0.5, col = "red", pch = 17)
legend("topleft",
       legend = c("Group 1", "Group 2"),
       col = c("blue", "red"),
       pch = c(19, 17))
mtext(expression(alpha + beta == gamma), side = 3, cex = 2)

####################### 6. MULTI-PLOT LAYOUTS #################################

# ---------- SLIDE: Combining plots: `par(mfrow)` ----------

par(mfrow = c(2, 3))
plot(x, y, main = "plot 1", col = "1")
plot(x, y, main = "plot 2", col = "2")
plot(x, y, main = "plot 3", col = "3")
plot(x, y, main = "plot 4", col = "4")
plot(x, y, main = "plot 5", col = "5")
plot(x, y, main = "plot 6", col = "6")
par(mfrow = c(1, 1))

# ---------- SLIDE: Combining plots: `layout` ----------

layout(matrix(c(1, 2, 4, 3), nrow = 2),
       width = c(2, 1), height = c(1, 2))
layout.show(4)

# ---------- SLIDE: Combining plots: `layout` example ----------

layout(matrix(c(1, 2, 4, 3), nrow = 2),
       width = c(2, 1), height = c(1, 2))
hist(x2, main = "1st plot")
plot(x, y, main = "2nd plot")
hist(rnorm(100), main = "3rd plot")

####################### 7. EXPORTING PLOTS ####################################

# ---------- SLIDE: Exporting plots: bitmap formats ----------

# Create output directory if needed
dir.create("outputs", showWarnings = FALSE)

png("outputs/plot1.png", width = 800, height = 600)
plot(x, y, main = "Scatter Plot", col = "blue")
dev.off()

# ---------- SLIDE: Exporting plots: vector formats ----------

# pdf saves each plot() call as a separate page
pdf("outputs/plot2.pdf", width = 8, height = 6)
plot(x, y, main = "Scatter Plot", col = "blue")   # page 1
hist(x2, main = "Histogram")                       # page 2
plot(x, y, main = "Line Plot", col = "red")        # page 3
dev.off()

####################### 8. REAL-DATA EXERCISES ################################

# ---------- SLIDE: `par` function ----------

# par() sets graphical parameters that apply to subsequent plots.
# Examples:
#   par(mfrow = c(2, 2))           multi-plot layout
#   par(mar = c(5, 4, 4, 2) + 0.1) plot margins
#   par(cex = 1.5)                  text/symbol magnification
#   par(lty = 2)                    dashed lines
# See ?par for the full list.

# ---------- SLIDE: Summary ----------

################################################################################
# EXERCISE: ChIP-seq Bedgraph Data
#
# Your task: create a script that
#   1. Loads "data/chip_seq.bedgraph" (tab-delimited, no header) with columns:
#      chr_name, start, end, score
#   2. Subsets rows for chromosome "chr1".
#   3. Creates a color vector: "black" for score <= 5, "red" for score > 5.
#   4. Plots start vs. score using type = "h".
#   5. Adds a legend.
#
# HINTS:
#   my_data <- read.table("data/chip_seq.bedgraph", header = FALSE, sep = "\t")
#   colnames(my_data) <- c("chr_name", "start", "end", "score")
#   chr1 <- my_data[my_data$chr_name == "chr1", ]
#   colors <- ifelse(chr1$score > 5, "red", "black")
################################################################################

# Step 1: Load the bedgraph data
# (your code here)

# Step 2: Assign column names
# (your code here)

# Step 3: Inspect the data (head, str)
# (your code here)

# Step 4: Subset to chr1
# (your code here)

# Step 5: Plot start vs. score
# (your code here)

# Step 6: Color by score threshold and re-plot
# (your code here)

# Step 7: How many chromosomes are in the dataset? Use table() to find out.

# Create multipanel plots for first 4 chromosomes (chr1 - chr4) using par(mfrow).
# use a loop to subset each chromosome and plot start vs. score with type = "h".
# (your code here)

################################################################################
# EXERCISE: Gene Location Histograms (BED file)
#
# Your task:
#   1. Load "data/genes.bed" (tab-delimited, no header, 6 columns:
#      chrom, start, end, name, score, strand).
#   2. Explore the number of genes per chromosome (table()).
#   3. Create a histogram of gene start positions for "chr12".
#   4. Create a 2x3 multiplot of histograms for chr1 through chr6.
################################################################################

# Step 1: Load genes.bed
genes <- read.table("data/genes.bed", header = FALSE,
                    col.names = c("chrom", "start", "end",
                                  "name", "score", "strand"))

# Step 2: Explore - number of genes per chromosome
n_genes_chr <- table(genes$chrom)
print(n_genes_chr)
barplot(n_genes_chr,
        main = "Number of Genes per Chromosome",
        xlab = "Chromosome", ylab = "Number of Genes")

# Calculate gene lengths
genes$length <- genes$end - genes$start
# TASK: Create a histogram of gene lengths

# Step 3: Histogram for chr12
chr12_genes <- genes[genes$chrom == "chr12", ]
# TASK: use hist() to plot the distribution of start positions
# Experiment with breaks = 100 or 300, colors, xlim.

# Step 4: Multiplot for chr1 - chr6
par(mfrow = c(2, 3))
chromosomes <- paste0("chr", 1:6)
common_breaks <- seq(min(genes$start), max(genes$start), length.out = 50)

for (chr in chromosomes) {
  chr_genes <- genes[genes$chrom == chr, ]
  # TASK: add hist() call here using chr_genes$start and common_breaks
}
par(mfrow = c(1, 1))