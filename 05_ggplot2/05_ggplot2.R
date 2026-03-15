################################################################################
# R SCRIPT FOR SESSION 5: INTRODUCTION TO ggplot2 GRAPHICS IN R
#
# This script covers:
# 1. Data, mapping, and layers — the ggplot2 grammar
# 2. Mapping aesthetics (color, size, shape)
# 3. Geoms and statistical transformations (geom_smooth)
# 4. Scales — color palettes, axis transformations, limits
# 5. Facets — facet_wrap() and facet_grid()
# 6. Coordinates and themes
# 7. Other plot types: bar, line, boxplot, histogram, density
# 8. Combining plots with gridExtra
# 9. Exporting with ggsave()
# Final exercise: gene density along chromosomes (BED file)
#
# DATA SOURCES USED:
# - mpg: built-in ggplot2 dataset (fuel economy data)
# - iris: built-in dataset (flower measurements)
# - mtcars: built-in dataset (car performance data)
# - Orange: built-in dataset (orange tree growth)
# - data/genes.bed: BED file with gene annotations
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################

library(ggplot2)
library(gridExtra)

# ---------- SLIDE: ggplot2 ----------

################################################################################
# SECTION 1: DATA, MAPPING, AND LAYERS
################################################################################

# ---------- SLIDE: Data ----------

# ggplot() takes the dataset as its first argument.
# Nothing is drawn yet — we need mapping and a geom.
ggplot(data = mpg)

# ---------- SLIDE: Mapping ----------

# aes() maps dataset columns to visual properties.
# Still no geom — still nothing drawn.
ggplot(mpg, aes(x = cty, y = hwy))

# ---------- SLIDE: Example: Scatter Plot ----------

ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point() +
  ggtitle("Scatter Plot: City vs Highway MPG") +
  xlab("City MPG") +
  ylab("Highway MPG")

# TASK 1:
# Modify the scatter plot above to:
#   - Change the point color to "darkgreen" (fixed, not mapped).
#   - Increase the point size.
# Hint: add color = ... and size = ... inside geom_point().


################################################################################
# SECTION 2: MAPPING AESTHETICS
################################################################################

# ---------- SLIDE: Mapping color to a variable ----------

ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  ggtitle("Scatter Plot: Color by Car Class") +
  xlab("City MPG") +
  ylab("Highway MPG")

# ---------- SLIDE: Fixed vs. mapped aesthetics ----------

# Mapped (inside aes): color varies by data
ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point()

# Fixed (outside aes): all points the same color
ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point(color = "darkgreen")

# TASK 2:
# Modify the color-by-class plot to:
#   - Also map 'displ' (engine displacement) to point size.
#   - Change the color scale with scale_color_brewer(palette = "Set1").
# Note: scale_color_brewer() works for discrete variables.


################################################################################
# SECTION 3: LAYERS — GEOMS AND STATISTICAL TRANSFORMATIONS
################################################################################

# ---------- SLIDE: Adding multiple layers ----------

ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  ggtitle("Scatter Plot with Linear Trend Line") +
  xlab("City MPG") +
  ylab("Highway MPG")

# TASK 3:
# Change the trendline method from "lm" to "loess" in geom_smooth().
# LOESS = locally estimated scatterplot smoothing — a flexible curve.
# Also try se = FALSE to remove the confidence band.


################################################################################
# SECTION 4: SCALES
################################################################################

# ---------- SLIDE: Example: Viridis color scale ----------

g1 <- ggplot(mpg, aes(cty, hwy, colour = class)) +
  geom_point() + ggtitle("Default color scale")

g2 <- ggplot(mpg, aes(cty, hwy, colour = class)) +
  geom_point() +
  scale_colour_viridis_d() + ggtitle("Viridis color scale")

grid.arrange(g1, g2, ncol = 2)

# ---------- SLIDE: Example: Log-scaled axis ----------

g1 <- ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point() + ggtitle("Linear scale")

g2 <- ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point() +
  scale_x_log10() + ggtitle("Log10 x-axis")

grid.arrange(g1, g2, ncol = 2)

# TASK 4:
# Modify one of the above plots to set custom axis limits:
#   - Use scale_x_continuous(limits = c(10, 30)) to restrict the x-axis.
#   - Use scale_y_continuous(limits = c(20, 45)) to restrict the y-axis.


################################################################################
# SECTION 5: FACETS
################################################################################

# ---------- SLIDE: facet_wrap example ----------

ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3) +
  facet_wrap(~ Species) +
  ggtitle("Faceted Scatter Plot: Iris Data") +
  xlab("Sepal Length") +
  ylab("Sepal Width")

# TASK 5:
# What happens when you use facet_wrap(~ Species, ncol = 1)?
# Try also nrow = 1. How does this change the layout?

# ---------- SLIDE: facet_grid example ----------

# facet_grid() creates a strict rows × columns grid using two variables.
ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  facet_grid(drv ~ cyl) +
  ggtitle("MPG by Drive Type and Cylinder Count") +
  xlab("City MPG") +
  ylab("Highway MPG")

# TASK 6:
# The facet_grid plot above uses fixed (shared) scales across all panels.
# Add scales = "free" inside facet_grid() — what changes?
# Also try scales = "free_x" and scales = "free_y" separately.
# When would free scales be misleading? When are they useful?


################################################################################
# SECTION 6: COORDINATES AND THEMES
################################################################################

# ---------- SLIDE: Coordinates ----------

ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Box Plot: MPG by Cylinder Count (Flipped)")

# TASK 7:
# Zoom into the boxplot above using coord_cartesian(ylim = c(15, 25)).
# Then try the same zoom with scale_y_continuous(limits = c(15, 25)).
# Compare the two results — how does the boxplot shape differ?
# Hint: scale_*() removes data outside the limits before computing stats;
#       coord_cartesian() only clips the viewport, keeping all data.

# ---------- SLIDE: Customizing with theme() ----------

ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.line = element_line(size = 0.75, color = "gray"),
    text = element_text(size = 12)
  ) +
  ggtitle("Scatter Plot with Custom Theme")

# TASK 8:
# Modify one of the above plots to:
#   - Switch to theme_dark().
#   - Adjust title and axis label font size using:
#     theme(text = element_text(size = 14, family = "Arial"))


################################################################################
# SECTION 7: OTHER PLOT TYPES
################################################################################

# ---------- SLIDE: Bar plot ----------

ggplot(mpg, aes(x = class)) +
  geom_bar(fill = "steelblue") +
  ggtitle("Bar Plot: Count of Car Classes") +
  xlab("Car Class") +
  ylab("Count")

# Bar plot with pre-computed counts (stat = "identity"):
x <- c(Mazda = 10, Toyota = 20, Honda = 30)
ggplot(data.frame(x = names(x), y = x), aes(x, y)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Brand") + ylab("Count")

# TASK 9:
# Map 'drv' (drive train: f/r/4) to the fill aesthetic in a bar plot of 'class':
#   ggplot(mpg, aes(x = class, fill = drv)) + geom_bar()
# Then switch to position = "fill" to show proportions instead of counts.
# Bonus: reorder the x-axis bars by total count:
#   aes(x = reorder(class, class, length), fill = drv)

# ---------- SLIDE: Line plot ----------

df_line <- data.frame(x = 1:10, y = (1:10)^2)
ggplot(df_line, aes(x = x, y = y)) +
  geom_line(color = "red", size = 1.2) +
  ggtitle("Line Plot: y = x²") +
  xlab("X") +
  ylab("Y")

# TASK 10:
# Add geom_point() to the line plot above to mark each data point on the line.

# ---------- SLIDE: Line plot — multiple lines ----------

# Dataset: Orange
# Built-in R data frame — 35 rows, 3 columns.
# Records of orange tree growth:
#   Tree          — ordered factor, tree ID (ordered by increasing max diameter)
#   age           — numeric, days since 1968-12-31
#   circumference — numeric, trunk circumference at breast height (mm)
# 5 trees measured at 7 time points each.

ggplot(Orange, aes(x = age, y = circumference, color = Tree)) +
  geom_line() +
  ggtitle("Orange Tree Growth") +
  xlab("Age (days)") + ylab("Circumference (mm)")

# ---------- SLIDE: Box plot ----------

ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(fill = "orange") +
  ggtitle("Box Plot: MPG by Cylinder Count") +
  xlab("Cylinders") +
  ylab("MPG")

# TASK 11:
# Replace geom_boxplot() with geom_violin(fill = "lightblue") on the same data.
# Then add a layer of raw points on top:
#   geom_jitter(width = 0.1, alpha = 0.4)
# What does the violin shape show that the boxplot does not?

# ---------- SLIDE: Histogram ----------

ggplot(iris, aes(x = Sepal.Length)) +
  geom_histogram(bins = 20, fill = "purple", alpha = 0.7) +
  ggtitle("Histogram of Sepal Length") +
  xlab("Sepal Length") +
  ylab("Frequency")

# ---------- SLIDE: Histogram — overlapping by species ----------

ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_histogram(alpha = 0.5, bins = 20, position = "identity") +
  ggtitle("Sepal Length by Species") +
  xlab("Sepal Length") + ylab("Count")

# TASK 12:
# Try different position settings in the histogram above: "stack", "dodge", "fill".
# Also try facet_wrap(~ Species) instead of using position.

# ---------- SLIDE: Histogram — faceted ----------

ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_histogram(alpha = 0.5, bins = 20) +
  facet_grid(Species ~ .) +
  ggtitle("Sepal Length — Faceted by Species")

# ---------- SLIDE: Density plot ----------

ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_density(alpha = 0.5) +
  ggtitle("Density Plot of Sepal Length by Species") +
  xlab("Sepal Length") +
  ylab("Density")

# TASK 13:
# In the density plot above, swap fill = Species for color = Species.
# How does the appearance change (filled area vs. outline only)?
# Now add fill = "grey80" as a fixed aesthetic inside geom_density() —
# what does it do when combined with the color mapping?

# TASK 14:
# Choose one plot type from this section (bar, line, box, histogram, density).
# Customize it: change colors, labels, alpha, add a theme, or use a different dataset.


################################################################################
# SECTION 8: COMBINING PLOTS
################################################################################

# ---------- SLIDE: grid.arrange ----------

p_density <- ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_density(alpha = 0.5) +
  facet_grid(Species ~ .)

p_line <- ggplot(Orange, aes(x = age, y = circumference, color = Tree)) +
  geom_line()

p_box <- ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(fill = "orange")

grid.arrange(p_density, p_line, p_box, ncol = 2)

# ---------- SLIDE: Custom layout with layout_matrix ----------

grid.arrange(p_density, p_line, p_box,
             layout_matrix = rbind(c(1, 1), c(2, 3)),
             heights = c(2, 1))

# TASK 15:
# Experiment with different layouts:
#   - Change nrow or ncol in grid.arrange().
#   - Try a different layout_matrix.
#   - Add one of your plots from earlier tasks.


################################################################################
# SECTION 9: EXPORTING PLOTS
################################################################################

# ---------- SLIDE: ggsave ----------

# ggsave() saves the last plot, or a named plot object.
combined_plot <- grid.arrange(p_density, p_line, p_box,
                              layout_matrix = rbind(c(1, 1), c(2, 3)),
                              heights = c(2, 1))

ggsave("combined_plot.png", combined_plot, width = 10, height = 8, dpi = 300)
ggsave("combined_plot.pdf", combined_plot, width = 10, height = 8)

# TASK 16:
# Export one of your individual plots to a PDF file.
# Try different formats (jpg, tiff) and adjust width/height.


################################################################################
# FINAL EXERCISE: GENE DENSITY ALONG CHROMOSOMES
################################################################################

# ---------- SLIDE: Exercise overview ----------

# Step 1: Import and prepare the data
genes <- read.table("data/genes.bed", header = FALSE)
colnames(genes) <- c("chrom", "start", "end", "name", "score", "strand")
head(genes)

# ---------- SLIDE: Step 2: Filter chromosomes ----------

genes_filtered <- subset(genes,
  chrom %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6"))
cat("Number of genes in chr1-chr6:", nrow(genes_filtered), "\n")

# ---------- SLIDE: Step 3: Bar plot of gene counts ----------

ggplot(genes_filtered, aes(x = chrom, fill = chrom)) +
  geom_bar() +
  ggtitle("Number of Genes per Chromosome (chr1 - chr6)") +
  xlab("Chromosome") + ylab("Number of Genes") +
  theme_minimal()

# ---------- SLIDE: Step 4: Faceted gene density ----------

# Histogram of gene start positions, one panel per chromosome.
# Experiment with bins, colors, or switch to geom_density().
ggplot(genes_filtered, aes(x = start)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
  facet_wrap(~ chrom, scales = "free_x") +
  ggtitle("Gene Start Position Distribution by Chromosome") +
  xlab("Genomic position") + ylab("Gene count") +
  theme_minimal()

# try to remove "scales = "free_x"" and see what happens to the x-axis across panels
# try facet_grid(~ chrom) instead of facet_wrap() and see how the layout changes
# try facet_grid(chrom ~ .)


# ---------- SLIDE: Step 5: Export ----------

p_final <- ggplot(genes_filtered, aes(x = start)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
  facet_wrap(~ chrom, scales = "free_x") +
  theme_minimal()

ggsave("gene_density.pdf", p_final, width = 12, height = 8)
