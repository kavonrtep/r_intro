# =====================================================
# Introduction to ggplot2 Graphics in R - Practice Script
# =====================================================
# This script introduces key ggplot2 concepts:
#   • Data: The dataset used for plotting.
#   • Mapping: How variables are assigned to aesthetics.
#   • Layers: Adding geometric objects and statistical transformations.
#   • Scales: Converting data values to visual properties.
#   • Facets: Splitting data into panels.
#   • Coordinates and Themes: Adjusting plot layout and non-data elements.
#
# Each section includes a complete example and a TASK for independent exploration.
# Run each section and then attempt the tasks to deepen your understanding.
#
# Load required libraries:
library(ggplot2)
library(gridExtra)  # For combining multiple plots

# -----------------------------------------------------
# Section 1: Data, Mapping, and Layers
# -----------------------------------------------------
# In ggplot2, every plot begins with the ggplot() function. The first argument is the dataset,
# and the mapping argument (created by aes()) defines how variables map to visual properties.
#
# Example: Create a basic scatter plot using the built-in mpg dataset.
?mpg
print(mpg)
ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point() +                               # Add points layer
  ggtitle("Scatter Plot: City vs Highway MPG") +
  xlab("City MPG") +
  ylab("Highway MPG")


# TASK 1:
# Modify the above scatter plot to:
#   - Change the point color to "darkgreen".
#   - Increase the point size
# hint : modify geom_point()


# -----------------------------------------------------
# Section 2: Mapping Aesthetics
# -----------------------------------------------------
# Mapping aesthetics assign variables to visual properties like color, size, or shape.
# Example: Map the 'class' variable to color.
ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  ggtitle("Scatter Plot: Color by Car Class") +
  xlab("City MPG") +
  ylab("Highway MPG")


# TASK 2:
# Modify the above plot to:
#   - Also map 'displ' (engine displacement) to the size of the points.
#   - Change the color scale using a different palette (e.g., scale_color_brewer()).

# Some scale_color_* functions are suitable for discrete variables,
# while others are for continuous variables!



# -----------------------------------------------------
# Section 3: Layers: Geoms and Statistical Transformations
# -----------------------------------------------------
# Layers transform the mapped data into visual elements.
# Geometric objects (geoms) like geom_point(), geom_line(), and geom_bar() display the data.
#
# Example: Add a smooth trendline layer to the scatter plot.
ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE) +  # Overlay linear model trendline with confidence interval
  ggtitle("Scatter Plot with Trend Line") +
  xlab("City MPG") +
  ylab("Highway MPG")


# TASK 3:
# Experiment by:
#   - Changing the trendline method to "loess" in geom_smooth().
# Loess stands for locally estimated scatterplot smoothing





# -----------------------------------------------------
# Section 4: Scales
# -----------------------------------------------------
# Scales control how data values are mapped to visual aesthetics.
# They can modify axis transformations, breaks, and legends.
#
# Example: Change the color scale to use the viridis palette for discrete variables.
ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  scale_color_viridis_d() +
  ggtitle("Scatter Plot with Viridis Color Scale") +
  xlab("City MPG") +
  ylab("Highway MPG")


# Example: Apply a log transformation to the x-axis.
ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point() +
  scale_x_log10() +
  ggtitle("Scatter Plot with Log-Scaled X-Axis") +
  xlab("City MPG (log scale)") +
  ylab("Highway MPG")


# TASK 4:
# Modify one of the above plots to:
#   - Set custom axis limits with scale_x_continuous() or scale_y_continuous().




# -----------------------------------------------------
# Section 5: Facets
# -----------------------------------------------------
# Faceting splits data into multiple panels based on one or more variables.
#
# Example: Create a faceted scatter plot of the iris dataset by species.
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3) +
  facet_wrap(~Species) +
  ggtitle("Faceted Scatter Plot: Iris Data") +
  xlab("Sepal Length") +
  ylab("Sepal Width")

# TASK 5:  What happen if you use facet_wrap with ncol=1 or nrow=1?


# using facet_grid() instead of facet_wrap() will create a grid of plots based on two variables.
ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  facet_grid(drv ~ cyl) +
  ggtitle("Faceted Scatter Plot: MPG by Drive and Cylinder Count") +
  xlab("City MPG") +
  ylab("Highway MPG")





# -----------------------------------------------------
# Section 6: Coordinates and Themes
# -----------------------------------------------------
# Coordinates change the plot’s view (e.g., flipping axes) while themes adjust non-data elements.
#
# Example: Create a box plot of mpg by cylinder count (using mtcars), flip the coordinates, and apply a minimal theme.
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Box Plot: MPG by Cylinder Count (Flipped)")


# Example: Customize a plot’s appearance with a custom theme.
ggplot(mpg, aes(x = cty, y = hwy, color = class)) +
  geom_point() +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.line = element_line(size = 0.75, color = "gray")
  ) +
  ggtitle("Scatter Plot with Custom Theme")


# TASK 6:
# Modify one of the above plots to:
#   - Change the background theme (try theme_dark()
#   - Adjust text sizes and fonts for titles and axis labels using theme() and
#  setting text = element_text(size = 12, family = "Arial")




# -----------------------------------------------------
# Section 7: Other Plot Types
# -----------------------------------------------------
# ggplot2 can create various plot types.
#
# Bar Plot Example: Plot the count of car classes.
ggplot(mpg, aes(x = class)) +
  geom_bar(fill = "steelblue") +
  ggtitle("Bar Plot: Count of Car Classes") +
  xlab("Car Class") +
  ylab("Count")


# Line Plot Example: Plot y = x^2.
df_line <- data.frame(x = 1:10, y = (1:10)^2)
ggplot(df_line, aes(x = x, y = y)) +
  geom_line(color = "red", size = 1.2) +
  ggtitle("Line Plot: y = x^2") +
  xlab("X") +
  ylab("Y")

# TASK - add points to the line plot above using geom_point()



# Box Plot Example: MPG by cylinder count.
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(fill = "orange") +
  ggtitle("Box Plot: MPG by Cylinder Count") +
  xlab("Cylinders") +
  ylab("MPG")


# Histogram Example: Histogram of iris Sepal.Length.
ggplot(iris, aes(x = Sepal.Length)) +
  geom_histogram(bins = 20, fill = "purple", alpha = 0.7) +
  ggtitle("Histogram of Sepal Length") +
  xlab("Sepal Length") +
  ylab("Frequency")

# TASK - create multiple histagrams using `Species`
# Try different fill colors or facet_wrap() by `Species`
# test geom_histogram position settings (dodge, fill, stack, identity)




# Density Plot Example: Density of iris Sepal.Length by species.
ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_density(alpha = 0.5) +
  ggtitle("Density Plot of Sepal Length by Species") +
  xlab("Sepal Length") +
  ylab("Density")


# TASK 7:
# Choose one plot type from above (bar, line, box, histogram, or density).
# Modify it by adding additional customization (colors, labels, transparency, etc.)
# or create a new plot using a dataset of your choice.

# -----------------------------------------------------
# Section 8: Combining Plots
# -----------------------------------------------------
# Use gridExtra to arrange multiple ggplot objects in a single layout.
#
# Example: Combine the density, line, and box plots.
p_density <- ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_density(alpha = 0.5) +
  ggtitle("Density Plot of Sepal Length by Species") +
  xlab("Sepal Length") +
  ylab("Density")
p_line <- ggplot(df_line, aes(x = x, y = y)) +
  geom_line(color = "red", size = 1.2) +
  ggtitle("Line Plot: y = x^2") +
  xlab("X") +
  ylab("Y")
p_box2 <- ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
    geom_boxplot(fill = "orange") +
    ggtitle("Box Plot: MPG by Cylinder Count") +
    xlab("Cylinders") +
    ylab("MPG")


combined_plot <- grid.arrange(p_density, p_line, p_box2, nrow = 2)



# TASK 8:
# Experiment with different layouts by:
#   - Changing the number of rows or columns (using nrow or ncol).
#   - Using layout_matrix to create custom arrangements.
#   - Adding additional plots from your previous tasks.

# -----------------------------------------------------
# Section 9: Exporting Plots
# -----------------------------------------------------
# The ggsave() function saves a ggplot object to a file.
#
# Example: Save the combined plot as both PNG and PDF files.
ggsave("combined_plot.png", combined_plot, width = 10, height = 8, device = "png")
ggsave("combined_plot.pdf", combined_plot, width = 10, height = 8, device = "pdf")

# TASK 9:
# Export one of your individual plots (e.g., the faceted plot or custom scatter plot) to a PDF file.
# Try different file formats (e.g., jpg, tiff) and adjust the dimensions as needed.


################################################
## Final Task
################################################
# =====================================================
# ggplot2 Visualization of Gene Densities along Chromosomes
# =====================================================
# This script guides you through:
#   1. Importing and preparing gene annotation data from a BED file.
#   2. Filtering the data to include only genes from chromosomes chr1 to chr6.
#   3. Creating a bar plot of gene counts per chromosome.
#   4. Creating a histogram (or density plot) of gene start positions for a single chromosome.
#   5. Creating a faceted plot (multiplot) for gene start positions across chromosomes chr1 - chr6.

#   - Experiment with different plot types (histogram vs. density).

------------------------------------------------
# Step 1: Import and Prepare the Data
# -----------------------------------------------------
# Import the BED file. BED files are tab-delimited and typically have no header.
genes <- read.table("data/genes.bed", header = FALSE)

# Assign meaningful column names based on the standard BED format.
colnames(genes) <- c("chrom", "start", "end", "name", "score", "strand")
# Inspect the first few rows of the data.
head(genes)

# -----------------------------------------------------
# Step 2: Filter Data for Chromosomes chr1 to chr6
# -----------------------------------------------------
# Filter out records that are not related to chromosomes chr1, chr2, chr3, chr4, chr5, or chr6.
genes_filtered <- subset(genes, chrom %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6"))
# Check the first few rows of the filtered data.
head(genes_filtered)
# Print the number of genes in the filtered dataset.
cat("Number of genes in chr1-chr6:", nrow(genes_filtered), "\n")

# -----------------------------------------------------
# Step 3: Bar Plot of Gene Counts per Chromosome
# -----------------------------------------------------
# Create a summary table of gene counts per chromosome.
gene_counts <- as.data.frame(table(genes_filtered$chrom))
colnames(gene_counts) <- c("chrom", "count")
# Create a bar plot using ggplot2.

ggplot(genes_filtered, aes(x = chrom, fill = chrom)) +
  geom_bar() +
  ggtitle("Number of Genes per Chromosome (chr1 - chr6)") +
  xlab("Chromosome") +
  ylab("Number of Genes") +
  theme_minimal()


# -----------------------------------------------------
# Step 4: Create a histogram of gene densities for each chromosome 1-6
# by plotting the distribution of gene start positions.
# use either geom_histogram() or geom_density() to visualize the distribution.
# used facet_wrap() to create separate panels for each chromosome.
# Experiment by changing the number of bins, colors, and labels. or bw parameter for density plot
# Try using geom_density() instead of geom_histogram() to visualize the distribution.
# - Modify the faceted plot by experimenting with different facet_wrap parameters (e.g., ncol, nrow).
# - Try using geom_density() for each panel instead of geom_histogram().

# Step - export visualization to a file
