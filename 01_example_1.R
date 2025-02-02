################################################################################
# R SCRIPT: Loading CSV Data, Simple Analysis, and Base Plotting
#
# This script demonstrates:
#  1. Reading a CSV file into a data frame.
#  2. Exploring the data with basic functions like head(), str(), and summary().
#  3. Performing simple analysis (e.g., calculating group means).
#  4. Creating plots using base graphics with custom colors.
#
# INSTRUCTIONS:
#   1. Save the CSV content above as "sample_data.csv" in your working directory.
#   2. Open this script in RStudio and run it section by section.
#   3. Modify the code or try additional analyses as a mini-exercise.
################################################################################

# ----------------------------- #
# 1. Loading the CSV Data File  #
# ----------------------------- #

# Read the CSV file into a data frame called 'df'
df <- read.csv("data/sample_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Explore the first few rows of the data
head(df)  # Displays the first 6 rows

# Check the structure of the data frame
str(df)

# Get a summary of the data
summary(df)

# ----------------------------- #
# 2. Simple Data Analysis       #
# ----------------------------- #

# Calculate the mean height for each Treatment group using aggregate()
mean_height_by_treatment <- aggregate(Height ~ Treatment, data = df, FUN = mean)
print(mean_height_by_treatment)

# Calculate the mean height for each Species
mean_height_by_species <- aggregate(Height ~ Species, data = df, FUN = mean)
print(mean_height_by_species)

# ----------------------------- #
# 3. Basic Plotting with Colors #
# ----------------------------- #

# Scatter Plot: Height vs. Leaf_Count, colored by Species.
# Define colors: assign "blue" for Species_A and "red" for Species_B.
# Create a vector of colors corresponding to each row in df.
point_colors <- ifelse(df$Species == "Species_A", "blue", "red")

# Create the scatter plot
plot(df$Height, df$Leaf_Count,
     col = point_colors,       # Colors based on species
     pch = 19,                 # Use solid circle symbols
     xlab = "Plant Height (cm)",
     ylab = "Leaf Count",
     main = "Scatter Plot: Plant Height vs. Leaf Count")

# Add a legend to explain the colors
legend("topleft", legend = c("Species_A", "Species_B"),
       col = c("blue", "red"), pch = 19)

# Boxplot: Compare plant Height across Treatment groups.
boxplot(Height ~ Treatment, data = df,
        col = c("lightblue", "lightgreen"),
        xlab = "Treatment Group",
        ylab = "Plant Height (cm)",
        main = "Boxplot of Plant Height by Treatment")

# How to save the plot as an image file:
png("output/plot_height_vs_leafcount.png", width = 800, height = 600)
plot(df$Height, df$Leaf_Count,
     col = point_colors, pch = 19,
     xlab = "Plant Height (cm)", ylab = "Leaf Count",
     main = "Scatter Plot: Plant Height vs. Leaf Count")
legend("topleft", legend = c("Species_A", "Species_B"),
       col = c("blue", "red"), pch = 19)
dev.off()


# if you need a vector graphics format (e.g., PDF), use pdf() and dev.off() instead: