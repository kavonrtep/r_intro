################################################################################
# R SCRIPT: Analysis and Plotting of ChIP-seq Bedgraph Data for chr1
#
# This script demonstrates:
#   1. Loading a bedgraph file (without a header) into a data frame.
#   2. Defining column names upon import.
#   3. Subsetting the data for chromosome "chr1".
#   4. Creating a color vector using the rep() function.
#   5. Setting the color to red for scores higher than 2.
#   6. Plotting the data with base R graphics.
#
# INSTRUCTIONS:
#   1. Save your bedgraph file (e.g., "chipseq.bedgraph") in your working directory.
#   2. The file is tab-delimited and does NOT include a header.
#   3. Run this script in RStudio, and modify as needed.
################################################################################

# Make sure your working directory is set to the location of this script
# Use getwd() to check your current working directory and setwd() to change it if needed.
# or use the RStudio interface to set the working directory to the folder containing this script.


# ----------------------------- #
# 1. Load the Bedgraph Data     #
# ----------------------------- #
# Read the bedgraph file into a data frame.
# Since the file has no header, set header = FALSE.
chipseq_data <- read.table("data/chip_seq.bedgraph", header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE)

# Define the column names manually.
colnames(chipseq_data) <- c("chr_name", "start", "end", "score")

# Inspect the first few rows and the structure of the data.
head(chipseq_data)   # Check the first 6 rows.
tail(chipseq_data)   # Check the last 6 rows.
str(chipseq_data)    # See the structure and data types.

# ----------------------------- #
# 2. Subset Data for Chromosome 1 #
# ----------------------------- #
# Create a new data frame containing only rows where chr_name equals "chr1".
chr1_data <- subset(chipseq_data, chr_name == "chr1")
# the above line is equivalent to:
# chr1_data <- chipseq_data[chipseq_data$chr_name == "chr1", ]

# Verify the subset.
head(chr1_data)
tail(chr1_data)
str(chr1_data)

# ----------------------------- #
# 3. Create a Color Vector       #
# ----------------------------- #
# Create a vector of colors with the same length as the number of rows in chr1_data.
# By default, all points will be "black".
colors <- rep("black", nrow(chr1_data))

# Now, for rows where the 'score' is greater than 5, set the corresponding color to "red".
colors[chr1_data$score > 5] <- "red"

# ----------------------------- #
# 4. Plotting the Data         #
# ----------------------------- #
# Create a scatter plot of the genomic start positions (x-axis) versus the ChIP-seq scores (y-axis).
# Points will be colored based on the color vector defined above.
plot(chr1_data$start, chr1_data$score,
     col = colors,        # Color points based on our color vector.
     pch = 19,            # Use solid circle symbols.
     xlab = "Genomic Start Position",
     ylab = "ChIP-seq Score",
     main = "ChIP-seq Scores on chr1\n(Red if score > 5)",
     type = "h",
)          # Use a histogram-like plot for better visualization.


# Add a legend to clarify the color coding.
legend("topright", legend = c("Score <= 5", "Score > 5"),
       col = c("black", "red"), pch = 19)

################################################################################

# NOTES:
#   - This script uses functions introduced in the introductory lecture such as
#     read.table(), subset(), rep(), and base plotting functions.
#   - Experiment by modifying the color threshold or exploring other plotting options.
#
#

# TASK:
# 1. Export plot to a PDF file named "chr1_chipseq_plot.pdf" in output directory.
# 2. Change the color threshold to 20 and observe how the plot changes.

