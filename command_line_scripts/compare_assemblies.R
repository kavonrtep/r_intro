#!/usr/bin/env Rscript
# assembly_stats_optparse.R
#
# This script calculates genome assembly statistics (N50 and L50) for two genomes.
# It uses the optparse library to parse command-line options.
#
# Usage:
#   Rscript assembly_stats_optparse.R --g1 <genome1.fai> --g2 <genome2.fai> [--outbase <basename>]
#
# --g1: Path to genome 1 FASTA index file (.fai)
# --g2: Path to genome 2 FASTA index file (.fai)
# --outbase: Base name for output files (default: "assembly_output")
#
# For each genome, the script reads the FASTA index file (assumed to be tab-delimited with contig length in the second column),
# calculates:
#   - Total assembly length,
#   - N50: The contig length at which the cumulative length reaches at least 50% of total assembly length,
#   - L50: The minimal number of contigs needed to reach that 50% threshold.
#
# It then creates a summary table (saved as CSV) and a single plot showing both genomes’ cumulative contig length curves.
# The plot is saved as a PNG file.
#
# Base graphics are used for plotting.

# Load required libraries

library(optparse)


# Define command-line options using optparse
option_list <- list(
  make_option(c("--g1"), type = "character", default = NULL,
              help = "Path to genome 1 FASTA index (.fai) file", metavar = "character"),
  make_option(c("--g2"), type = "character", default = NULL,
              help = "Path to genome 2 FASTA index (.fai) file", metavar = "character"),
  make_option(c("--outbase"), type = "character", default = "assembly_output",
              help = "Base name for output files [default = %default]", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that both genome files
if (is.null(opt$g1) || is.null(opt$g2)) {
  print_help(opt_parser)
  stop("Both --g1 and --g2 arguments must be provided.", call. = FALSE)
}

# Save the output base name for later use
outbase <- opt$outbase

# Function to calculate assembly statistics from a FASTA index file
calculate_assembly_stats <- function(fai_file) {
  # Read the FASTA index file (tab-delimited, no header expected)
  fai <- read.table(fai_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  # Check that the file has at least two columns
  if (ncol(fai) < 2) {
    stop(paste("File", fai_file, "does not have at least 2 columns"))
  }

  # Extract contig lengths (assumed to be in the second column)
  contig_lengths <- as.numeric(fai[[2]])

  # Calculate total assembly length
  total_length <- sum(contig_lengths, na.rm = TRUE)

  # Sort contig lengths in descending order
  sorted_lengths <- sort(contig_lengths, decreasing = TRUE)

  # Compute the cumulative sum of sorted lengths
  cumsum_lengths <- cumsum(sorted_lengths)

  # Identify the index where the cumulative sum reaches at least 50% of the total length
  L50 <- which(cumsum_lengths >= total_length / 2)[1]

  # N50 is the length of the contig at this index
  N50 <- sorted_lengths[L50]

  # Return statistics along with sorted values (for plotting)
  list(total_length = total_length, N50 = N50, L50 = L50,
       sorted_lengths = sorted_lengths, cumsum = cumsum_lengths)
}

# Calculate statistics for Genome 1 and Genome 2
stats_g1 <- calculate_assembly_stats(opt$g1)
stats_g2 <- calculate_assembly_stats(opt$g2)

# Create a summary data frame for both genomes
summary_df <- data.frame(
  Genome = c("Genome 1", "Genome 2"),
  Total_Length = c(stats_g1$total_length, stats_g2$total_length),
  N50 = c(stats_g1$N50, stats_g2$N50),
  L50 = c(stats_g1$L50, stats_g2$L50)
)

# Save the summary table as a CSV file
summary_filename <- paste0(outbase, "_summary.csv")
write.csv(summary_df, file = summary_filename, row.names = FALSE)
cat("Assembly statistics saved to", summary_filename, "\n")

# Create a single plot for both genomes: cumulative contig length vs contig index

# Determine maximum number of contigs to plot (for setting the x-axis limits)
max_contigs <- max(length(stats_g1$cumsum), length(stats_g2$cumsum))
max_length <- max(stats_g1$total_length, stats_g2$total_length)

# Set output file for the plot (PNG)
plot_filename <- paste0(outbase, "_plot.png")
png(plot_filename, width = 800, height = 600)

# Plot Genome 1 cumulative sum
plot(stats_g1$cumsum, type = "b", pch = 19, col = "blue",
     xlab = "Contig Index", ylab = "Cumulative Length",
     xlim = c(1, max_contigs), ylim = c(0, max_length),
     main = "Cumulative Contig Length for Both Genomes")
# Add Genome 1 horizontal line at 50% threshold

# Add Genome 2 cumulative sum to the same plot using lines()
lines(stats_g2$cumsum, type = "b", pch = 19, col = "darkgreen")
# Add Genome 2 horizontal line at its 50% threshold

# Add a legend
legend("bottomright", legend = c("Genome 1", "Genome 2"),
       col = c("blue", "darkgreen"), pch = 19, lty = 1, bty = "n")

dev.off()
cat("Plot saved to", plot_filename, "\n")
