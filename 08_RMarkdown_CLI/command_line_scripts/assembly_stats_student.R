#!/usr/bin/env Rscript
# assembly_stats_student.R
#
# TASK: Complete this script so that it:
#   1. Reads any number of FASTA index (.fai) files from positional arguments
#   2. Computes N50, L50, and total length for each assembly
#   3. Saves a summary CSV table
#   4. Plots cumulative contig length curves using ggplot2
#
# Usage (once complete):
#   Rscript assembly_stats_student.R --outbase results/comparison asm1.fai asm2.fai asm3.fai
#   Rscript assembly_stats_student.R --outbase results/comparison data/assemblies/*.fai

library(optparse)
library(ggplot2)

# --- Command-line options ---
# Named options (--flag) are defined in option_list.
# Positional arguments (the .fai files) are passed after all flags and
# collected automatically when positional_arguments = TRUE.

option_list <- list(
  make_option(c("--outbase"), type = "character", default = "assembly_output",
              help = "Base name for output files [default = %default]",
              metavar = "character")
)

opt_parser <- OptionParser(
  usage       = "Usage: %prog [options] file1.fai file2.fai ...",
  option_list = option_list
)

# parse_args() returns a list with two elements when positional_arguments = TRUE:
#   $options — named flags (e.g. opt$outbase)
#   $args    — character vector of positional file paths
args       <- parse_args(opt_parser, positional_arguments = TRUE)
opt        <- args$options
fai_files  <- args$args

# Stop early with a helpful message if no .fai files were provided
if (length(fai_files) < 1) {
  print_help(opt_parser)
  stop("At least one .fai file must be provided as a positional argument.", call. = FALSE)
}


# --- PROVIDED: Assembly statistics function (do not modify) ---
# Given a .fai file path, returns a list with:
#   $total_length   — total assembly size in bp
#   $N50            — N50 contig length
#   $L50            — number of contigs needed to reach 50% of total length
#   $sorted_lengths — contig lengths sorted descending
#   $cumsum         — cumulative sum of sorted_lengths
calculate_assembly_stats <- function(fai_file) {
  fai <- read.table(fai_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(fai) < 2) stop(paste("File", fai_file, "does not have at least 2 columns"))

  contig_lengths <- as.numeric(fai[[2]])
  total_length   <- sum(contig_lengths, na.rm = TRUE)
  sorted_lengths <- sort(contig_lengths, decreasing = TRUE)
  cumsum_lengths <- cumsum(sorted_lengths)
  L50 <- which(cumsum_lengths >= total_length / 2)[1]
  N50 <- sorted_lengths[L50]

  list(total_length = total_length, N50 = N50, L50 = L50,
       sorted_lengths = sorted_lengths, cumsum = cumsum_lengths)
}

# --- TASK 2: Calculate stats for each assembly ---
# Loop over fai_files and call calculate_assembly_stats() for each.
# Store results in a list, e.g.:
#   stats_list <- lapply(fai_files, calculate_assembly_stats)
# Use basename(fai_file) as a human-readable assembly name.

# stats_list <- ...


# --- TASK 3: Build a summary data frame and save as CSV ---
# Create a data frame with one row per assembly containing:
#   Assembly, Total_Length, N50, L50
# Hint: use sapply() to extract values from stats_list
# Save with write.csv(..., row.names = FALSE)
# Filename: paste0(opt$outbase, "_summary.csv")


# --- TASK 4: Build a tidy data frame for ggplot2 ---
# For each assembly, create a data frame with columns:
#   Assembly — assembly name (use basename of the file)
#   Index    — contig index (1, 2, 3, ...)
#   Cumsum   — cumulative length at that index
# Combine all assemblies with rbind() or do.call(rbind, ...)
# Hint: stats_list[[i]]$cumsum gives the cumulative lengths for assembly i


# --- TASK 5: Plot with ggplot2 and save ---
# Plot cumulative contig length (y) vs contig index (x),
# with one line per assembly colored by the Assembly column.
# Use geom_line(), labs(), and theme_bw().
# Save with ggsave(), filename: paste0(opt$outbase, "_plot.png")
