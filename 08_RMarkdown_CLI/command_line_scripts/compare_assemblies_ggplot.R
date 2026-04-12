#!/usr/bin/env Rscript
# compare_assemblies_ggplot.R
#
# Calculates genome assembly statistics (N50 and L50) for two assemblies
# and plots cumulative contig length curves using ggplot2.
#
# Usage:
#   Rscript compare_assemblies_ggplot.R --g1 <genome1.fai> --g2 <genome2.fai> [--outbase <basename>]
#
# --g1: Path to genome 1 FASTA index file (.fai)
# --g2: Path to genome 2 FASTA index file (.fai)
# --outbase: Base name for output files (default: "assembly_output")

library(optparse)
library(ggplot2)

option_list <- list(
  make_option(c("--g1"), type = "character", default = NULL,
              help = "Path to genome 1 FASTA index (.fai) file", metavar = "character"),
  make_option(c("--g2"), type = "character", default = NULL,
              help = "Path to genome 2 FASTA index (.fai) file", metavar = "character"),
  make_option(c("--outbase"), type = "character", default = "assembly_output",
              help = "Base name for output files [default = %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$g1) || is.null(opt$g2)) {
  print_help(opt_parser)
  stop("Both --g1 and --g2 arguments must be provided.", call. = FALSE)
}

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

stats_g1 <- calculate_assembly_stats(opt$g1)
stats_g2 <- calculate_assembly_stats(opt$g2)

# Save summary CSV
summary_df <- data.frame(
  Genome       = c("Genome 1", "Genome 2"),
  Total_Length = c(stats_g1$total_length, stats_g2$total_length),
  N50          = c(stats_g1$N50, stats_g2$N50),
  L50          = c(stats_g1$L50, stats_g2$L50)
)
summary_filename <- paste0(opt$outbase, "_summary.csv")
write.csv(summary_df, file = summary_filename, row.names = FALSE)
cat("Assembly statistics saved to", summary_filename, "\n")

# Build a tidy data frame for ggplot2
plot_df <- rbind(
  data.frame(Genome = "Genome 1",
             Index  = seq_along(stats_g1$cumsum),
             Cumsum = stats_g1$cumsum),
  data.frame(Genome = "Genome 2",
             Index  = seq_along(stats_g2$cumsum),
             Cumsum = stats_g2$cumsum)
)

p <- ggplot(plot_df, aes(x = Index, y = Cumsum, color = Genome)) +
  geom_line() +
  labs(
    title = "Cumulative Contig Length",
    x     = "Contig Index",
    y     = "Cumulative Length (bp)",
    color = NULL
  ) +
  theme_bw()

plot_filename <- paste0(opt$outbase, "_plot.png")
ggsave(plot_filename, plot = p, width = 8, height = 5)
cat("Plot saved to", plot_filename, "\n")
