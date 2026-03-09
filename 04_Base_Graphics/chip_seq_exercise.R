################################################################################
# EXERCISE: Visualizing ChIP-seq Bedgraph Data
#
# The bedgraph file contains ChIP-seq enrichment scores across the genome.
# Each row describes a genomic interval (chromosome, start, end) with an
# associated score representing signal enrichment.
#
# File: data/chip_seq.bedgraph (tab-delimited, no header)
# Columns: chr_name, start, end, score
################################################################################

# Step 1: Load the bedgraph data
my_data <- read.table("data/chip_seq.bedgraph", header = FALSE, sep = "\t")

# Step 2: Assign column names
colnames(my_data) <- c("chr_name", "start", "end", "score")

# Step 3: Inspect the data
head(my_data)
str(my_data)

# Step 4: Subset to chr1
chr1 <- my_data[my_data$chr_name == "chr1", ]

# Step 5: Basic plot - start position vs. score
plot(chr1$start, chr1$score,
     type = "h",
     xlab = "Position (bp)",
     ylab = "Enrichment score",
     main = "ChIP-seq Signal on chr1")

# Step 6: Color by score threshold and re-plot
# Regions with score > 5 are highlighted in red as potential peaks
colors <- ifelse(chr1$score > 5, "red", "black")

plot(chr1$start, chr1$score,
     type = "h",
     col = colors,
     xlab = "Position (bp)",
     ylab = "Enrichment score",
     main = "ChIP-seq Signal on chr1\n(red = score > 5)")

legend("topright",
       legend = c("score <= 5", "score > 5"),
       col = c("black", "red"),
       lty = 1)

# Step 7: How many chromosomes are in the dataset?
chr_counts <- table(my_data$chr_name)
print(chr_counts)

# Multipanel plot for chr1 - chr4
par(mfrow = c(2, 2))

for (chr in paste0("chr", 1:4)) {
  chr_data <- my_data[my_data$chr_name == chr, ]
  colors <- ifelse(chr_data$score > 5, "red", "black")

  plot(chr_data$start, chr_data$score,
       type = "h",
       col = colors,
       xlab = "Position (bp)",
       ylab = "Enrichment score",
       main = chr)
}

par(mfrow = c(1, 1))