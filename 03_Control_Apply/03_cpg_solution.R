library(stringi)

dna_sequence <- readRDS("data/dna_sequence.rds")

window_size <- 100
n_windows <- nchar(dna_sequence) - window_size + 1
cpg_counts <- numeric(n_windows)

for (i in seq_len(n_windows)) {
  window <- substr(dna_sequence, i, i + window_size - 1)
  cpg_counts[i] <- stri_count_fixed(window, "CG")
}

plot(cpg_counts, type = "l",
     xlab = "Position (bp)", ylab = "CpG count (per 100 bp)",
     main = "CpG Island Detection by Sliding Window",
     col = "darkblue")
abline(h = mean(cpg_counts), col = "red", lty = 2)
legend("topright", legend = "genome average", col = "red", lty = 2)