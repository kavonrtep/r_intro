################################################################################
# R SCRIPT FOR SESSION 2: IMPORTING & MANIPULATING EXPERIMENTAL TABULAR DATA
#
# This script covers:
#
# 1. Introduction to the read.table Function Family:
#    - The read.table function is fully configurable for importing tabular data.
#    - Specialized functions such as read.csv and read.delim are built on read.table
#      with preset defaults:
#         * read.csv uses sep = "," and header = TRUE.
#         * read.delim uses sep = "\t" and header = TRUE.
#
# 2. For each dataset in the "data" directory, we will:
#      a) Import the data.
#      b) Inspect the data (using head, tail, str, summary, and checking for missing values).
#      c) Create a simple plot for exploratory analysis.
#
# The datasets include:
#   - DNase.tsv: ELISA assay data for DNase in rat serum.
#   - ChickWeight.csv: Data on the effect of diet on chick growth.
#   - ChickWeight_feed.txt and ChickWeight_feed.xlsx: Feed supplement experiment data.
#   - genes.bed: Genomic coordinates data for genes (BED format).
#
# INSTRUCTIONS:
#   1. Ensure the following files are placed in the "data" folder:
#         DNase.tsv, ChickWeight.csv, ChickWeight_feed.txt,
#         ChickWeight_feed.xlsx, and genes.bed.
#   2. Run each section step by step to explore the data and its corresponding plot.
################################################################################

########################## 1. DNase.tsv #######################################

# --- Import DNase.tsv ---
# DNase.tsv is a tab-delimited file with 176 rows and 3 columns:
#  - Run: Factor indicating the assay run number.
#  - conc: Numeric vector of concentrations.
#  - density: Numeric vector of optical density readings.

# Check how file data/DNase.tsv look like in text editor

dnase <- read.delim("data/DNase.tsv", header = TRUE)
# --- Inspect DNase Data ---
head(dnase)
tail(dnase)
str(dnase)
dim(dnase)
summary(dnase)

# What happens if we use read.csv instead of read.delim?
dnase2 <- read.csv("data/DNase.tsv", header = TRUE)
head(dnase2)

# What happen with header = FALSE?, try it and see the result
dnase3 <- read.delim("data/DNase.tsv", header = FALSE)
head(dnase3)
str(dnase3)

# --- Plot DNase Data (Base R Scatter Plot) ---
# plot function is a base R function for creating simple plots like scatter plots or line plots
# basic syntax: plot(x, y, main = "Title", xlab = "X-axis label", ylab = "Y-axis label", col = "color")
# Plot optical density (density) vs. concentration (conc)

plot(dnase$conc, dnase$density,
     main = "DNase ELISA Assay",
     xlab = "Concentration",
     ylab = "Optical Density",
     col = "blue")


########################## 2. ChickWeight.csv ##################################

# --- Import ChickWeight.csv ---
# ChickWeight.csv contains 578 rows and 4 columns:
#   - weight: Body weights of chicks.
#   - Time: Time points of measurement.
#   - Chick: Chick IDs.
#   - Diet: Diet types.

# Inspect the file data/ChickWeight.csv in a text editor
chick_weight <- read.csv("data/ChickWeight.csv", header = TRUE)
# Convert Chick and Diet to factors
chick_weight$Chick <- as.factor(chick_weight$Chick)
chick_weight$Diet <- as.factor(chick_weight$Diet)

# --- Inspect ChickWeight Data ---
head(chick_weight)
tail(chick_weight)
str(chick_weight)
summary(chick_weight)



# --- Plot ChickWeight Data ---
# Using base R: Scatter plot of Weight vs. Time with colors based on Diet
plot(chick_weight$weight, chick_weight$Time,
     main = "Chick Weight vs. Time (Base R)",
     xlab = "Weight",
     ylab = "Time",
     pch = 16, col = chick_weight$Diet)

legend("topright", legend = levels(chick_weight$Diet),
       pch = 16, col = levels(chick_weight$Diet))

# --- Additional Plot using ggplot2 for ChickWeight Data ---
# This plot creates a faceted scatter plot with a regression line.
# Ensure the ggplot2 package is installed.
# install.packages("ggplot2")  # Uncomment if needed.
library(ggplot2)
ggplot(data = chick_weight, aes(x = weight, y = Time)) +
  geom_point(aes(color = Diet), size = 2) +
  facet_wrap(~ Diet) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  ggtitle("Chick Weight vs. Time by Diet (ggplot2)") +
  xlab("Weight") +
  ylab("Time")

########################## 3. ChickWeight_feed.txt #############################

# --- Import ChickWeight_feed.txt ---
# This is a space-delimited text file containing feed supplement data.
# Columns: weight (numeric), feed (factor)
feed_txt <- read.table("data/ChickWeight_feed.txt", header = TRUE,
                       stringsAsFactors = TRUE)


# --- Inspect Feed Supplement Data (TXT) ---
head(feed_txt)
tail(feed_txt)
str(feed_txt)
summary(feed_txt)


# --- Plot Feed Supplement Data (TXT) ---
# Boxplot of weight by feed type
boxplot(weight ~ feed, data = feed_txt,
        main = "Chick Weight by Feed Type (TXT Data)",
        xlab = "Feed Type",
        ylab = "Weight"
)

########################## 4. ChickWeight_feed.xlsx ############################

# --- Import ChickWeight_feed.xlsx ---
# An Excel file with feed supplement data. Make sure the readxl package is installed.
# install.packages("readxl")  # Uncomment if needed.
library(readxl)
feed_xlsx <- read_excel("data/ChickWeight_feed.xlsx", sheet = 1)
# feed_xlsx is not a data.frame but a tibble
# tibble is a enhanced version of data.frame, it is more user friendly, provides better printing and is more efficient


# --- Inspect Feed Supplement Data (Excel) ---
print("ChickWeight_feed.xlsx Data:")
head(feed_xlsx)
tail(feed_xlsx)
str(feed_xlsx)
summary(feed_xlsx)
print(feed_xlsx)
print(feed_xlsx, n = 20)


# --- Plot Feed Supplement Data (Excel) ---
# Boxplot of weight by feed type for the Excel version
boxplot(weight ~ feed, data = feed_xlsx,
        main = "Chick Weight by Feed Type (Excel Data)",
        xlab = "Feed Type",
        ylab = "Weight")

# make a boxplot with ggplot2
ggplot(data = feed_xlsx, aes(x = feed, y = weight)) +
  geom_boxplot(fill = "lightblue")



########################## 5. genes.bed #########################################

# --- Import genes.bed ---
# BED files are tab-delimited and usually have no header.
genes <- read.table("data/genes.bed", header = FALSE)
head(genes)
# column names are V1, V2, V3, V4, V5, V6 when we use header = FALSE

# Assign column names based on standard BED format (adjust if necessary)
colnames(genes) <- c("chrom", "start", "end", "name", "score", "strand")
head(genes)

# alternativelly header can be defines while reading the file
genes <- read.table("data/genes.bed", header = FALSE,
                    col.names = c("chrom", "start", "end", "name", "score", "strand"))
# --- Inspect genes.bed Data ---
head(genes)

# how many genes are in the file?, how many genes are in each chromosome?

n_genes <- nrow(genes)
print(n_genes)

n_genes_chr <- table(genes$chrom)
print(n_genes_chr)
plot(n_genes_chr)

# Calculate gene length (end - start) and add as a new column
genes$length <- genes$end - genes$start

# --- Plot: Histogram of Gene Lengths
hist(genes$length,
       breaks = 100,
       main = "Histogram of Gene Lengths on Chromosome 1",
       xlab = "Gene Length (bp)",
       col = "orange", xlim = c(0, 20000))

chr12_genes <- genes[genes$chrom == "chr12", ]
hist(chr12_genes$start, breaks = 200, col = "blue", main = "Histogram of Gene Start Positions on Chromosome 1")

########################## 6. DATA FRAMES & TIBBLES ###############################

# Convert ChickWeight data frame to a tibble for improved printing.
# install.packages("tibble")  # Uncomment if needed.
library(tibble)
chick_weight_tb <- as_tibble(chick_weight)
print(chick_weight_tb)

# TASK: Compare the printed output of the data frame and the tibble.

########################## 7. DATA EXPORT ######################################

# For demonstration, export the feed_txt data as a CSV file and a TSV file, with and without row names.
# Export file to work_dir folder
write.csv(feed_txt, file = "work_dir/ChickWeight_feed.csv", row.names = FALSE)
write.csv(feed_txt, file = "work_dir/ChickWeight_feed_with_row_names.csv", row.names = TRUE)
write.table(feed_txt, file = "work_dir/ChickWeight_feed.tsv", sep = "\t", row.names = FALSE)




########################## 8. HANDLING MISSING DATA ##############################
# inspect the data file DNase_missing_values.tsv in a text editor
dnase <- read.delim("data/DNase_missing_values_2.tsv", header = TRUE)

str(dnase)
summary(dnase)
# density column contains missing values - NA
mean(dnase$density)  # This will not produce the correct mean due to missing values
mean(dnase$density, na.rm = TRUE)  # Calculate mean without missing values
# we can use na.rm = TRUE in many functions like mean, median, sum, etc

# or we can remove rows with missing values completely
dnase_no_na <- na.omit(dnase)  # Remove rows with missing values
str(dnase_no_na)
str(dnase)

# inspect the data file DNase_missing_values.tsv in a text editor
dnase2 <- read.delim("data/DNase_missing_values.tsv", header = TRUE)
str(dnase2)

# if NA values are represented by a different character, we can specify it using
# the na.strings argument!

########################## 9. Exporting R Objects ################################

## Saving and loading R objects with saveRDS and readRDS
## this is a way to save and load R objects (data frames, lists, etc) in a binary format
## objects saved with saveRDS can be loaded with readRDS

# Save the dnase_no_na object as a .rds file
suns_spots <- readRDS("data/sun_splots.rds")
summary(suns_spots)
suns_spots

# export dnase_no_na object as rds file
saveRDS(dnase_no_na, file = "work_dir/dnase_no_na.rds")
# load is back to R
dnase_no_na_loaded <- readRDS("work_dir/dnase_no_na.rds")





