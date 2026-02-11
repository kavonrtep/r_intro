################################################################################
# R SCRIPT FOR SESSION 2: CONTROL STRUCTURES, VECTORIZED COMPUTATION,
# APPLY FAMILY OF FUNCTIONS, AND WRITING YOUR OWN FUNCTIONS
#
# This script covers:
# 1. Control Structures: Conditional statements (if, else if, else) and loops
#    (for, while, repeat) with best practices.
#
# 2. Vectorized Computation:
#    - Performing element-wise operations on vectors.
#    - Understanding recycling rules.
#    - Comparison with loops for performance.
#
# 3. The apply Family of Functions:
#    - Examples with apply, lapply, sapply, tapply, and mapply.
#    - When and why to use each function.
#
# 4. Writing Your Own Functions:
#    - How to define functions with default arguments, argument matching,
#      and returning values.
#
# DATA SOURCES USED:
# - sample_data.csv: A simple dataset with plant measurements.
# - feature_counts.csv: Expression data (tab-delimited) with gene names as row names,
#   and samples as column names.
# - feature_annotation.csv: Feature annotation for the expression data.
#
# INSTRUCTIONS:
# 1. Place sample_data.csv, feature_counts.csv, and feature_annotation.csv in your
#    working directory (or update the file paths accordingly).
# 2. Run the script section by section to understand each concept.
################################################################################

####################### 1. CONTROL STRUCTURES ##############################

# --- CONDITIONAL STATEMENTS ---
# Example using if, else if, and else.
number <- 15

if (number < 0) {
  # This block executes if number is negative.
  print("The number is negative.")
} else {
  # This block executes if number is positive.
  print("The number 0 or positive.")
}

# Comparison operators: <, >, <=, >=, ==, !=
# examples
10 < 20
10 == 10
# common mistake: using = instead of == in comparisons

10 != 10
10 <= 10

# Logical operators: AND (&), OR (|), NOT (!)
# examples
TRUE & FALSE
TRUE | FALSE
!TRUE

# --- LOOPS ---

# Example 1: For loop to iterate over a vector
# Create a vector of numbers
num_vector <- c(2, 4, 6, 8, 10)
# Initialize an empty vector to store squares
squares <- numeric(length(num_vector))

# For loop to calculate square of each element
for (i in seq_along(num_vector)) {
  squares[i] <- num_vector[i]^2
  # Print each computation step (for educational purposes)
  cat("Square of", num_vector[i], "is", squares[i], "\n")
}

# Example 2: While loop
# Calculate the sum of numbers until the sum exceeds 20
counter <- 1
sum_val <- 0
while (sum_val <= 20) {
  sum_val <- sum_val + counter
  cat("Adding", counter, "gives sum =", sum_val, "\n")
  counter <- counter + 1
}

# Example 3: Repeat loop with an exit condition
# Use a repeat loop to compute factorial of a number
n <- 5
factorial_result <- 1
i <- 1
repeat {
  factorial_result <- factorial_result * i
  i <- i + 1
  # Exit the loop when i exceeds n
  if (i > n) break
}
cat("Factorial of", n, "is", factorial_result, "\n")

# Best practice reminder:
# Avoid overcomplicating loops when vectorized solutions exist, as they are usually more concise and faster.

################### 2. VECTORIZED COMPUTATION ##############################

# --- VECTORIZED OPERATIONS ---
# Create a simple numeric vector
v <- c(1, 2, 3, 4, 5)

# Vectorized computation: Square each element
v_squared <- v^2
cat("Vectorized squares:", v_squared, "\n")

# --- COMPARISON WITH LOOPS ---
# Using a loop to compute squares (less efficient for large vectors)
v_squared_loop <- numeric(length(v))
for (i in seq_along(v)) {
  v_squared_loop[i] <- v[i]^2
}
cat("Loop computed squares:", v_squared_loop, "\n")


#################### 3. APPLY FAMILY OF FUNCTIONS ###########################

# --- APPLY ---
# Example: Compute the mean of each row in a matrix.
# Create a synthetic matrix (5 rows, 4 columns)
set.seed(123)  # for reproducibility
synthetic_matrix <- matrix(rnorm(20), nrow = 5, ncol = 4)
cat("Synthetic Matrix:\n")
print(synthetic_matrix)
# Apply the mean function to each row (MARGIN = 1)
row_means <- apply(synthetic_matrix, 1, mean)
cat("Row means using apply:", row_means, "\n")

# --- LAPPLY & SAPPLY ---
# Create a list of numeric vectors
num_list <- list(a = 1:5, b = seq(2, 10, by = 2), c = rnorm(4))
# lapply: returns a list with the mean of each vector
list_means <- lapply(num_list, mean)
cat("Means computed with lapply:\n")
print(list_means)
# sapply: simplifies the result to a vector (if possible)
vector_means <- sapply(num_list, mean)
cat("Means computed with sapply:", vector_means, "\n")

# --- TAPPLY ---
# Use sample_data.csv to demonstrate tapply
# sample_data.csv is assumed to have columns: ID, Species, Height, Leaf_Count, Treatment
sample_data <- read.csv("data/sample_data.csv", header = TRUE)
# Display first few rows of sample_data
head(sample_data)
# Calculate average Height for each Treatment group using tapply
avg_height_by_treatment <- tapply(sample_data$Height, sample_data$Treatment, mean)
cat("Average Height by Treatment group:\n")
print(avg_height_by_treatment)

# --- MAPPLY ---
# mapply can apply a function to multiple arguments in parallel.
# Example: Create a function to compute the power of a number with a given exponent.
power_func <- function(x, y) {
  return(x^y)
}
# Use mapply to apply the function to vectors of bases and exponents
bases <- c(2, 3, 4, 5)
exponents <- c(3, 2, 2, 1)
result_power <- mapply(power_func, bases, exponents)
cat("Result of mapply (power computation):\n")
print(result_power)

# --- APPLY FUNCTIONS ON EXPRESSION DATA ---
# Read the expression dataset: feature_counts.csv (tab-delimited, gene names as row names)
# Make sure the file is in your working directory.
feature_counts <- read.delim("data/feature_counts.csv", header = TRUE)
# Inspect the first few rows of feature_counts
head(feature_counts)
# Use apply to compute the total expression (row sums) for each gene across samples
gene_totals <- apply(feature_counts, 1, sum)

hist(log10(gene_totals), main = "Histogram of Total Expression per Gene", xlab = "Total Expression")


# Read the corresponding feature annotation file (assumed file name: feature_annotation.csv)
feature_annotation <- read.delim("data/feature_annotation.csv", header = TRUE, sep = "\t")
# Inspect feature annotation
head(feature_annotation)

#################### 4. WRITING YOUR OWN FUNCTIONS ##########################

# --- FUNCTION DEFINITION: Basic Example ---
# Define a function to compute the mean of a numeric vector while removing NA values.
my_mean <- function(x, na.rm = TRUE) {
  # Check if the input is numeric
  if (!is.numeric(x)) {
    stop("Input vector 'x' must be numeric.")
  }
  # Compute the mean with or without removing NA values based on the argument
  result <- mean(x, na.rm = na.rm)
  return(result)
}

# Test the custom function with a vector containing NA
test_vector <- c(2, 4, NA, 8, 10)
mean_value <- my_mean(test_vector)
cat("Custom function my_mean output:", mean_value, "\n")

# --- FUNCTION WITH DEFAULT ARGUMENTS ---
# Define a function to normalize a vector (scale values between 0 and 1)
normalize_vector <- function(x, na.rm = TRUE) {
  if (!is.numeric(x)) {
    stop("Input vector 'x' must be numeric.")
  }
  if (na.rm) {
    x <- na.omit(x)
  }
  # Compute minimum and maximum values
  x_min <- min(x)
  x_max <- max(x)
  # Normalize each element
  normalized <- (x - x_min) / (x_max - x_min)
  return(normalized)
}

# Test the normalization function on a synthetic vector
synthetic_vector <- c(5, 10, 15, 20, 25)
normalized_vector <- normalize_vector(synthetic_vector)
cat("Normalized vector:\n")
print(normalized_vector)

# --- FUNCTION UTILIZING APPLY ---
# Create a function that takes an expression data frame and returns normalized expression
# for each gene across all samples using row-wise normalization.
normalize_expression <- function(expr_data) {
  # Use apply to normalize each gene (each row)
  normalized_data <- t(apply(expr_data, 1, normalize_vector))
  return(normalized_data)
}

# Apply the function on the feature_counts dataset
normalized_feature_counts <- normalize_expression(feature_counts)
cat("Normalized feature counts (first few rows):\n")
print(head(normalized_feature_counts))
heatmap(as.matrix(feature_counts[1:100, 1:10]))

test_differential_expression <- function(expr_data, groups) {
  expr_data <- as.matrix(expr_data)
  group_means <- apply(expr_data, 1, function(x) {
    tapply(x, groups, mean)
  })
  # t test for each gene
  p_value <- apply(expr_data,1, function(x) {
    res <- t.test(x ~ groups)
    return(res$p.value)
  })
 results <- data.frame(t(group_means), p_value)
}

DE_results <- test_differential_expression(feature_counts, feature_annotation$gender)
head(DE_results)

# make volcano plot
plot(
  log2(DE_results$male/DE_results$female),
  -log10(DE_results$p_value),
  xlab = "-log10(p-value)", ylab = "log2(fold change)",
  main = "Volcano plot of differential expression",
  col = "#00000060"
)
