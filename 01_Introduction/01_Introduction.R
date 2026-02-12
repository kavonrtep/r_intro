################################################################################
# R SCRIPT FOR SESSION 1: INTRODUCTION TO R FUNDAMENTALS & BEST PRACTICES
#
# This script covers:
#
# 1. RStudio Interface & Workflow
# 2. Language Fundamentals
#    - Basic syntax (expressions, statements, operators)
#    - Variables and assignment (<-, =)
#    - Naming conventions for variables and functions
#    - Basic mathematical operations
#
# 3. Data Structures Overview
#    - Atomic vectors (numeric, character, logical, etc.) and subsetting
#    - Lists and nested lists
#    - Factors for categorical data
#    - Matrices and arrays (including indexing and dimension reduction)
#    - Data frames vs. tibbles: similarities and differences
#
# 4. Vectorized Operations & Sequence Generation
#    - Element-wise operations on vectors and matrices
#    - Generating sequences with :, seq(), seq_along(), seq_len(), rep()
#
# 5. Finding Help and Documentation
#    - Using ? and help() to find documentation
#    - Searching topics with ??
#
# 6. Working with Packages
#    - Installing packages (install.packages, BiocManager for Bioconductor)
#    - Loading and updating packages (library, require, update.packages)
#
# 7. Additional Topics & Best Practices
#    - Working with files and directories (getwd, setwd)
#    - Basic plotting, debugging tips, and code organization
#
# INSTRUCTIONS:
#   1. Open this script in RStudio.
#   2. Run each code section step by step to see outputs in the console.
#   3. Look for "TASK" comments and complete the mini-exercises.
#
# Suggestions for additional topics in future sessions:
#   - Control Structures: if/else, for loops, while loops, and vectorized alternatives.
#   - Writing and Debugging Functions: creating, testing, and documenting your own functions.
#   - Data Import/Export: reading from files (CSV, Excel) and writing output.
#   - Advanced Data Manipulation: an introduction to dplyr, tidyr, and other tidyverse packages.
################################################################################

########################### 0. RSTUDIO & WORKFLOW ###############################
# A. RStudio Interface & Workflow:
#    - Console: where commands are executed.
#    - Script Editor: where you write and save your code.
#    - Environment/History Pane: lists variables and past commands.
#    - Plots/Files/Help Pane: displays plots, files, and help documentation.
#
# B. Working with Files:
#    - Use getwd() to see your current working directory.
#    - Use setwd("your/path/here") to change your working directory.
#
# TASK (Optional):
#   1. Run getwd() to check your current directory.
#   2. (If desired) Change it using setwd().

getwd()
# setwd("C:/your/path/here")  # Uncomment and modify this line to change your working directory.
# Alternatively, use the RStudio interface to set the working directory from file explorer.

######################## 1. LANGUAGE FUNDAMENTALS ##############################

# --- 1.1 Basic Syntax, Operators, and Statements ---
2 + 2                     # Simple arithmetic expression.
print(2 + 2)              # Using print() to display the result.

# TASK 1: Modify the expression below to multiply 8 by 5 instead of 8 by 4.
8 * 4  # Modify this expression accordingly.

# --- 1.2 Variables, Assignment & Naming Conventions ---
x <- 10        # Preferred assignment operator in R.
y = 5          # Also valid, but the "<-" operator is more common in the R community.

# Best practices:
# - Use descriptive names (e.g., total_score instead of ts).
# - Use lowercase letters and underscores to separate words.
# - Function names often use a verb-noun pattern (e.g., calculate_mean).

z <- x + y     # Variables can be used in arithmetic operations.
z              # Simply typing the variable name prints its value.

# TASK 2: Create a variable named "height" with the value 1.85 and print it.

# --- 1.3 Basic Mathematical Operations ---
sum_xy   <- x + y   # Addition.
diff_xy  <- x - y   # Subtraction.
prod_xy  <- x * y   # Multiplication.
quot_xy  <- x / y   # Division.
expo_xy  <- x^2     # Exponentiation.

print(sum_xy)
print(diff_xy)
print(prod_xy)
print(quot_xy)
print(expo_xy)

# Using built-in functions:
sqrt_x <- sqrt(x)   # Square root of x.
log_x  <- log(x)    # Natural logarithm of x.
exp_x  <- exp(y)    # Exponential function applied to y.

# TASK 3: Print sqrt_x, log_x, and exp_x. Then, try changing x or y and re-run.

######################## 2. DATA STRUCTURES OVERVIEW ###########################

# --- 2.1 Atomic Vectors ---
# Vectors are the simplest R data structures and can store elements of a single type.
numeric_vec <- c(1, 2, 3, 4)                      # Numeric vector.
char_vec    <- c("apple", "banana", "cherry", "cherry")  # Character vector.
log_vec     <- c(TRUE, FALSE, TRUE)               # Logical vector.

print(numeric_vec)
print(char_vec)
print(log_vec)

# TASK 4: Create a vector "my_vec" with three numbers of your choice.

# --- 2.1.1 Indexing and Subsetting Vectors ---
# R uses 1-based indexing.
numeric_vec[1]          # First element of numeric_vec.
numeric_vec[1:3]        # Elements from first to third.
numeric_vec[c(1, 4)]    # First and fourth elements.

# Logical indexing:
is_large <- numeric_vec > 2
is_large                # Logical vector indicating which elements are greater than 2.
numeric_vec[is_large]   # Returns elements where the condition is TRUE.
numeric_vec[numeric_vec > 2]  # Concise logical indexing.
numeric_vec[-2]         # Returns all elements except the second.

# Additional functions for accessing elements:
head(numeric_vec)       # Returns the first few elements of numeric_vec.
tail(numeric_vec)       # Returns the last few elements of numeric_vec.

# --- 2.2 Lists and Nested Lists ---
# Lists can store elements of different types.
my_list <- list(name = "Sample", numbers = numeric_vec, flag = TRUE)
print(my_list)
print(my_list$name)            # Access element using the $ operator.
print(my_list[["numbers"]])    # Alternative access using double brackets.

# TASK 5: Add a new element "my_char" to my_list that stores char_vec.

# --- 2.3 Factors ---
# Factors represent categorical data with predefined levels.
fruit_factor <- factor(char_vec)
print(levels(fruit_factor))     # Unique levels (categories) in the factor.
print(fruit_factor)
print(as.integer(fruit_factor)) # Underlying integer codes.

# TASK 6: Create a factor "condition_factor" with levels "Control" and "Treatment"
# and assign it a small vector (e.g., containing "Control", "Treatment", "Control").

# --- 2.4 Matrices and Arrays ---
# Matrices are 2-dimensional arrays; arrays can have more dimensions.
my_matrix <- matrix(1:9, nrow = 3, ncol = 3)
print(my_matrix)

# Creating a matrix filled by rows:
my_matrix2 <- matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE)
print(my_matrix2)

# Accessing elements by row and column indices:
print(my_matrix[1, 2])  # Element at row 1, column 2.
print(my_matrix[2, ])   # Entire second row.
print(my_matrix[, 3])   # Entire third column.

# Dimension reduction: R returns the simplest possible structure by default.
row2 <- my_matrix[2, ]
row2
# Use drop = FALSE to preserve matrix structure:
row2_as_matrix <- my_matrix[2, , drop = FALSE]
row2_as_matrix

# TASK 7: Store the second row of my_matrix in a variable and then print
# the element at row 2, column 3. Experiment with drop = FALSE.

# Matrix transposition:
transposed_matrix <- t(my_matrix)
print(my_matrix)
print(transposed_matrix)

# --- 2.5 Data Frames and Tibbles ---
# Data frames are lists of equal-length vectors (similar to spreadsheets).
my_df <- data.frame(
  ID       = 1:4,
  fruit    = c("apple", "banana", "cherry", "cherry"),
  quantity = c(10, 5, 7, 5)
)
print(my_df)

# Access columns in various ways:
print(my_df$fruit)         # Using the $ operator.
print(my_df[["fruit"]])    # Using double brackets.
print(my_df[, "fruit"])    # Using the column name.
print(my_df[, 2])          # Using the column index.

# Adding a new column to the data frame:
my_df$price <- c(2.5, 1.8, 3.2, 3.5)
print(my_df)

# subsetting data frames
print(my_df[my_df$quantity > 5, ])
print(subset(my_df, quantity > 5))


# TASK 8: Use str(my_df) to view the structure of the data frame;
# then, print the "fruit" column using the $ operator.

# TASK 8.1 Select all rows where fruit is "cherry" and print the result.

################### 3. VECTORIZED OPERATIONS & SEQUENCE GENERATION ###############

# Vectorized operations allow operations on whole vectors without explicit loops.
# Example: measurements for two conditions.
condition_A <- c(2.1, 3.4, 6.5, 7.9, 10.2)
condition_B <- c(1.5, 2.3, 5.9, 8.0, 9.8)

# Element-wise arithmetic operations:
sum_AB   <- condition_A + condition_B   # Sum corresponding elements.
diff_AB  <- condition_A - condition_B   # Difference for each element.
ratio_AB <- condition_A / condition_B   # Ratio of corresponding elements.

print(sum_AB)
print(diff_AB)
print(ratio_AB)

# TASK 9: Multiply condition_A and condition_B element-wise and print the result.

# Calculate summary statistics:
mean_A <- mean(condition_A)
sd_A   <- sd(condition_A)
print(mean_A)
print(sd_A)

mean_sum  <- mean(sum_AB)
mean_diff <- mean(diff_AB)

# TASK 10: Print mean_sum and mean_diff.
# Also, calculate the mean for condition_B and compare the results.

# --- Matrix Operations ---
set.seed(123)  # For reproducibility.
count_matrix <- matrix(sample(10:100, 12, replace = TRUE), nrow = 3, ncol = 4)
print(count_matrix)

# Arithmetic on matrices:
print(count_matrix + 5)   # Add 5 to every element.
print(count_matrix * 2)   # Multiply every element by 2.

another_matrix <- matrix(sample(10:100, 12, replace = TRUE), nrow = 3, ncol = 4)
print(count_matrix + another_matrix)

# Recycling example:
print(count_matrix + c(0.1, 0.2, 0.3))
print(count_matrix + c(0.1, 0.2, 0.3, 0.4))
print(count_matrix + c(0.1, 0.2, 0.3, 0.4, 0.5))

# TASK 11: Subtract another_matrix from count_matrix and print the result.
# Identify which elements become negative and consider why.

# Row-wise and Column-wise operations:
print(rowSums(count_matrix))
print(colSums(count_matrix))
print(rowMeans(count_matrix))
print(colMeans(count_matrix))

# TASK 12: Save rowSums to a new variable and print it.
# Also, save colSums, rowMeans, and colMeans to separate variables and print them.

# --- Matrix Multiplication (Optional) ---
mat1 <- count_matrix          # 3x4 matrix.
mat2 <- matrix(1:8, nrow = 4, ncol = 2)  # 4x2 matrix.
mat_mult <- mat1 %*% mat2
print(mat_mult)
# TASK 13 (Challenge): Print mat_mult and explain in comments why its dimensions are 3x2.

# --- Generating Sequences ---
# Colon operator: creates a sequence of integers.
sequence1 <- 1:10
print(sequence1)
# TASK 14: Create and print a sequence from 5 to 15 using the colon operator.

# seq() function: allows more control over sequence generation.
sequence2 <- seq(from = 0, to = 1, by = 0.2)
print(sequence2)
sequence3 <- seq(from = 0, to = 1, length.out = 6)
print(sequence3)
# TASK 15: Generate a sequence from 10 to 20 with a step size of 2 using seq().

# seq_along() function: generates a sequence along the length of an object.
example_vector <- c("a", "b", "c", "d")
sequence4 <- seq_along(example_vector)
print(sequence4)
# TASK 16: Create a vector of 5 colors and use seq_along() to generate their indices.

# seq_len() function: creates a sequence of a specified length.
sequence5 <- seq_len(7)
print(sequence5)
# TASK 17: Create a sequence of length 10 using seq_len() and print it.

# --- Using the rep() Function ---
# The rep() function replicates the values in a vector.
repeated_seq <- rep(1:3, times = 3)  # Replicates 1, 2, 3 three times.
print(repeated_seq)
# TASK (Optional): Modify the rep() function to create a different pattern.

######################## 4. FINDING HELP & DOCUMENTATION #######################

# --- Basic Help ---
?mean        # Opens the help page for the mean() function.
help(mean)   # Alternative method to access documentation.

# --- Searching for Help with ?? ---
# ??anova   # Searches for help related to analysis of variance.
# TASK 18: Run ?plot and read the documentation in the Help pane.
# TASK 19: Try using ??plot (or ??factor) to explore related topics.

######################## 5. WORKING WITH PACKAGES ###############################

# --- Installing Packages ---
# install.packages("tidyverse")
#
# For Bioconductor packages:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GenomicRanges")

# --- Loading and Updating Packages ---
# library(tidyverse)
# update.packages()

# --- Checking Package Libraries ---
print(.libPaths())

# TASK 20: (Requires an internet connection)
#   1. Uncomment the install.packages() call below to install "stringr".
#   2. Load it using library(stringr) and check its help page using ?str_detect.
# install.packages("stringr")
# library(stringr)
# ?str_detect

#################### 6. ADDITIONAL TOPICS & BEST PRACTICES #######################

# --- Working with Files & Directories ---
print(getwd())

# --- Basic Plotting ---
x_values <- 1:40
y_values <- x_values^2  # Squaring each x value.
plot(x_values, y_values, type = "b", col = "blue",
     main = "Simple Plot of x vs. x^2",
     xlab = "x", ylab = "x squared")
# TASK 21: Modify the plot by changing the color and type to see how the plot changes.

# --- Debugging Tips ---
# When you encounter errors:
#   - Read the error message carefully.
#   - Use traceback() immediately after an error to see the call stack.
#   - Use browser() inside functions to debug step-by-step.
#
# Example (uncomment to test):
# faulty_function <- function() {
#   stop("This is an intentional error for demonstration.")
# }
# faulty_function()
# traceback()


################################################################################
#        Geting data from Excel or LibreOffice Calc to R
################################################################################
# Open file data/toothgrowth.xlsx in Libre Office Calc
# this file contains data on tooth growth in guinea pigs
# under different supplementations # with vitamin C
# in form of Orange Juice or Ascorbic Acid


# save it in text/csv format, select field delimiter. It could be
# comma or Tab delimited, In this example we will use Tab delimited:
# import the data into R
# Read the tab-delimited file into a data frame.
# - header = TRUE tells R that the first row contains column names.
# - sep = "\t" specifies that the fields are separated by tabs.
tooth_growth <- read.table("exported_data.txt",
                   header = TRUE,    # First row has column names.
                   sep = "\t",       # Tab-delimited file.
                   stringsAsFactors = FALSE)  # Do not convert text to factors.

# Check the first few rows to ensure the data is imported correctly.
head(tooth_growth)
str(tooth_growth)
summary(tooth_growth)
# plot dose vs length
# Create a scatter plot with:
# - x-axis: dose
# - y-axis: length
plot(tooth_growth$dose, tooth_growth$len,
     xlab = "Dose (mg/day)",
     ylab = "Tooth Length",
     main = "Tooth Growth in Guinea Pigs")








# Plotting demonstration 1: see script 01_example_1.R for a detailed example.

# Plotting demonstration 2: see script 01_example_2.R for a detailed example.




