################################################################################
# R SCRIPT FOR SESSION 1: INTRODUCTION TO R FUNDAMENTALS
#
# This script covers:
#   1. Language Fundamentals
#       - Basic syntax (expressions, statements, operators)
#       - Variables: naming conventions, assignment operators (<-, =)
#       - Basic mathematical operations
#
#   2. Data Structures Overview
#       - Atomic vectors (numeric, character, logical, etc.)
#       - Indexing and subsetting (vectors, lists, matrices) with numeric or logical indices
#       - Lists and nested lists
#       - Factors and their role in categorical data
#       - Matrices and arrays (dimension reduction)
#       - Data frames vs. tibbles
#
#   3. Finding Help and Documentation
#       - Using ? and help()
#       - Vignettes and package manuals
#       - Online resources (CRAN, Bioconductor, RStudio Community)
#
#   4. Working with Packages
#       - Installing packages (install.packages, BiocManager)
#       - Loading and updating packages (library, require, update.packages)
#       - Package directories and version management
#
# INSTRUCTIONS:
#   1. Open this script in RStudio.
#   2. Run each code section step by step to see outputs in the console.
#   3. Look for "TASK" comments and complete the mini-exercises.
################################################################################

# ================================ 1. BASICS ===================================

# --- 1.1 Basic Syntax, Operators, and Statements ---

# Below is a simple arithmetic expression in R.
# You can select the line and press Ctrl + Enter (Windows) or Cmd + Enter (Mac)
# to run it in RStudio.

2 + 2

# You can also wrap expressions in the print() function to display them:
print(2 + 2)

# By default, R prints the result of any expression you evaluate in the console.
# There's no need for an explicit statement terminator like semicolons
# (though you can use them, it's not common in R).

# TASK 1: Change the expression below to multiply 8 by 5 instead of 8 by 4.
8 * 4


# --- 1.2 Variables and Assignment ---

# In R, you can assign values to variables using the "<-" or "=" operator.
x <- 10       # This is the most common way in R
y = 5         # This also works, but "<-" is preferred in most style guides

# Now we can operate on these variables:
z <- x + y
z

# Variable naming conventions:
# - Start with a letter (a-z, A-Z) or dot (.), but not a number.
# - Avoid using reserved words like if, else, for, in, etc.
# - Use meaningful names for readability.

# TASK 2: Create a new variable called "height" that stores the value 1.85.
# Then print it to the console.


# --- 1.3 Basic Mathematical Operations ---

# R can do a variety of math operations out of the box:
sum_xy <- x + y
diff_xy <- x - y
prod_xy <- x * y
quot_xy <- x / y
expo_xy <- x^2

print(sum_xy)   # Summation
print(diff_xy)  # Difference
print(prod_xy)  # Product
print(quot_xy)  # Division
print(expo_xy)  # Exponentiation

# R also has built-in mathematical functions like sqrt(), log(), exp(), sin(), etc.
sqrt_x <- sqrt(x)
log_x <- log(x)      # natural log by default
exp_x <- exp(y)

# TASK 3: Print out sqrt_x, log_x, and exp_x. Then try changing x or y
# to see how the results change.


# ================================ 2. DATA STRUCTURES ==========================

# --- 2.1 Atomic Vectors ---

# The simplest R data structure is a vector, which can contain elements of the
# same type (numeric, character, logical, etc.).

numeric_vec <- c(1, 2, 3, 4)                   # numeric vector
char_vec <- c("apple", "banana", "cherry")     # character vector
log_vec <- c(TRUE, FALSE, TRUE)                # logical vector

print(numeric_vec)
print(char_vec)
print(log_vec)

# TASK 4: Create a vector named my_vec that contains three numbers of your choice.


# --- 2.1.1 Indexing and Subsetting Vectors ---

# R indexing is 1-based (the first element is at index 1).
# Numeric indexing:
numeric_vec[1]          # first element
numeric_vec[1:3]        # first three elements
numeric_vec[c(1, 4)]    # the 1st and 4th elements

# Logical indexing:
# - If you provide a logical vector of the same length, TRUE indicates that
#   you want that element, FALSE indicates that you don't.
is_large <- numeric_vec > 2
is_large
numeric_vec[is_large]   # returns elements where numeric_vec > 2

# Negative indexing removes elements:
numeric_vec[-2]         # all elements except the 2nd


# --- 2.2 Lists and Nested Lists ---

# Lists can hold objects of different types and different lengths.
my_list <- list(name = "Sample", numbers = numeric_vec, flag = TRUE)
my_list

# You can access elements of a list with:
my_list$name        # $ operator
my_list[["numbers"]]

# Notice the difference between using single brackets [] vs. double brackets [[]]:
my_list[1]          # returns a sub-list containing the first element
my_list[[1]]        # returns the actual first element (e.g., "Sample")

# TASK 5: Add another item to my_list with the name "my_char" that stores char_vec.
# Hint: you can do something like: my_list$my_char <- char_vec


# --- 2.3 Factors ---

# Factors are used to represent categorical data (levels) in R.
# They store data as integers with labels.
fruit_factor <- factor(char_vec)
fruit_factor

# Check the internal integer representation:
as.integer(fruit_factor)

# Notice that "apple", "banana", and "cherry" are treated as distinct levels.

# TASK 6: Create a factor named condition_factor with levels "Control" and "Treatment",
# and assign it a small vector of values like c("Control", "Treatment", "Control").


# --- 2.4 Matrices and Arrays ---

# A matrix is a 2D data structure of the same type.
my_matrix <- matrix(1:9, nrow = 3, ncol = 3)
my_matrix

# You can also create arrays with more dimensions:
my_array <- array(1:12, dim = c(2, 2, 3))
my_array

# Indexing works similarly, but now we use row and column indices (or more):
my_matrix[1, 2]   # element in row 1, col 2
my_matrix[2, ]    # entire 2nd row
my_matrix[, 3]    # entire 3rd column

# Dimension reduction: When you subset rows or columns, R drops that dimension
# by default. For example:
row2 <- my_matrix[2, ]
# row2 is now a vector, not a matrix. If you want to preserve the 2D structure:
row2_as_matrix <- my_matrix[2, , drop = FALSE]

# TASK 7: Print the second row of my_matrix (store it in a new variable).
# Then print the element in row 2, column 3 of my_matrix.
# Try using drop = FALSE and see the difference.


# --- 2.5 Data Frames and Tibbles ---

# A data frame is a list of vectors of equal length (like a spreadsheet).
my_df <- data.frame(
  ID = 1:3,
  fruit = c("apple", "banana", "cherry"),
  quantity = c(10, 5, 7)
)
my_df

# Tibbles (from the tibble package) are a modern version of data frames
# with some nicer printing features, but structurally they're similar.

# To create a tibble, you need the tibble or tidyverse package installed.
# We'll do the installation demonstration in the "Working with Packages" section.
# For now, here's how you'd create one if the package is installed:
# library(tibble)
# my_tibble <- tibble(
#   ID = 1:3,
#   fruit = c("apple", "banana", "cherry"),
#   quantity = c(10, 5, 7)
# )
# my_tibble

# TASK 8: Print the structure of my_df using the str() function.
# Then use the $ operator to print the "fruit" column only.


# ====================== 3. FINDING HELP AND DOCUMENTATION =====================

# You can use ? or help() to get help on a function.
# For example, try:
?mean
# or
help(mean)

# Many packages include "vignettes", which are longer-form documentation.
# For example, if you have "dplyr" installed, you can see:
# vignette("dplyr")

# Online resources include:
# - CRAN: https://cran.r-project.org/
# - Bioconductor: https://bioconductor.org/
# - RStudio Community: https://community.rstudio.com/

# TASK 9: Run ?plot and read the documentation in the Help pane.

# In addition to ?function_name (e.g., ?mean), you can search for topics with ??.
# For instance:
# ??anova
#
# This will search for all help pages that mention "anova" in the title or
# description across installed packages.
#
# TASK 6:
#   1. Try using ??plot (or ??factor) to see what results come up.
#   2. Check some of the help pages that appear.


# ============================ 4. WORKING WITH PACKAGES ========================

# --- 4.1 Installing Packages ---

# To install packages from CRAN, use install.packages("package_name").
# For example, to install the "tidyverse" meta-package:
# install.packages("tidyverse")

# Bioconductor packages are installed via BiocManager:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GenomicRanges")

# NOTE: Installing packages requires an internet connection and may take time.

# --- 4.2 Loading and Updating Packages ---

# After installation, you load a package with library() or require().
# library(tidyverse)

# To update installed packages, use:
# update.packages()

# --- 4.3 Package Directories and Version Management ---

# R installs packages into a library directory. You can check which directories
# are used for libraries with .libPaths().

# TASK 10: If you have an internet connection, try installing the "stringr" package
# (part of the tidyverse) by removing the comment below and running it.
# install.packages("stringr")

# Then load it:
# library(stringr)

# TASK 11: Use ?str_detect (from stringr) to open the help and see how to use
# the function. Try a simple example.

################################################################################
# END OF SCRIPT
#
# Additional Tips:
#   - Keep scripts organized with headings and comments.
#   - Experiment with each example and see what happens if you change values.
#   - Refer to official documentation when you encounter a new function.
#
################################################################################
