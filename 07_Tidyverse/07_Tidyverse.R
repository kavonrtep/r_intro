################################################################################
# R SCRIPT FOR SESSION 7: DATA MANIPULATION WITH TIDYVERSE
#
# This script covers:
# 1. Tidyverse overview and the readr package
# 2. The pipe operator %>%
# 3. Selecting columns with select()
# 4. Filtering rows with filter()
# 5. Creating new columns with mutate()
# 6. Grouping and summarizing with group_by() and summarize()
# 7. Reshaping data with pivot_longer() and pivot_wider()
# 8. Joining data frames
# 9. String manipulation with stringr
# 10. Regular expressions
#
# DATA SOURCES USED:
# - mtcars: built-in R dataset (car performance data)
# - iris: built-in R dataset (flower measurements)
# - Small tibble: simulated gene expression data (created inline)
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################

library(tidyverse)

# ---------- SLIDE: Tidyverse ----------

################################################################################
# SECTION 1: TIDYVERSE AND READR
################################################################################

# ---------- SLIDE: Tidyverse `reader` package ----------

# readr provides fast, user-friendly data import functions.
# All return a tibble — a modern data frame.
#
# read_csv("file.csv")         # comma-separated
# read_tsv("file.tsv")         # tab-separated
# read_delim("file.txt", ";")  # custom delimiter
# read_fwf("file.txt")         # fixed-width
# read_table("file.txt")       # whitespace-separated

# ---------- SLIDE: Sample data ----------

# Create a small tibble simulating a gene expression dataset.
df <- tibble(
  `Gene Description` = c("Protein kinase C (alpha)",
                         "Ribosomal protein L2",
                         "ATP synthase", "Control gene"),
  `Gene Accession Number` = c("Gene1", "Gene2", "Gene3", "AFFX-Control"),
  Sample1 = c(5.2, 3.4, 7.8, 1.2),
  call1 = c("P", "A", "P", "A"),
  Sample2 = c(5.5, 3.8, 8.0, 1.3),
  call2 = c("P", "A", "P", "A")
)

# ---------- SLIDE: Sample data (cont'd) ----------

# Inspect dimensions (number of rows and columns)
print(dim(df))       # 4 rows, 6 columns

# Print column names
print(colnames(df))

# Show the first few rows of the dataset
print(head(df))

################################################################################
# SECTION 2: THE PIPE OPERATOR
################################################################################

# ---------- SLIDE: Using the Pipe Operator: `%>%` ----------

# The %>% operator passes the result of one function into the next.
s <- rnorm(100)
m <- mean(s)

# Using the pipe operator
m <- s %>% mean()
# or
m <- rnorm(100) %>% mean()
# is equivalent to
m <- mean(rnorm(100))

# ---------- SLIDE: Using the Pipe Operator (%>%) ----------

# Here, we select the "Gene Description" column and then display its head.
df %>%
  select(`Gene Description`) %>%
  head()

# This is equivalent to
head(select(df, `Gene Description`))

################################################################################
# SECTION 3: SELECT
################################################################################

# ---------- SLIDE: `select()` function (dplyr package) ----------

# select() picks columns from a data frame.
# Key helpers:
#   starts_with(), ends_with(), contains(), matches(), one_of()
#   everything()  — include all remaining columns
#   -column       — exclude a column
# See ?select for more details.

# ---------- SLIDE: `select()` function examples ----------

# Select columns by name
df %>%
  select(`Gene Description`, `Gene Accession Number`) %>%
  head()

# Remove columns using `-`
df %>%
  select(-`Gene Description`, -`Gene Accession Number`) %>%
  head()

# Remove columns using `-` and contains()
df %>%
  select(-contains("call")) %>%
  head()

# TASK 1:
# Remove only the columns whose names start with "call" using starts_with().
# Hint: select(-starts_with("call"))


################################################################################
# SECTION 4: FILTER
################################################################################

# ---------- SLIDE: Filtering rows with `filter()` function ----------

# filter() keeps rows that match a condition.
# Logical operators: ==, !=, >, <, >=, <=, %in%, &, |, !
# Use between(x, lo, hi) to filter within a range.

kinase_df <- df %>%
  filter(str_detect(`Gene Description`, "kinase"))
print(kinase_df)

# TASK 2:
# Filter the df tibble to keep only rows where Sample1 > 4.
# Hint: filter(Sample1 > 4)


################################################################################
# SECTION 5: MUTATE
################################################################################

# ---------- SLIDE: Creating new columns with `mutate()` function ----------

# mutate() adds new columns or modifies existing ones.
# Syntax: mutate(data, new_column = expression)

df_extended <- df %>%
  mutate(IsControl = str_detect(`Gene Accession Number`, "^AFFX"))
print(df_extended)

# ---------- SLIDE: Using `mutate()` to create new columns ----------

# Calculate the average expression across Sample1 and Sample2.
df_extended <- df_extended %>%
  mutate(AvgExpr = rowMeans(select(., Sample1, Sample2), na.rm = TRUE))
print(df_extended)

# TASK 3:
# Add a column "HighExpr" that is TRUE when AvgExpr > 5.
# Hint: mutate(HighExpr = AvgExpr > 5)


################################################################################
# SECTION 6: GROUP_BY AND SUMMARIZE
################################################################################

# ---------- SLIDE: Grouping and summarizing with `group_by()` and `summarize()` ----------

# group_by() adds grouping metadata — the data frame looks the same,
# but subsequent operations are applied per group.
# Use ungroup() to remove grouping.

by_cyl <- mtcars %>% group_by(cyl)
# Grouping does not change the data frame (it just adds a grouping attribute)
print(by_cyl)

# ---------- SLIDE: Grouping changes the behavior of functions ----------

# summarize() calculates summary statistics for each group.
by_cyl <- mtcars %>% group_by(cyl) %>%
  summarize(mean_mpg = mean(mpg))
print(by_cyl)
by_cyl_gear <- mtcars %>% group_by(cyl, gear) %>%
  summarize(mean_mpg = mean(mpg))
print(by_cyl_gear)

# ---------- SLIDE: Grouping (cont'd) ----------

control_summary <- df_extended %>%
  group_by(IsControl) %>%
  summarize(mean_expr = mean(AvgExpr), sd_expr = sd(AvgExpr))
print(control_summary)

# TASK 4:
# Using mtcars, calculate the median horsepower (hp) for each number of cylinders.
# Hint: mtcars %>% group_by(cyl) %>% summarize(median_hp = median(hp))


################################################################################
# SECTION 7: RESHAPING DATA
################################################################################

# ---------- SLIDE: Reshaping data with `pivot_longer()` and `pivot_wider()` functions ----------

# pivot_longer(): wide → long (many columns → key-value pairs)
# pivot_wider():  long → wide (key-value pairs → many columns)
# Wide: one row per gene, one column per sample
# Long: one row per gene-sample combination

# ---------- SLIDE: Wide format example ----------

wide_df <- tibble(
  ID = c(1, 2, 3, 4),
  Gene = c("Gene1", "Gene2", "Gene3", "Gene4"),
  Expr1 = c(5.2, 3.4, 7.8, 1.2),
  Expr2 = c(5.5, 3.8, 8.0, 1.3),
  Expr3 = c(5.5, 3.8, 8.0, 1.3),
  Expr4 = c(5.5, 3.8, 8.0, 1.3),
  Expr5 = c(5.5, 3.8, 8.0, 1.3)
)
print(wide_df)

# ---------- SLIDE: Long format example ----------

# pivot_longer() converts wide data to long data.
# cols:      columns to pivot
# names_to:  name of the new column storing original column names
# values_to: name of the new column storing values
long_df <- wide_df %>%
  pivot_longer(cols = starts_with("Expr"),
               names_to = "Sample",
               values_to = "Expression")
print(long_df)

# ---------- SLIDE: Converting long data to wide data ----------

# pivot_wider() converts long data back to wide data.
# names_from:  column whose values become new column names
# values_from: column whose values fill the new columns
wide_df <- long_df %>%
  pivot_wider(names_from = "Sample",
              values_from = "Expression")
print(wide_df)

# TASK 5:
# Start from long_df. Add a column "log2Expr" = log2(Expression), then
# pivot back to wide format so each "Expr*" column contains log2 values.
# Hint: mutate first, then pivot_wider()


################################################################################
# SECTION 8: JOINING DATA FRAMES
################################################################################

# ---------- SLIDE: Data joining with `left_join()`, `right_join()`, `inner_join()`, and `full_join()` ----------

# left_join():  all rows from the left table + matched rows from right
# right_join(): all rows from the right table + matched rows from left
# inner_join(): only rows with matching values in both tables
# full_join():  all rows from both tables

df1 <- tibble(ID = c(1, 2, 3), Value = c(10, 20, 30))
df2 <- tibble(ID = c(2, 3, 4), Value = c(200, 300, 400))
joined_df <- df1 %>% left_join(df2, by = "ID")
print(joined_df)

joined_df <- df1 %>% right_join(df2, by = "ID")
print(joined_df)

# ---------- SLIDE: Data joining (cont'd) ----------

df1 <- tibble(ID = c(1, 2, 3), Value = c(10, 20, 30))
df2 <- tibble(ID = c(2, 3, 4), Value = c(200, 300, 400))

joined_df <- df1 %>% inner_join(df2, by = "ID")
print(joined_df)

joined_df <- df1 %>% full_join(df2, by = "ID")
print(joined_df)

# TASK 6:
# Create a third tibble df3 with columns ID and Label:
#   df3 <- tibble(ID = c(1, 3, 5), Label = c("A", "B", "C"))
# Join df1 with df3 using full_join(). Which rows get NA? Why?


################################################################################
# SECTION 9: STRING MANIPULATION WITH STRINGR
################################################################################

# ---------- SLIDE: String manipulation with functions from `stringr` package ----------

# Common stringr functions:
#   str_detect()     — detect a pattern (returns TRUE/FALSE)
#   str_replace()    — replace first match
#   str_split()      — split string by pattern
#   str_sub()        — extract substring by position
#   str_trim()       — remove leading/trailing whitespace
#   str_to_lower(), str_to_upper() — change case
#   str_c()          — concatenate strings
#   str_remove(), str_remove_all() — remove pattern

# ---------- SLIDE: `stringr` package examples ----------

library(stringr)
# Detect patterns in strings
str_detect("Hello, world!", "world")
# Replace patterns in strings
str_replace("Hello, world!", "world", "universe")
# Split strings based on patterns
str_split("Hello, world!", ",")
# Extract substrings
str_sub("Hello, world!", 3, 7)
# Remove leading and trailing whitespace
str_trim("  Hello, world!  ")
# Remove patterns from strings
str_remove("Hello, world!", "world")

# TASK 7:
# Given the vector genes <- c("TP53_human", "BRCA1_human", "EGFR_mouse"),
# use str_remove() to strip the "_human" and "_mouse" suffixes.
# Hint: str_remove(genes, "_.*") removes everything from "_" onward.


################################################################################
# SECTION 10: REGULAR EXPRESSIONS
################################################################################

# ---------- SLIDE: Regular expressions : defining search patterns ----------

# Basic regex syntax:
#   .       — any character (except newline)
#   ^       — start of string
#   $       — end of string
#   []      — any character inside brackets
#   |       — alternation (OR)
#   *       — 0 or more occurrences
#   +       — 1 or more occurrences
#   ?       — 0 or 1 occurrence
#   {n}     — exactly n occurrences
#   {n,}    — n or more occurrences
#   {n,m}   — between n and m occurrences
#   \\      — escape character
#   \\d     — any digit
#   \\w     — any word character (letter, digit, underscore)

# ---------- SLIDE: Regular expressions examples ----------

# 1. Dot (.)
# Matches any character except a newline.
str_detect("cat", "c.t")       # TRUE: 'c' + any character ('a') + 't'
str_detect("ct", "c.t")        # FALSE: missing a character between 'c' and 't'

# 2. Start (^) and End ($)
# ^ matches the start of a string, $ matches the end.
str_detect("Hello world", "^Hello")  # TRUE: string starts with "Hello"
str_detect("Hello world", "world$")  # TRUE: string ends with "world"

# 3. Character Class ([])
# Matches any single character listed inside the brackets.
str_extract("hello", "[aeiou]")       # "e": extracts the first vowel found
str_extract_all("hello", "[aeiou]")

# ---------- SLIDE: Regular expressions examples (cont'd) ----------

# 4. Alternation (|)
str_detect("I like cats", "cat|dog")    # TRUE: "cat" is found
str_detect("I like dogs", "cat|dog")    # TRUE: "dog" is found
str_detect("I like birds", "cat|dog")   # FALSE: neither found

# 5. Asterisk (*)
# Matches zero or more occurrences of the preceding element.
str_detect("ac", "ab*c")    # TRUE: zero 'b's between 'a' and 'c'
str_detect("abbc", "ab*c")  # TRUE: two 'b's between 'a' and 'c'

# 6. Plus (+)
# Matches one or more occurrences of the preceding element.
str_detect("abc", "ab+c")   # TRUE: one 'b' present
str_detect("abbc", "ab+c")  # TRUE: two 'b's present
str_detect("ac", "ab+c")    # FALSE: requires at least one 'b'

# ---------- SLIDE: Regular expressions examples (cont'd) ----------

# 7. Question Mark (?)
str_detect("color", "colou?r")    # TRUE
str_detect("colour", "colou?r")   # TRUE

# 8. Exact Count ({n})
str_detect("aaa", "a{3}")   # TRUE: exactly three 'a's
str_detect("aa", "a{3}")    # FALSE: only two 'a's

# 9. Minimum Count ({n,})
str_detect("aa", "a{2,}")    # TRUE: at least two 'a's present
str_detect("a", "a{2,}")     # FALSE: only one 'a'

# ---------- SLIDE: Regular expressions examples (cont'd) ----------

# 10. Range Count ({n,m})
str_detect("aa", "a{2,4}")    # TRUE: two 'a's is within the range
str_detect("aaaa", "a{2,4}")  # TRUE: four 'a's is within the range
str_detect("aaaaa", "a{2,4}") # FALSE: five 'a's exceeds the range

# 11. Escape Character (\\)
str_detect("a.b", "a\\.b")   # TRUE: matches literal "a.b"
str_detect("aXb", "a\\.b")   # FALSE: "X" is not a literal period

# ---------- SLIDE: Regular expressions examples (cont'd) ----------

# 12. Digit (\\d)
str_extract("Order number: 12345", "\\d+")  # "12345"

# 13. Word Character (\\w)
str_extract_all("Hello, world!", "\\w+")    # c("Hello", "world")

# ---------- SLIDE: Testing regular expressions ----------

# Use online tools like https://regex101.com/ to test regular expressions.

# TASK 8:
# Given: accessions <- c("AFFX-BioB-3_at", "Gene_1234", "AFFX-BioC-5_at", "Kinase_99")
# Use str_detect() with a regex to extract only the AFFX control probes.
# Hint: str_detect(accessions, "^AFFX")