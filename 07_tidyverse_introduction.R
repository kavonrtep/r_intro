# ---
# title: "Data manipulation and cleaning with Tidyverse"
# format:
#   revealjs:
#     self-contained: true
# editor: visual
# ---
# <style>
# .reveal {
#   font-size: 160%;
# }
# </style>

# ## Tidyverse
# - **Tidyverse** is a collection of R packages designed for data science.
# - It provides a consistent set of functions and syntax for data manipulation and visualization.
# - Key packages include `dplyr`, `ggplot2`, `tidyr`, and `readr`.
#   - `dplyr`: data manipulation
#   - `ggplot2`: data visualization
#   - `tidyr`: data cleaning
#   - `readr`: data import
# - Use `install.packages("tidyverse")` to install the package.
# - Use `library(tidyverse)` to load the packages.

# ## Tidyverse `reader` package
# - The `readr` package provides functions for reading data into R.
# - The functions are similar to base R functions like `read.csv()`, but faster and more user-friendly.
# - Common functions:
#   - `read_csv()`: reads comma-separated files.
#   - `read_tsv()`: reads tab-separated files.
#   - `read_delim()`: reads files with custom delimiters.
#   - `read_fwf()`: reads fixed-width files.
#   - `read_table()`: reads tabular data.
# - Returns a `tibble`, a modern version of a data frame.

# ## Sample data
library(tidyverse)

# 1. Basic Data Inspection
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

# ## Sample data (cont'd)
# Inspect dimensions (number of rows and columns)
print(dim(df))  # e.g., 4 rows, 6 columns

# Print column names
print(colnames(df))

# Show the first few rows of the dataset
print(head(df))

# ## Using the Pipe Operator: `%>%`
# - The `%>%` operator passes the result of one function into the next.
s <- rnorm(100)
m <- mean(s)

# Using the pipe operator
m <- s %>% mean()
# or
m <- rnorm(100) %>% mean()
# is equivalent to
m <- mean(rnorm(100))

# ## Using the Pipe Operator (%>%)
# - The `%>%` operator passes the result of one function into the next.
# - Here, we select the "Gene Description" column and then display its head.
df %>%
  select(`Gene Description`) %>%
  head()

# This is equivalent to
head(select(df, `Gene Description`))

# ## `select()` function (dplyr package)
# - The `select()` function is used to select columns from a data frame.
# - It can be used to reorder columns or to exclude columns.
# - Syntax: `select(data, column1, column2, ...)`
# - Use `:` to select a range of columns.
# - Use `-` to exclude columns.
# - Use `everything()` to include all columns.
# - Use `starts_with()`, `ends_with()`, `contains()`, `matches()`, `one_of()` to select
#   columns based on name patterns.
# - See `?select` for more details.

# ## `select()` function examples
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

# ## Filtering rows with `filter()` function
# - `filter()` is used to select rows based on conditions.
# - Syntax: `filter(data, condition)`
# - Use logical operators like `==`, `!=`, `>`, `<`, `>=`, `<=`, `%in%`, `&`, `|`, `!`
#   to define conditions.
# - Use `between()` to filter rows within a range.
kinase_df <- df %>%
  filter(str_detect(`Gene Description`, "kinase"))
print(kinase_df)

# ## Creating new columns with `mutate()` function
# - `mutate()` is used to create new columns or modify existing columns.
# - Syntax: `mutate(data, new_column = expression)`
df_extended <- df %>%
  mutate(IsControl = str_detect(`Gene Accession Number`, "^AFFX"))
print(df_extended)

# ## Using `mutate()` to create new columns
# - Calculate the average expression across Sample1 and Sample2.
# - We select only the numeric sample columns.
df_extended <- df_extended %>%
  mutate(AvgExpr = rowMeans(select(., Sample1, Sample2), na.rm = TRUE))
print(df_extended)

# ## Grouping and summarizing with `group_by()` and `summarize()`
# - `group_by()` is used to group data by one or more variables.
# - Syntax: `group_by(data, column1, column2, ...)`
# - `summarize()` is used to calculate summary statistics for each group.
# - Syntax: `summarize(data, new_column = function(column))`
# - Un-group data using `ungroup()`.
by_cyl <- mtcars %>% group_by(cyl)
# Grouping does not change the data frame (it just adds a grouping attribute)
print(by_cyl)

# ## Grouping changes the behavior of functions
# - When a data frame is grouped, functions like `summarize()` calculate summary
#   statistics for each group.
# - Here, we calculate the mean `mpg` for each number of cylinders.
by_cyl <- mtcars %>% group_by(cyl) %>%
  summarize(mean_mpg = mean(mpg))
print(by_cyl)
by_cyl_gear <- mtcars %>% group_by(cyl, gear) %>%
  summarize(mean_mpg = mean(mpg))
print(by_cyl_gear)

# ## Grouping (cont'd)
control_summary <- df_extended %>%
  group_by(IsControl) %>%
  summarize(mean_expr = mean(AvgExpr), sd_expr = sd(AvgExpr))
print(control_summary)

# ## Reshaping data with `pivot_longer()` and `pivot_wider()` functions
# - `pivot_longer()` is used to convert wide data to long data.
# - `pivot_wider()` is used to convert long data to wide data.
# - Wide data has multiple columns for each variable.
# - Long data has multiple rows for each variable.

# ## Wide format example
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

# ## Long format example
# - Function `pivot_longer()` converts wide data to long data.
# - Syntax: `pivot_longer(data, cols, names_to, values_to)`
# - `cols`: columns to convert.
# - `names_to`: name of the new column that stores original column names.
# - `values_to`: name of the new column that stores the values.
long_df <- wide_df %>%
  pivot_longer(cols = starts_with("Expr"),
               names_to = "Sample",
               values_to = "Expression")
print(long_df)

# ## Converting long data to wide data
# - Function `pivot_wider()` converts long data to wide data.
# - Syntax: `pivot_wider(data, names_from, values_from)`
# - `names_from`: column to create new columns.
# - `values_from`: column to fill new columns.
wide_df <- long_df %>%
  pivot_wider(names_from = "Sample",
              values_from = "Expression")
print(wide_df)


## Data joining with `left_join()`, `right_join()`, `inner_join()`, and `full_join()`
# - `left_join()`: returns all rows from the left table and the matched rows from the right table.
# - `right_join()`: returns all rows from the right table and the matched rows from the left table.
# - `inner_join()`: returns only the rows that have matching values in both tables.
# - `full_join()`: returns all rows when there is a match in either table.

df1 <- tibble(ID = c(1, 2, 3), Value = c(10, 20, 30))
df2 <- tibble(ID = c(2, 3, 4), Value = c(200, 300, 400))
joined_df <- df1 %>% left_join(df2, by = "ID")
print(joined_df)

joined_df <- df1 %>% right_join(df2, by = "ID")
print(joined_df)

df1 <- tibble(ID = c(1, 2, 3), Value = c(10, 20, 30))
df2 <- tibble(ID = c(2, 3, 4), Value = c(200, 300, 400))

joined_df <- df1 %>% inner_join(df2, by = "ID")
print(joined_df)

joined_df <- df1 %>% full_join(df2, by = "ID")
print(joined_df)








# ## String manipulation with functions from `stringr` package
# - The `stringr` package provides functions for string manipulation.
# - Common functions include:
#   - `str_detect()`: detect patterns in strings.
#   - `str_replace()`: replace patterns in strings.
#   - `str_split()`: split strings based on patterns.
#   - `str_sub()`: extract substrings.
#   - `str_trim()`: remove leading/trailing whitespace.
#   - `str_to_lower()`, `str_to_upper()`: convert strings to lower or upper case.
#   - `str_c()`: concatenate strings.
#   - `str_remove()`, `str_remove_all()`: remove patterns from strings.

# ## `stringr` package examples
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

# ## Regular expressions : defining search patterns
# - Basic regex syntax:
#   - `.`: matches any character.
#   - `^`: matches the start of a string.
#   - `$`: matches the end of a string.
#   - `[]`: matches any character inside the brackets.
#   - `|`: matches either the expression before or after the pipe.
#   - `*`: matches zero or more occurrences of the preceding character.
#   - `+`: matches one or more occurrences of the preceding character.
#   - `?`: matches zero or one occurrence of the preceding character.
#   - `{n}`: matches exactly n occurrences.
#   - `{n,}`: matches n or more occurrences.
#   - `{n,m}`: matches between n and m occurrences.
#   - `\\`: escape character.
#   - `\\d`: matches any digit.
#   - `\\w`: matches any word character.

# ## Regular expressions examples
# 1. Dot (.)
# Matches any character except a newline.
str_detect("cat", "c.t")      # TRUE: 'c' + any character ('a') + 't'
str_detect("ct", "c.t")     # FALSE: missing a character between 'c' and 't'

# 2. Start (^) and End ($)
# ^ matches the start of a string, $ matches the end.
str_detect("Hello world", "^Hello")  # TRUE: string starts with "Hello"
str_detect("Hello world", "world$")   # TRUE: string ends with "world"

# 3. Character Class ([])
# Matches any single character listed inside the brackets.
str_extract("hello", "[aeiou]")       # "e": extracts the first vowel found
str_extract_all("hello", "[aeiou]")

# ## Regular expressions examples (cont'd)
# 4. Alternation (|)
# Matches either the expression before or after the pipe.
str_detect("I like cats", "cat|dog")    # TRUE: "cat" is found
str_detect("I like dogs", "cat|dog")    # TRUE: "dog" is found
str_detect("I like birds", "cat|dog")   # FALSE: neither "cat" nor "dog" is found

# 5. Asterisk (*)
# Matches zero or more occurrences of the preceding element.
str_detect("ac", "ab*c")    # TRUE: zero 'b's between 'a' and 'c'
str_detect("abbc", "ab*c")  # TRUE: one or more 'b's between 'a' and 'c'

# 6. Plus (+)
# Matches one or more occurrences of the preceding element.
str_detect("abc", "ab+c")   # TRUE: one 'b' present
str_detect("abbc", "ab+c")  # TRUE: two 'b's present
str_detect("ac", "ab+c")    # FALSE: requires at least one 'b'

# ## Regular expressions examples (cont'd)
# 7. Question Mark (?)
# Matches zero or one occurrence of the preceding element.
# Example: "colou?r" matches both "color" and "colour"
str_detect("color", "colou?r")    # TRUE
str_detect("colour", "colou?r")   # TRUE

# 8. Exact Count ({n})
# Matches exactly n occurrences of the preceding element.
str_detect("aaa", "a{3}")   # TRUE: exactly three 'a's
str_detect("aa", "a{3}")    # FALSE: only two 'a's

# 9. Minimum Count ({n,})
# Matches n or more occurrences.
str_detect("aa", "a{2,}")    # TRUE: at least two 'a's present
str_detect("a", "a{2,}")     # FALSE: only one 'a'

# ## Regular expressions examples (cont'd)
# 10. Range Count ({n,m})
# Matches between n and m occurrences (inclusive).
str_detect("aa", "a{2,4}")    # TRUE: two 'a's is within the range
str_detect("aaaa", "a{2,4}")  # TRUE: four 'a's is within the range
str_detect("aaaaa", "a{2,4}") # FALSE: five 'a's exceeds the range

# 11. Escape Character (\\)
# Use double backslashes to escape special characters.
# For example, to match a literal period (.) you must escape it.
str_detect("a.b", "a\\.b")   # TRUE: matches literal "a.b"
str_detect("aXb", "a\\.b")   # FALSE: "X" is not a literal period

# ## Regular expressions examples (cont'd)
# 12. Digit (\\d)
# Matches any digit (0-9).
str_extract("Order number: 12345", "\\d+")  # "12345": extracts one or more digits

# 13. Word Character (\\w)
# Matches any word character (letters, digits, underscore).
str_extract_all("Hello, world!", "\\w+")   # c("Hello", "world"): extracts words

# ## Testing regular expressions
# - Use online tools like https://regex101.com/ to test regular expressions.
