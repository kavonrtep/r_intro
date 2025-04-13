## Follow-Up Assignment: Enhancing Assembly Statistics Calculation

### Overview

In this assignment, you will modify an existing R script
`command_line_scripts/compare_assemblies.R`(originally designed to process
two FASTA index files) so that it accepts **any number**  of FASTA files representing
genome assemblies. For each FASTA file, you will compute key assembly statistics such as:

- **Total Assembly Length:**  Sum of the lengths of all sequences.
- **N50:**  The contig (sequence) length at which the cumulative sum of lengths reaches at
  least 50% of the total assembly length.
- **L50:**  The minimal number of contigs needed to reach that 50% threshold.

You will extract sequence lengths directly from the FASTA files using the **Biostrings**
package. In addition, you will enhance the command-line interface using the **optparse**
library so that the list of files is provided as positional arguments after any
options—making it as simple as using something like `*.fasta` on the command line.
Finally, you will produce a comparative plot of the cumulative contig lengths for each
assembly.


---

### Task 1: Adjusting Command-Line Options

2. **Use Positional Arguments for the FASTA Files**

Update the script so that FASTA files are provided as positional arguments.
  Command line then looks like:

```bash
./compare_assemblies.R --outbase my_output assembly1.fasta assembly2.fasta assembly3.fasta ...
```

The optparse library in R simplifies the creation of command-line scripts by allowing you
to define and parse both options (flags starting with `--` or `-`) and positional arguments (
extra arguments listed after the options). With positional arguments, your script can
accept file names or other inputs directly from the command line, making it easier to
handle multiple values (like a list of files using wildcards such as `*.fasta`) without
needing to explicitly separate them in an option. This setup helps streamline the way
users interact with your script, providing clear help messages and robust error checking
when required inputs are missing.

Example of code using positional arguments:

```bash
# Load the optparse library
library(optparse)

# Define command-line options.
# Here we define one option: --outbase, which specifies a base name for output files.
option_list <- list(
  make_option(c("--outbase"), type = "character", default = "output",
              help = "Base name for output files [default = %default]", metavar = "character")
)

# Create an OptionParser object.
# The usage string shows that after options, positional arguments (like input files) must be provided.
opt_parser <- OptionParser(usage = "Usage: %prog [options] input1 input2 ...", option_list = option_list)

# Parse the command-line arguments.
# The result 'args' contains both 'options' (the named arguments) and 'args' (the positional arguments).
args <- parse_args(opt_parser, positional_arguments = TRUE)

# Extract the parsed options (e.g., the output base name)
opt <- args$options

# Extract the positional arguments (for example, input file names)
input_files <- args$args

# For demonstration, print the parsed options and positional arguments
cat("Output base name:", opt$outbase, "\n")
cat("Input files:", paste(input_files, collapse = ", "), "\n")

````
Inspect the code snippet above to see how the `optparse` library is used to define options and positional argument and adjust
the code of original script to accept positional arguments. 

<details>
<summary>💡 Hint</summary>
Your code snippet might look like:

```R
library(optparse)

# Define command-line options only for output basename
option_list <- list(
  make_option(c("--outbase"), type = "character", default = "assembly_output",
              help = "Base name for output files [default = %default]", metavar = "character")
)

# Update usage to reflect that the FASTA files are positional arguments
opt_parser <- OptionParser(usage = "Usage: %prog [options] fasta_file1 fasta_file2 ...", option_list = option_list)

# Parse arguments: the result will include both options and positional arguments
args <- parse_args(opt_parser, positional_arguments = TRUE)

# Extract the options and the positional arguments (FASTAs)
opt <- args$options
fasta_files <- args$args

# Ensure at least one FASTA file is provided
if (length(fasta_files) < 1) {
  print_help(opt_parser)
  stop("At least one FASTA file must be provided as a positional argument.", call. = FALSE)
}
```

</details>


---

### Task 2: Using Biostrings to Read FASTA Files

2. **Load and Utilize the Biostrings Package**

- **Step:**  Update the script to load full FASTA files and extract sequence lengths,
  rather than reading precomputed FASTA index files.
- **Instructions:**
    - Load the **Biostrings**  package.
    - For each provided FASTA file (from the positional arguments), use the function
      `readDNAStringSet()` to read the sequences.
    - Use the `width()` function to extract the lengths of all sequences in the file.

<details>
<summary>💡 Hint</summary>

```R
library(Biostrings)

# For a given FASTA file (fa_file):
seqs <- readDNAStringSet(fa_file)
contig_lengths <- width(seqs)
```
</details>

---

### Task 3: Updating the Assembly Statistics Calculation

2. **Calculate Assembly Metrics for Multiple Assemblies**

- **Step:**  Adapt your assembly statistics function to process FASTA files, extracting
  sequence lengths from each file instead of reading an index.
- **Instructions:**
    - For each FASTA file, calculate:
        - **Total Assembly Length:**  Sum of all sequence (contig) lengths.
        - **N50:**  Determine the contig length at which the cumulative sorted contig
          lengths (in descending order) reaches or exceeds 50% of the total length.
        - **L50:**  Count the minimal number of contigs needed to reach that threshold.
    - Store both the calculated statistics and the sorted lengths (for plotting purposes)
      for each assembly.

<details>
<summary>💡 Hint</summary> 

```R
calculate_assembly_stats <- function(fa_file) {
  seqs <- readDNAStringSet(fa_file)
  contig_lengths <- width(seqs)
  total_length <- sum(contig_lengths)

  sorted_lengths <- sort(contig_lengths, decreasing = TRUE)
  cumsum_lengths <- cumsum(sorted_lengths)

  L50 <- which(cumsum_lengths >= total_length / 2)[1]
  N50 <- sorted_lengths[L50]

  list(total_length = total_length, N50 = N50, L50 = L50,
       sorted_lengths = sorted_lengths, cumsum = cumsum_lengths)
}
```
</details>

- **Tip:**  Loop through the vector `fasta_files` so that each file’s statistics are
  computed.

---

### Task 4: Creating a Summary Table

2. **Assemble a Comparative Summary**

- **Step:**  Create a summary data frame that collects the assembly name (or file name),
  total assembly length, N50, and L50 for each provided FASTA file.

- **Instructions:**
    - After processing all files, combine the results into a single data frame.
    - Save this summary as a CSV file using the `write.csv()` function.

<details>
<summary>💡 Hint</summary>

```R

Your summary data frame could be structured as:

```R
summary_df <- data.frame(
  Assembly = sapply(fasta_files, basename),
  Total_Length = sapply(stats_list, function(x) x$total_length),
  N50 = sapply(stats_list, function(x) x$N50),
  L50 = sapply(stats_list, function(x) x$L50)
)
write.csv(summary_df, file = paste0(opt$outbase, "_summary.csv"), row.names = FALSE)
```
</details>
---

### Task 5: Plotting Cumulative Contig Lengths

2. **Plot Comparison Across Assemblies**

- **Step:**  Create one plot that overlays the cumulative contig length curves for each
  assembly.

- **Instructions:**

    - For each assembly, use the sorted contig lengths and cumulative sum computed previously.
    - Initialize a plot with proper x-axis (contig index) and y-axis (cumulative length)
      limits to encompass all provided assemblies.
    - Overlay the cumulative plots—using different symbols or line types—and add a legend
      identifying each assembly.

