## Introduction to R and RStudio 

### Why R for Biologists 
- Brief history and strengths of R in data science and bioinformatics
- Installing R and RStudio (or other IDEs)
- https://posit.co/download/rstudio-desktop/
### Navigating RStudio 
- Overview of the RStudio interface: console, editor, environment, plots, and help panes
- https://docs.posit.co/ide/user/ide/guide/ui/ui-panes.html
- https://www.youtube.com/watch?v=T44GWKlFFRI&t=6s
- Customizing the workspace (themes, layouts, etc.)
- Setting the working directory

### Language Fundamentals 
- https://cran.r-project.org/doc/manuals/R-intro.html#Getting-help
- Basic syntax (expressions, statements, operators)
- Variables: naming conventions, assignment operators (`<-`, `=`)
- Basic mathematical operations

### Data Structures Overview 

- Atomic vectors (numeric, character, logical, etc.) (https://cran.r-project.org/doc/manuals/R-intro.html#Simple-manipulations-numbers-and-vectors)
- Factors and their role in categorical data
- Matrices and arrays (https://cran.r-project.org/doc/manuals/R-intro.html#Arrays-and-matrices)
- Lists and nested lists (https://cran.r-project.org/doc/manuals/R-intro.html#Lists-and-data-frames)
- Data frames (https://cran.r-project.org/doc/manuals/R-intro.html#Data-frames)

### Finding Help and Documentation 
 
- Using `?` and `help()` (https://cran.r-project.org/doc/manuals/R-intro.html#Getting-help)
- Vignettes and package manuals (example: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- Online resources (CRAN, Bioconductor, RStudio Community)
  - https://cran.r-project.org/
  - https://www.bioconductor.org/
  - https://community.rstudio.com/

### Working with Packages 
 
- Installing packages (`install.packages`, BiocManager for Bioconductor)
 
- Loading and updating packages (`library`, `require`, `update.packages`)

- Package directories and version management

---

## Importing and Manipulating Tabular Data 

### Reading Data 
 
- Reading CSV, TSV, and other delimited files (`read.csv`, `read.delim`, `readr` package)
 
- Reading Excel files (`readxl` package)

- Reading data from the clipboard and online resources

### Data Frames and Tibbles 

- Converting between data frames and tibbles
 
- Inspecting data (head, tail, `str`, `summary`)

### Data Export 

- Writing CSV, TSV, and Excel files

- Exporting to other formats (TXT, RDS, RData)

### Handling Missing Data 
 
- Identifying missing values (`NA`, `NaN`)

- Strategies for imputation or filtering


---


## Control Structures and Functions 

### Control Structures 
 
- Conditional statements: `if`, `else if`, `else`
 
- Loops: `for`, `while`, `repeat`

- Best practices for readable and efficient looping

### Vectorized Computation 

- Element-wise operations on vectors

- Recycling rules for vectors

- Comparison with loops for performance

### apply Family of Functions 
 
- `apply`, `lapply`, `sapply`, `tapply`, `mapply`

- When to use each function

- Advantages in code simplicity and speed

### Writing Your Own Functions 

- Function definition syntax

- Default arguments and argument matching

- Return values

- Documentation with Roxygen2 (intro)

### Debugging 

- Common errors and warnings
 
- Using `traceback()`, `browser()`, `debug()`, and `recover()`

- Strategies for troubleshooting


---


## Base Graphics and Plotting 

### Basic Plots 

- Scatter plots, line plots, bar charts, histograms, boxplots
 
- Using `plot()`, `barplot()`, `hist()`, etc.

### Customization 

- Titles, axis labels, legends, colors, point shapes, line types

- Annotation and text placement
 
- Combining multiple plots on one canvas (`par`, `layout`)

### Exporting Plots 

- Saving plots to files (PDF, PNG, TIFF)

- Adjusting resolution and dimensions


---


## ggplot2 and the Tidyverse Ecosystem 

### Tidy Data Principles 

- What is “tidy” data?

- Long vs. wide formats

### ggplot2 Grammar of Graphics 
 
- Structure of a `ggplot` command (data, aesthetics, geometries)

- Layers, scales, themes
 
- Common geoms: `geom_point`, `geom_line`, `geom_bar`, `geom_histogram`, `geom_boxplot`, etc.

### Faceting and Custom Themes 
 
- `facet_wrap` and `facet_grid` for multi-panel plots
 
- Built-in themes (`theme_bw`, `theme_minimal`)

- Creating and sharing custom themes

### Extensions and Additional Plot Types 
 
- `ggrepel` for improved labeling
 
- `ggforce` for advanced geometries

- Other specialized geoms (violin plots, ridgeline plots, etc.)


---


## Data Cleaning and Manipulation with Tidyverse 

### Key tidyr Functions 
 
- `gather`/`pivot_longer` and `spread`/`pivot_wider`
 
- `separate` and `unite`

- Handling missing data

### Key dplyr Verbs 
 
- `select`, `filter`, `mutate`, `arrange`, `summarize`, `group_by`
 
- Piping with `%>%`
 
- Joins (`inner_join`, `left_join`, `right_join`, `full_join`)

### stringr Essentials 

- Pattern matching with regex

- Searching and replacing text

- String manipulation (splitting, substring extraction, etc.)

### Real-World Data Wrangling Examples 

- Combining multiple data sources

- Reshaping datasets for analysis


---


## RMarkdown for Reporting 

### Document Structure 

- YAML header basics (title, author, date)

- Code chunks and chunk options

- Inline R code, text formatting

### Exporting Reports 

- HTML, PDF, Word

- Parameterized reports

- Setting up bibliographies and citations

### Writing Reports with knitr 

- Controlling output with chunk options

- Caching and code efficiency

- Embedding plots, tables, and images

### Command Line Execution 

- Rendering RMarkdown from the terminal

- Automated workflows (Makefiles, scripts)


---


## Introduction to Bioconductor 

### What is Bioconductor? 
 
- Bioconductor vs. CRAN, installation with `BiocManager`

- Overview of core Bioconductor packages

### The Biostrings Package 

- Working with DNA, RNA, and protein sequences

- Basic sequence manipulation (subsetting, complementing, reversing)

- Pattern searching (regular expressions, fuzzy matching)

### Common Use Cases 

- Calculating GC content, sequence translation, etc.

- Alignments (brief introduction)


---


## Working with Genomic Intervals using GenomicRanges 

### Introduction to GRanges Objects 
 
- Creating Genomic Ranges (`GRanges`, `IRanges`)

- Storing metadata (scores, sample info)

### Interval Operations 

- Overlaps, subsetting by range, finding nearest features

- Merging and resizing ranges

### Annotations 

- Annotating genomic features (promoters, exons, transcripts)
 
- Using `AnnotationHub` and `TxDb` packages


---


## Visualizing Genomic Data 

### ggplot2 for Genomic Data 

- Using genomic ranges in ggplot2

- Plotting coverage, read depth, and other region-based metrics

### Introduction to ggbio 

- Creating genome browser-like tracks

- Customizing track layouts (ideograms, coverage, gene models)

### Creating Density Plots and Heatmaps 

- Heatmaps of read coverage or expression data
 
- Using `pheatmap` or `ComplexHeatmap`

### Integration with Other Tools 

- Combining ggbio/ggplot2 outputs

- Annotating with external data sources


---


## Gene Expression Analysis with RNA-Seq 

### RNA-Seq Workflow Overview 

- Experimental design, read alignment, and counting

- Raw count data vs. normalized data

### edgeR and limma 

- Installing and loading packages

- Data input, filtering, and normalization (CPM, TMM)

- Model design and differential expression testing

### Exploratory Data Analysis 

- Principal Component Analysis (PCA)

- Clustering and heatmaps

- Data transformations (log2, variance-stabilizing)

### Hypothesis Testing and Multiple Testing Correction 

- Statistical concepts (p-values, FDR)

- Visualizing results with MA plots, volcano plots

### Results Interpretation 

- Functional annotation (brief intro to GO, KEGG)

- Reporting findings in tables and plots


---


## Phylogenetic Analysis and Visualization 

### Overview of Phylogenetic Packages 
 
- Installing `ape`, `treeio`, and related packages

- Basic concepts of phylogenetics (trees, distance matrices)

### Building Phylogenetic Trees 

- Importing and exporting tree files (Newick, Nexus)

- Distance methods vs. character-based methods (brief overview)
 
- Tree construction functions in `ape`

### Tree Visualization and Annotation 
 
- Plotting phylogenetic trees with `plot.phylo`

- Customizing labels, colors, and styles
 
- Using `ggtree` from the `treeio` package for advanced visualization

### Molecular Evolution Analyses (Optional/Extended Topic) 

- Estimating divergence times

- Checking models of evolution

- Basic comparative methods (e.g., phylogenetic independent contrasts)


---


## Additional Topics / Capstone Projects (Optional) 

- Integrating R with command-line tools (system commands, shell scripts)

- Accessing databases and APIs (NCBI, UniProt)

- Interactive dashboards with Shiny (for advanced users)

- Batch processing and reproducible pipelines (Snakemake, Nextflow)

- Group project ideas or case studies (e.g., analyzing a public RNA-Seq dataset from GEO)
