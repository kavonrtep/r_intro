## **Introduction to R and RStudio** 1.1 **Why R for Biologists** 
- Brief history and strengths of R in data science and bioinformatics
- Installing R and RStudio (or other IDEs)

### **Navigating RStudio** 
- Overview of the RStudio interface: console, editor, environment, plots, and help panes
- Customizing the workspace (themes, layouts, etc.)
- Setting the working directory

### Language Fundamentals 
- Basic syntax (expressions, statements, operators)
- Variables: naming conventions, assignment operators (`<-`, `=`)
- Basic mathematical operations

### Data Structures Overview 
- Atomic vectors (numeric, character, logical, etc.)
- Lists and nested lists
- Factors and their role in categorical data
- Matrices and arrays
- Data frames vs. tibbles: similarities and differences

### Finding Help and Documentation  
- Using `?` and `help()`
- Vignettes and package manuals
- Online resources (CRAN, Bioconductor, RStudio Community)

### Working with Packages  
- Installing packages (`install.packages`, BiocManager for Bioconductor)
- Loading and updating packages (`library`, `require`, `update.packages`)
- Package directories and version management


---

2. **Importing and Manipulating Tabular Data** 2.1 **Reading Data**  
- Reading CSV, TSV, and other delimited files (`read.csv`, `read.delim`, `readr` package)
 
- Reading Excel files (`readxl` package)

- Reading data from the clipboard and online resources
2.2 **Data Frames and Tibbles** 
- Converting between data frames and tibbles
 
- Inspecting data (head, tail, `str`, `summary`)
2.3 **Data Export** 
- Writing CSV, TSV, and Excel files

- Exporting to other formats (TXT, RDS, RData)
2.4 **Handling Missing Data**  
- Identifying missing values (`NA`, `NaN`)

- Strategies for imputation or filtering


---

3. **Control Structures and Functions** 3.1 **Control Structures**  
- Conditional statements: `if`, `else if`, `else`
 
- Loops: `for`, `while`, `repeat`

- Best practices for readable and efficient looping
3.2 **Vectorized Computation** 
- Element-wise operations on vectors

- Recycling rules for vectors

- Comparison with loops for performance
3.3 `apply` Family of Functions**  
- `apply`, `lapply`, `sapply`, `tapply`, `mapply`

- When to use each function

- Advantages in code simplicity and speed
3.4 **Writing Your Own Functions** 
- Function definition syntax

- Default arguments and argument matching

- Return values

- Documentation with Roxygen2 (intro)
3.5 **Debugging** 
- Common errors and warnings
 
- Using `traceback()`, `browser()`, `debug()`, and `recover()`

- Strategies for troubleshooting


---

4. **Base Graphics and Plotting** 4.1 **Basic Plots** 
- Scatter plots, line plots, bar charts, histograms, boxplots
 
- Using `plot()`, `barplot()`, `hist()`, etc.
4.2 **Customization** 
- Titles, axis labels, legends, colors, point shapes, line types

- Annotation and text placement
 
- Combining multiple plots on one canvas (`par`, `layout`)
4.3 **Exporting Plots** 
- Saving plots to files (PDF, PNG, TIFF)

- Adjusting resolution and dimensions


---

5. **ggplot2 and the Tidyverse Ecosystem** 5.1 **Tidy Data Principles** 
- What is “tidy” data?

- Long vs. wide formats
5.2 **ggplot2 Grammar of Graphics**  
- Structure of a `ggplot` command (data, aesthetics, geometries)

- Layers, scales, themes
 
- Common geoms: `geom_point`, `geom_line`, `geom_bar`, `geom_histogram`, `geom_boxplot`, etc.
5.3 **Faceting and Custom Themes**  
- `facet_wrap` and `facet_grid` for multi-panel plots
 
- Built-in themes (`theme_bw`, `theme_minimal`)

- Creating and sharing custom themes
5.4 **Extensions and Additional Plot Types**  
- `ggrepel` for improved labeling
 
- `ggforce` for advanced geometries

- Other specialized geoms (violin plots, ridgeline plots, etc.)


---

6. **Data Cleaning and Manipulation with Tidyverse** 6.1 Key `tidyr` Functions**  
- `gather`/`pivot_longer` and `spread`/`pivot_wider`
 
- `separate` and `unite`

- Handling missing data
6.2 Key `dplyr` Verbs**  
- `select`, `filter`, `mutate`, `arrange`, `summarize`, `group_by`
 
- Piping with `%>%`
 
- Joins (`inner_join`, `left_join`, `right_join`, `full_join`)
6.3 `stringr` Essentials** 
- Pattern matching with regex

- Searching and replacing text

- String manipulation (splitting, substring extraction, etc.)
6.4 **Real-World Data Wrangling Examples** 
- Combining multiple data sources

- Reshaping datasets for analysis


---

7. **RMarkdown for Reporting** 7.1 **Document Structure** 
- YAML header basics (title, author, date)

- Code chunks and chunk options

- Inline R code, text formatting
7.2 **Exporting Reports** 
- HTML, PDF, Word

- Parameterized reports

- Setting up bibliographies and citations
7.3 **Writing Reports with knitr** 
- Controlling output with chunk options

- Caching and code efficiency

- Embedding plots, tables, and images
7.4 **Command Line Execution** 
- Rendering RMarkdown from the terminal

- Automated workflows (Makefiles, scripts)


---

8. **Introduction to Bioconductor** 8.1 **What is Bioconductor?**  
- Bioconductor vs. CRAN, installation with `BiocManager`

- Overview of core Bioconductor packages
8.2 **The Biostrings Package** 
- Working with DNA, RNA, and protein sequences

- Basic sequence manipulation (subsetting, complementing, reversing)

- Pattern searching (regular expressions, fuzzy matching)
8.3 **Common Use Cases** 
- Calculating GC content, sequence translation, etc.

- Alignments (brief introduction)


---

9. **Working with Genomic Intervals using GenomicRanges** 9.1 **Introduction to GRanges Objects**  
- Creating Genomic Ranges (`GRanges`, `IRanges`)

- Storing metadata (scores, sample info)
9.2 **Interval Operations** 
- Overlaps, subsetting by range, finding nearest features

- Merging and resizing ranges
9.3 **Annotations** 
- Annotating genomic features (promoters, exons, transcripts)
 
- Using `AnnotationHub` and `TxDb` packages


---

10. **Visualizing Genomic Data** 10.1 **ggplot2 for Genomic Data** 
- Using genomic ranges in ggplot2

- Plotting coverage, read depth, and other region-based metrics
10.2 **Introduction to ggbio** 
- Creating genome browser-like tracks

- Customizing track layouts (ideograms, coverage, gene models)
10.3 **Creating Density Plots and Heatmaps** 
- Heatmaps of read coverage or expression data
 
- Using `pheatmap` or `ComplexHeatmap`
10.4 **Integration with Other Tools** 
- Combining ggbio/ggplot2 outputs

- Annotating with external data sources


---

11. **Gene Expression Analysis with RNA-Seq** 11.1 **RNA-Seq Workflow Overview** 
- Experimental design, read alignment, and counting

- Raw count data vs. normalized data
11.2 **edgeR and limma** 
- Installing and loading packages

- Data input, filtering, and normalization (CPM, TMM)

- Model design and differential expression testing
11.3 **Exploratory Data Analysis** 
- Principal Component Analysis (PCA)

- Clustering and heatmaps

- Data transformations (log2, variance-stabilizing)
11.4 **Hypothesis Testing and Multiple Testing Correction** 
- Statistical concepts (p-values, FDR)

- Visualizing results with MA plots, volcano plots
11.5 **Results Interpretation** 
- Functional annotation (brief intro to GO, KEGG)

- Reporting findings in tables and plots


---

12. **Phylogenetic Analysis and Visualization** 12.1 **Overview of Phylogenetic Packages**  
- Installing `ape`, `treeio`, and related packages

- Basic concepts of phylogenetics (trees, distance matrices)
12.2 **Building Phylogenetic Trees** 
- Importing and exporting tree files (Newick, Nexus)

- Distance methods vs. character-based methods (brief overview)
 
- Tree construction functions in `ape`
12.3 **Tree Visualization and Annotation**  
- Plotting phylogenetic trees with `plot.phylo`

- Customizing labels, colors, and styles
 
- Using `ggtree` from the `treeio` package for advanced visualization
12.4 **Molecular Evolution Analyses**  (Optional/Extended Topic)
- Estimating divergence times

- Checking models of evolution

- Basic comparative methods (e.g., phylogenetic independent contrasts)


---

13. **Additional Topics / Capstone Projects**  (Optional)
- Integrating R with command-line tools (system commands, shell scripts)

- Accessing databases and APIs (NCBI, UniProt)

- Interactive dashboards with Shiny (for advanced users)

- Batch processing and reproducible pipelines (Snakemake, Nextflow)

- Group project ideas or case studies (e.g., analyzing a public RNA-Seq dataset from GEO)


---



