## `DNase.tsv`

Tab delimited file 176 rows and 3 columns of data obtained
during development of an ELISA assay for the recombinant protein
DNase in rat serum.

columns:

- `Run` - a factor indicating the assay run number
- `conc` - a numeric vector of concentrations
- `density` - a numeric vector of optical density readings

## ChickWeight.csv

The ChickWeight dataset contains 578 rows and 4 columns of data from an experiment on the
effect of diet on early growth of chicks. The columns are:

- `weight` - a numeric vector of body weights of the chicks
- `Time` - a numeric vector of the time points when the measurements were taken
- `Chick` - a factor indicating the chick ID
- `Diet` - a factor indicating the diet type

## ChickWeight_feed.txt, ChickWeight_feed.xlsx

An experiment was conducted to measure and compare the
effectiveness of various feed supplements on the growth rate of
chickens. File with txt suffix is space delimited, xlsx is an excel file. Columns include:

- `weight` - a numeric vector of body weights of the chicks
- `feed` - a factor indicating the feed type


## genes.bed

A BED file containing genomic coordinates of genes in the C. australis genome. 
BED is tab delimited file, format specification - https://www.ensembl.org/info/website/upload/bed.html
Columns include:
1. `chrom` - chromosome name
2. `start` - start position of the feature
3. `end` - end position of the feature
4. `name` - gene name (optional)
5. `score` - score (optional)
6. `strand` - strand (optional)


## feature_counts.csv

A tab-delimited file containing gene expression counts from an RNA-seq experiment.
Columns names are sample names and row names are gene names. Values are integer counts.
This dataset is accompanied by a feature annotation file `feature_annotation.csv`. Annotation file contains gene annotations including a `gender` column.

```txt
        num.tech.reps population      study gender
NA06985             1        CEU Montgomery female
NA06986             1        CEU Montgomery   male
NA06994             1        CEU Montgomery   male
NA07000             1        CEU Montgomery female
NA07037             1        CEU Montgomery female
NA07051             1        CEU Montgomery   male
```
