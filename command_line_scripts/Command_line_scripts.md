### Introduction to Command-Line Scripting in R on Linux 


In Linux you can write R scripts that run directly from the command line. To do this, include a 'shebang' line at the top of your R script file:

```r
#!/usr/bin/env Rscript
```

This line tells the system to use the Rscript executable (found in your PATH) to run the script. Once you have your script (for example, `calculate_N50_L50.R` or `comapare_assemblies.R`), you need to make it executable. You can do this by running in the terminal:

```bash
chmod +x calculate_N50_L50.R
chmod +x comapare_assemblies.R
```

Many scripts support a help option (e.g. using the `-h` flag) so that if you run the
script with `-h` or `--help`, it will print usage instructions. This is often implemented
in R using a package like **optparse** .

Below are examples of how to run these scripts from the command line:

Example for the First Script: `calculate_N50_L50.R`
Suppose this script is designed to take a FASTA index file (a `.fai` file) as a positional
argument and calculate assembly statistics (N50 and L50). With your data stored in the
`data/assemblies` folder and a file named `asm1.fasta.fai`, you could run in the terminal:


```bash
./calculate_N50_L50.R data/assemblies/asm1.fasta.fai output_stats.txt
```

Here, the first argument is the input FASTA index file and the second argument is the
output file where the assembly statistics will be written.

Example for the Second Script: `compare_assemblies.R`
This script uses the **optparse**  library to handle named options. It expects two input
FASTA index files for two genomes (using options like `--g1` and `--g2`), and it
calculates the assembly statistics for each genome and produces a cumulative plot.
Assuming your data files in `data/assemblies` are named as follows:

- Genome 1: `asm1.fasta.fai`
- Genome 2: `asm3.fasta.fai` (or you could choose `asm4.fasta.fai`)

You can run the script like this:

```bash
./compare_assemblies.R --g1 data/assemblies/asm1.fasta.fai --g2 data/assemblies/asm3.fasta.fai --outbase my_output
```

In this example:
 
- `--g1` specifies the FASTA index for Genome 1.
- `--g2` specifies the FASTA index for Genome 2.
- `--outbase` sets the base name for all output files (e.g., the CSV summary and plot image).

If you need help using the script, running it with the `-h` or `--help` option will print the help message (assuming the script is programmed to support that):

```bash
./comapare_assemblies.R -h
```

