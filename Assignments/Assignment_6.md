# Phylogenetic Analysis of Human HIV and SIV in Monkeys phylogenetics tree

Human Immunodeficiency Viruses (HIV-1 and HIV-2) are closely related to Simian Immunodeficiency Viruses (SIVs) found in various African primate species. Understanding the evolutionary relationships between HIV and SIV strains can reveal the origins of HIV and its transmission pathways from primates to humans. Human HIV include HIV-1 and HIV-2 variant. Your task will be to visualize HIV sequences from human and SIV from monkeys and extract clade with HIV-1 and HIV-2 sequences and identify corresponding monkeys SIVs


In this assignment you will:

 
- **Load**  and **root**  an NHX phylogeny of HIV-env sequences from monkeys.
- **Identify**  and **color**  the HIV-1 and HIV-2 clades on the full tree.
- **Extract**  each clade into its own subtree.
- **Parse**  out and **list**  the monkey host taxa from each clade.

All work will use the **ape**  package and base R.  The tree file is now at `data/phylogenetics/HIV_env_monkeys.nhx`.




### Task 1: Load, Root, and Inspect the Tree 

 
**Install and load**  **ape**  if needed:

```r
install.packages("ape")
library(ape)
```
 
**Read**  the NHX tree and **root**  it on the tip
`"SIV-MON;_Mona_monkey;_AY340701"` (the outgroup)

 
**Inspect**  the rooted tree:


### Task 2: Identify HIV-1 and HIV-2 Tip Sets 
How many HIV-1 and HIV-2 tips are there?


### Task 3: Plot and Color Entire Clades 

**Plot**  the rooted tree, coloring **tip labels**  and **branches**  by clade:
 
**Hint:** 
 
  - `edge.color` takes a vector as long as the number of edges; here we assign each branch the color of its child tip or node.
 

### Task 4: Extract Subtrees for HIV-1 and HIV-2 
**Find**  the MRCA node for each clade:

 
**Extract**  each subtree

**Plot**  the two subtrees side by side:


### Task 5: List SIV tips for each clade


