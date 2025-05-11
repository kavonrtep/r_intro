# ---
# title: "Using phylogenetic trees in R"
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

#
#
# # Phylogenetic Analysis in R with `ape` and `treeio`
#
# Phylogenetic trees are key to understanding evolutionary relationships. In R, the ape package (Analysis of Phylogenetics and Evolution) and the Bioconductor treeio package provide powerful tools to import, manipulate, and visualize phylogenetic trees.
#
# These packages can be used to:
#
# - Load phylogenetic trees into R (from formats like Newick/Nexus) using ape (and treeio).
# - Create simple example trees
# - Perform common tree manipulations: rooting/unrooting, ladderizing, pruning (dropping tips), and basic traversal of tree structure.
# - Visualize trees using base ape plotting (via plot.phylo) and enhanced graphics with ggtree (built on treeio).
# - Annotate trees with metadata (e.g. traits, taxonomy, geography) and visualize these data on the tree (e.g. coloring tips by trait).

# ## Packages installation:

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ape", quietly=TRUE))
  install.packages("ape", update=FALSE)
if (!requireNamespace("treeio", quietly=TRUE))
  BiocManager::install("treeio", update=FALSE)
if (!requireNamespace("ggtree", quietly=TRUE))
  BiocManager::install("ggtree", update=FALSE)
if (!requireNamespace("phangorn", quietly=TRUE))
  install.packages("phangorn", update=FALSE)

# ## Load packages

library(ape)
library(phangorn)
library(treeio)
library(ggtree)

# ## `ape` package

# `ape` provides functions to read and write tree files, to store trees in R as objects of class `phylo`, and to manipulate and analyze these trees. With `ape`, you can perform tasks like computing distances, ancestral state reconstruction, and comparative analysis, but in this lesson we focus on tree handling and basic plotting. Ape’s base plotting functions (like `plot.phylo`) produce static phylograms or cladograms using base R graphics.

# ## Functions to read and write trees

# `read.tree()`: Reads a tree from a file in Newick
# - Newick format uses parentheses to represent the tree structure, with branch lengths and node labels.
# - Example: `((A:0.1,B:0.2):0.3,C:0.4);`

newick_tree <- "((A:0.1,B:0.2):0.3,C:0.4);"
tree <- read.tree(text = newick_tree)
plot(tree)
axis(1)

# ## Plotting trees with `ape`

# - `ape` provides a function `plot.phylo()` - it can be invoked by `plot()` function on `phylo` objects.

tr <- read.tree(text = "((Pan:5,Human:5):2,Gorilla:7);")
par(mfrow=c(2,3))
plot(tr, type = "cladogram", main = "cladogram")
plot(tr, type = "unrooted", main="unrooted")
plot(tr, type = "fan", main = "fan")
plot(tr, type = "radial", main = "radial")
plot(tr, type = "tidy", main = "tidy")

# ## Loading tree from file

tr_dentist <- read.tree("./data/phylogenetics/HIV_dentist.nhx")
plot(tr_dentist)

# ## Properties of `phylo` class

tr <- rtree(5) # random tree with 5 tips
print(tr)
class(tr)
str(tr)
names(tr)

# ## Properties of `phylo` class

# - `edge`: a matrix with two columns, each row represents an edge in the tree.

plot.phylo(tr, type = "tidy", show.tip.label = TRUE, cex = 0.5)
# add node labels to identify internal nodes
print(tr$edge)

# ## Properties of `phylo` class

# - `tip.label`: a character vector with the labels of the tips (leaves) of the tree.
# - `edge.length`: a numeric vector with the lengths of the edges in the tree.

print(tr$tip.label)
print(tr$edge.length)
plot(tr, type ='tidy', show.tip.label = TRUE, cex = 0.5)
axis(1)

# ## Properties of `phylo` class

tr$tip.label  <- as.character(1:Ntip(tr))
tr$node.label <- as.character((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)))
cat("Edge matrix:\n")
print(tr$edge)
plot(
  tr,
  type            = "phylogram",     # style of tree
  use.edge.length = FALSE,           # ignore branch lengths here
  show.tip.label  = TRUE,            # draw tr$tip.label
  show.node.label = TRUE,            # draw tr$node.label
  label.offset    = 0.05,            # space between tips/nodes and text
  cex             = 1,               # text size
  edge.color      = "darkgray",      # branch color
  main            = "Tree with Numbered Tips, Nodes, and Edges"
)
edgelabels(
  text = seq_len(nrow(tr$edge)),   # 1,2,…,#edges
  frame = "none",                  # no box
  adj   = 0.5,                     # center on branch midpoint
  cex   = 1,
  col   = "red"
)

# ## Manipulating trees

# - `is.rooted()`: Check if the tree is rooted.
# - `unroot()`: Unroot a tree.
# - `root()`: Root a tree at a specified node or tip.

tips <- c("Outgroup", paste0("Species", 1:5))
set.seed(123)
unrooted_tree <- rtree(6, tip.label = tips, rooted = FALSE)
is.rooted(unrooted_tree)
plot(unrooted_tree, main = "Unrooted tree")

# ## Rooting a tree

outgroup <- which(unrooted_tree$tip.label == "Outgroup")
rooted_tree <- root(unrooted_tree, outgroup = "Outgroup", resolve.root = TRUE)
is.rooted(rooted_tree)
plot(rooted_tree, main = "Rooted tree")

# ## Rooting a tree

# reroot the tree at a different node
rooted_tree <- root(unrooted_tree, outgroup = "Species1", resolve.root = TRUE)
plot(rooted_tree, main = "Rooted tree at Species1")

# ## Ladderizing a tree

# - Ladderizing a tree means ordering the branches so that the tree has a tidy, ladder-like appearance (one side of each split has the smaller clade).
# - This doesn’t change the evolutionary relationships at all—it only affects the ordering of tips.

tree2 <- rtree(8)
par(mfrow=c(1,2))
plot(tree2, main = "Original ordering")
plot(ladderize(tree2), main = "After ladderize()")

# ## Pruning and dropping tips

# - Often, you might want to prune a tree to remove certain tips (species) – for example, if you want to focus on a subset or if some tips have no data for your analysis.
# - You can use the `drop.tip()` function to remove tips from a tree.

tr <- rtree(10)
par(mfrow=c(1,2))
plot(tr, main = "Original tree")
# Drop tips "t1" and "t2"
pruned_tree <- drop.tip(tr, c("t1", "t2"))
plot(pruned_tree, main = "Pruned tree")

# ## Extracting clades

# Sometimes you want to isolate a subtree (clade) or inspect specific parts of a tree. Ape has tools to help with that:
# - `getMRCA(tree, tips=c(...))`: finds the most recent common ancestor (MRCA) node of a given set of tips.
# - `extract.clade(tree, node)`: extracts the subtree/clade descended from a given internal node.
# - `nodepath(tree, from=NULL, to=NULL)`: returns the sequence of node indices along the path between two nodes.

# ## Understanding node numbering

# - Tip nodes are numbered 1 to N (where N is number of tips)
# - Internal nodes are numbered (N+1) to (N + Nnode) where Nnode is number of internal nodes.
# - The `tree$edge` matrix has two columns: `parent` and `child` for each edge.
# - Internal nodes are typically referred to by number.

set.seed(123)
tr <- rtree(3)
tr$tip.label  <- as.character(1:Ntip(tr))
tr$node.label <- as.character((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)))
plot(
  tr,
  type            = "phylogram",
  use.edge.length = FALSE,
  show.tip.label  = TRUE,
  show.node.label = TRUE,
  label.offset    = 0.05,
  cex             = 1,
  edge.color      = "darkgray",
  main            = "Tree with Numbered Tips, Nodes, and Edges"
)
edgelabels(
  text = seq_len(nrow(tr$edge)),
  frame = "none",
  adj   = 0.5,
  cex   = 1,
  col   = "red"
)

# ## Getting MRCA

mrca_node <- getMRCA(tr, c("1", "2"))
print(mrca_node)
plot(tr, show.tip.label = TRUE, show.node.label = TRUE,
     type = "phylogram", use.edge.length = FALSE, tip.color = c(2,2,1))

# ## Extracting clades descended from MRCA

clade <- extract.clade(tr, mrca_node)
plot(clade, show.tip.label = TRUE, show.node.label = TRUE,
     type = "phylogram", use.edge.length = FALSE)

# ## Extracting clades descended from MRCA

tr_dentist_rooted <- root(tr_dentist, outgroup = which(tr_dentist$tip.label %in% "isolate_ELI/19-301"), resolve.root = TRUE)
MRCA_node <- getMRCA(tr_dentist_rooted, grep("Dentist", tr_dentist_rooted$tip.label))
tr_dentist_clade <- extract.clade(tr_dentist_rooted, MRCA_node)
par(mfrow=c(1,2))
color <- grepl("Dentist", tr_dentist_rooted$tip.label) + 1
plot(tr_dentist_rooted, type = "phylogram", use.edge.length = FALSE, tip.color = color)
plot(tr_dentist_clade, type = "phylogram", use.edge.length = FALSE)

# ## Basic plotting with `plot.phylo`

# The simplest way to plot a tree in R is: `plot(my_tree)`

data(bird.orders)
plot(bird.orders, cex=0.6, main="Bird Orders (phylogram)")

# `plot.phylo` variations

plot(bird.orders, type="fan", cex=0.7, tip.color="blue", main="Bird Orders (fan layout)")
plot(bird.orders, type="unrooted", cex=0.7, tip.color="blue", main="Bird Orders (unrooted layout)")

# ## Adding annotations in base plots

# Ape provides functions to add to an existing tree plot:
# - `tiplabels()`: add symbols/text at tips
# - `nodelabels()`: add labels at internal nodes
# - `edgelabels()`: label branches

plot(bird.orders, cex=0.6, main="Bird Orders (phylogram)")
nodelabels(frame = "circle", cex = 0.5, col = "red")

# ## Setting colors for tip labels

tip_colors <- ifelse(bird.orders$tip.label %in% c("Anseriformes", "Galliformes"), "red", "black")
plot(bird.orders, cex=1, main="Bird Orders (phylogram)", tip.color=tip_colors)

# ## Highlighting a clade by MRCA

mrca_node <- getMRCA(bird.orders, c("Anseriformes", "Struthioniformes"))
clade <- extract.clade(bird.orders, mrca_node)
tip_labels_for_clade <- clade$tip.label
tip_colors <- ifelse(bird.orders$tip.label %in% tip_labels_for_clade, "red", "black")
plot(bird.orders, cex=1, main="Bird Orders (phylogram)", tip.color=tip_colors)

# ## Visualizing trees with `ggtree`

nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

ggtree(tree, color="firebrick", size=2, linetype="dotted")

# ## ggtree layouts

set.seed(2017-02-16)
tree50 <- rtree(50)
ggtree(tree50)
ggtree(tree50, layout="roundrect")
ggtree(tree50, layout="slanted")
ggtree(tree50, layout="ellipse")
ggtree(tree50, layout="circular")
ggtree(tree50, layout="fan", open.angle=120)
ggtree(tree50, layout="equal_angle")
ggtree(tree50, layout="daylight")
ggtree(tree50, branch.length='none')
ggtree(tree50, layout="ellipse", branch.length="none")
ggtree(tree50, branch.length='none', layout='circular')
ggtree(tree50, layout="daylight", branch.length = 'none')

# ## Phylogenetic Tree Annotation

# - `geom_tiplab()`: add tip labels
# - `geom_nodelab()`: add node labels
# - `geom_treescale()`: add scale bar

ggtree(tree) +
  geom_tiplab(size=3, color="blue") +
  geom_treescale(x=0.5, y=0.5, width=0.1)

# ## Customizing tree appearance

# Use ggplot2 themes and layers (`theme_tree()`, `xlim()`, `ylim()`, `geom_point2()`, `geom_text2()`).

# ## Highlighting or annotating clades

# e.g., highlight node 30:
# ggtree(bird.orders) +
#   geom_tiplab() + geom_treescale() +
#   geom_hilight(node=30, fill="goldenrod", alpha=0.3) + xlim(0, 50)

# ## Annotation Tree with Metadata

# - Traits (body size, habitat)
# - Taxonomic info
# - Geography or host information
# - Temporal data
# - Statistical results mapped onto branches or nodes

# ## Combining Tree with Tip Data

tree_example <- rtree(5, tip.label = LETTERS[1:5])
trait_data <- data.frame(
  label = LETTERS[1:5],
  Status = c("Endangered", "Not Endangered", "Not Endangered", "Endangered", "Endangered")
)
trait_data

# ## Attaching data to the tree

p <- ggtree(tree_example) %<+% trait_data +
  geom_tiplab(aes(color = Status), size = 4) +
  theme_tree() +
  scale_color_manual(values = c("red", "green"))
p

p + geom_tippoint(aes(color = Status), size = 3)

# ## Attaching data to the tree - continuous trait example

tree_example2 <- rtree(20, tip.label = LETTERS[1:20])
trait_data2 <- data.frame(
  label = LETTERS[1:20],
  TraitValue = rnorm(20, mean = 5, sd = 2)
)
p2 <- ggtree(tree_example2) %<+% trait_data2 +
  geom_tippoint(aes(color = TraitValue), size = 6)
p2

# ## Attaching heatmap to the tree

beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactors=FALSE)
p3 <- ggtree(beast_tree, mrsd="2013-01-01") +
    geom_tiplab(size=2, align=TRUE, linesize=.5) +
    theme_tree2()
p3

head(genotype)
gheatmap(p3, genotype, offset=8, width=0.6,
        colnames=FALSE, legend_title="genotype") +
    scale_x_ggtree()
