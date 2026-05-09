################################################################################
# R SCRIPT FOR SESSION 12: PHYLOGENETIC ANALYSIS
#
# This script covers:
# 1. Loading and inspecting phylogenetic trees with ape
# 2. Tree manipulations: rooting, ladderizing, pruning, clade extraction
# 3. Base ape plotting (plot.phylo, tiplabels, nodelabels)
# 4. ggtree: ggplot2-based tree visualization and layouts
# 5. Annotating trees with metadata using %<+%
# 6. Time-calibrated trees with mrsd
# 7. Attaching heatmaps with gheatmap (pivot_wider → gheatmap workflow)
#
# DATA SOURCES USED:
# - data/phylogenetics/HIV_dentist.nhx: HIV phylogeny from the Florida dentist case
# - data/phylogenetics/HIV_env_monkeys.nhx: HIV/SIV env gene phylogeny
# - bird.orders: built-in ape dataset
# - sample.nwk: bundled with treeio package
# - MCC_FluA_H3.tree + Genotype.txt: bundled with ggtree package
#
# INSTRUCTIONS:
# Run the script section by section, following along with the slides.
# Each "SLIDE:" comment marks where to switch to the next slide.
################################################################################

################################################################################
# SECTION 1: PACKAGES
################################################################################

# ---------- SLIDE: Packages installation ----------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE))
  install.packages("ape", update = FALSE)
if (!requireNamespace("treeio", quietly = TRUE))
  BiocManager::install("treeio", update = FALSE)
if (!requireNamespace("ggtree", quietly = TRUE))
  BiocManager::install("ggtree", update = FALSE)
if (!requireNamespace("phangorn", quietly = TRUE))
  install.packages("phangorn", update = FALSE)

# ---------- SLIDE: Load packages ----------

library(ape)
library(phangorn)
library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(tibble)

################################################################################
# SECTION 2: APE BASICS
################################################################################

# ---------- SLIDE: ape package ----------

# ape provides: read.tree/write.tree, rtree, root/drop.tip/ladderize, plot.phylo

# ---------- SLIDE: Functions to read and write trees ----------

# read.tree() reads Newick format
# Newick uses parentheses: ((A:0.1,B:0.2):0.3,C:0.4);
newick_tree <- "((A:0.1,B:0.2):0.3,C:0.4);"
tree <- read.tree(text = newick_tree)
plot(tree)
axis(1)

# ---------- SLIDE: Plotting trees with ape ----------

tr <- read.tree(text = "((Pan:5,Human:5):2,Gorilla:7);")
par(mfrow = c(2, 3))
plot(tr, type = "cladogram", main = "cladogram")
plot(tr, type = "unrooted",  main = "unrooted")
plot(tr, type = "fan",       main = "fan")
plot(tr, type = "radial",    main = "radial")
plot(tr, type = "tidy",      main = "tidy")

# ---------- SLIDE: Loading tree from file ----------

tr_dentist <- read.tree("data/phylogenetics/HIV_dentist.nhx")
plot(tr_dentist)

# ---------- SLIDE: Properties of phylo class ----------

tr <- rtree(5)   # random tree with 5 tips
print(tr)
class(tr)
str(tr)
names(tr)

# ---------- SLIDE: Properties of phylo class - edge ----------

# edge: matrix where each row is one branch (parent node -> child node)
plot.phylo(tr, type = "tidy", show.tip.label = TRUE, cex = 1)
print(tr$edge)

# ---------- SLIDE: Properties of phylo class - tip.label and edge.length ----------

# tip.label: names of the leaves; edge.length: branch lengths
print(tr$tip.label)
print(tr$edge.length)
plot(tr, type = "tidy", show.tip.label = TRUE, cex = 1)
axis(1)

# ---------- SLIDE: Properties of phylo class - numbered diagram ----------

tr$tip.label  <- as.character(1:Ntip(tr))
tr$node.label <- as.character((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))
cat("Edge matrix:\n")
print(tr$edge)
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
  text  = seq_len(nrow(tr$edge)),
  frame = "none",
  adj   = 0.5,
  cex   = 1,
  col   = "red"
)

################################################################################
# SECTION 3: TREE MANIPULATION
################################################################################

# ---------- SLIDE: Manipulating trees ----------

tips <- c("Outgroup", paste0("Species", 1:5))
set.seed(123)
unrooted_tree <- rtree(6, tip.label = tips, rooted = FALSE)
is.rooted(unrooted_tree)
plot(unrooted_tree, main = "Unrooted tree")

# ---------- SLIDE: Rooting a tree ----------

rooted_tree <- root(unrooted_tree, outgroup = "Outgroup", resolve.root = TRUE)
is.rooted(rooted_tree)
plot(rooted_tree, main = "Rooted tree")

# ---------- SLIDE: Rooting a tree ----------

# reroot at a different tip
rooted_tree <- root(unrooted_tree, outgroup = "Species1", resolve.root = TRUE)
plot(rooted_tree, main = "Rooted tree at Species1")

# ---------- SLIDE: Ladderizing a tree ----------

tree2 <- rtree(8)
par(mfrow = c(1, 2))
plot(tree2,             main = "Original ordering")
plot(ladderize(tree2),  main = "After ladderize()")

# ---------- SLIDE: Pruning and dropping tips ----------

tr <- rtree(10)
par(mfrow = c(1, 2))
plot(tr, main = "Original tree")
pruned_tree <- drop.tip(tr, c("t1", "t2"))
plot(pruned_tree, main = "Pruned tree")

################################################################################
# SECTION 4: CLADE EXTRACTION
################################################################################

# ---------- SLIDE: Extracting clades ----------

# getMRCA(): most recent common ancestor node of a set of tips
# extract.clade(): subtree descended from that node
# nodepath(): path of nodes from root to a tip

# ---------- SLIDE: Understanding node numbering ----------

# Tips are numbered 1..N; internal nodes (N+1)..(N+Nnode)
set.seed(123)
tr <- rtree(3)
tr$tip.label  <- as.character(1:Ntip(tr))
tr$node.label <- as.character((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))
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
  text  = seq_len(nrow(tr$edge)),
  frame = "none",
  adj   = 0.5,
  cex   = 1,
  col   = "red"
)

# ---------- SLIDE: Getting MRCA ----------

mrca_node <- getMRCA(tr, c("1", "2"))
print(mrca_node)
plot(tr, show.tip.label = TRUE, show.node.label = TRUE,
     type = "phylogram", use.edge.length = FALSE, tip.color = c(2, 2, 1))

# ---------- SLIDE: Extracting clades descended from MRCA ----------

clade <- extract.clade(tr, mrca_node)
plot(clade, show.tip.label = TRUE, show.node.label = TRUE,
     type = "phylogram", use.edge.length = FALSE)

# ---------- SLIDE: Extracting clades descended from MRCA ----------

# Real-world example: dentist HIV tree
tr_dentist_rooted <- root(
  tr_dentist,
  outgroup     = which(tr_dentist$tip.label %in% "isolate_ELI/19-301"),
  resolve.root = TRUE
)
MRCA_node       <- getMRCA(tr_dentist_rooted, grep("Dentist", tr_dentist_rooted$tip.label))
tr_dentist_clade <- extract.clade(tr_dentist_rooted, MRCA_node)
par(mfrow = c(1, 2))
color <- grepl("Dentist", tr_dentist_rooted$tip.label) + 1
plot(tr_dentist_rooted, type = "phylogram", use.edge.length = FALSE, tip.color = color)
plot(tr_dentist_clade,  type = "phylogram", use.edge.length = FALSE)

################################################################################
# SECTION 5: BASE APE PLOTTING
################################################################################

# ---------- SLIDE: Basic plotting with plot.phylo ----------

data(bird.orders)
plot(bird.orders, cex = 0.6, main = "Bird Orders (phylogram)")

# ---------- SLIDE: Basic plotting with plot.phylo ----------

plot(bird.orders, type = "fan",      cex = 0.7, tip.color = "blue", main = "Bird Orders (fan layout)")

# ---------- SLIDE: Basic plotting with plot.phylo ----------

plot(bird.orders, type = "unrooted", cex = 0.7, tip.color = "blue", main = "Bird Orders (unrooted layout)")

# ---------- SLIDE: Adding annotations in base plots ----------

# nodelabels(), tiplabels(), edgelabels() add overlays to an existing plot
plot(bird.orders, cex = 0.6, main = "Bird Orders (phylogram)")
nodelabels(frame = "circle", cex = 0.5, col = "red")

# ---------- SLIDE: Setting colors for tip labels ----------

tip_colors <- ifelse(bird.orders$tip.label %in% c("Anseriformes", "Galliformes"), "red", "black")
plot(bird.orders, cex = 1, main = "Bird Orders (phylogram)", tip.color = tip_colors)

# ---------- SLIDE: Setting colors for tip labels ----------

# Color the entire clade defined by two taxa
mrca_node         <- getMRCA(bird.orders, c("Anseriformes", "Struthioniformes"))
clade             <- extract.clade(bird.orders, mrca_node)
tip_labels_clade  <- clade$tip.label
tip_colors        <- ifelse(bird.orders$tip.label %in% tip_labels_clade, "red", "black")
plot(bird.orders, cex = 1, main = "Bird Orders (phylogram)", tip.color = tip_colors)

# TASK 1:
# Plot the bird.orders tree as a fan layout.
# Color tips of the clade whose MRCA is shared by "Passeriformes" and "Piciformes" in blue,
# all other tips in black.
# Hint: use getMRCA(), extract.clade(), and the tip.color argument of plot().

################################################################################
# SECTION 6: GGTREE
################################################################################

# ---------- SLIDE: Visualizing trees with ggtree ----------

nwk  <- system.file("extdata", "sample.nwk", package = "treeio")
tree <- read.tree(nwk)

ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

# ---------- SLIDE: Visualizing trees with ggtree ----------

ggtree(tree, color = "firebrick", size = 2, linetype = "dotted")

# ---------- SLIDE: ggtree layouts ----------

set.seed(2017 - 02 - 16)
tree50 <- rtree(50)
ggtree(tree50)
ggtree(tree50, layout = "roundrect")
ggtree(tree50, layout = "slanted")
ggtree(tree50, layout = "ellipse")
ggtree(tree50, layout = "circular")
ggtree(tree50, layout = "fan", open.angle = 120)
ggtree(tree50, layout = "equal_angle")
ggtree(tree50, layout = "daylight")
ggtree(tree50, branch.length = "none")
ggtree(tree50, layout = "ellipse", branch.length = "none")
ggtree(tree50, branch.length = "none", layout = "circular")
ggtree(tree50, layout = "daylight", branch.length = "none")

# ---------- SLIDE: Phylogenetic Tree Annotation ----------

ggtree(tree) +
  geom_tiplab(size = 3, color = "blue") +
  geom_treescale(x = 0.5, y = 0.5, width = 0.1)

# ---------- SLIDE: Customizing tree appearance ----------

# ggtree is ggplot2-based: use theme_tree(), xlim(), geom_nodepoint(), etc.

# ---------- SLIDE: Highlighting or annotating clades ----------

ggtree(bird.orders) +
  geom_tiplab() +
  geom_treescale() +
  geom_hilight(node = 30, fill = "goldenrod", alpha = 0.3) +
  ggplot2::xlim(0, 50)

################################################################################
# SECTION 7: ANNOTATING TREES WITH METADATA
################################################################################

# ---------- SLIDE: Annotation Tree with Metadata ----------

# Common metadata types: traits, taxonomy, geography, temporal data, statistics

# ---------- SLIDE: Combining Tree with Tip Data ----------

tree_example <- rtree(5, tip.label = LETTERS[1:5])
trait_data <- data.frame(
  label  = LETTERS[1:5],
  Status = c("Endangered", "Not Endangered", "Not Endangered", "Endangered", "Endangered")
)
trait_data

# ---------- SLIDE: Attaching data to the tree ----------

# %<+% operator attaches a data frame to a ggtree plot by matching the "label" column
p <- ggtree(tree_example) %<+% trait_data +
  geom_tiplab(aes(color = Status), size = 4) +
  theme_tree() +
  scale_color_manual(values = c("red", "green"))
p

# ---------- SLIDE: Attaching data to the tree ----------

p + geom_tippoint(aes(color = Status), size = 3)

# ---------- SLIDE: Attaching data to the tree - continuous trait example ----------

tree_example2 <- rtree(20, tip.label = LETTERS[1:20])
trait_data2 <- data.frame(
  label      = LETTERS[1:20],
  TraitValue = rnorm(20, mean = 5, sd = 2)
)
p2 <- ggtree(tree_example2) %<+% trait_data2 +
  geom_tippoint(aes(color = TraitValue), size = 6)
p2

# TASK 2:
# Using the bird.orders tree, create a data frame assigning each tip to one of two groups:
# "Waterbirds" for Anseriformes, Gaviiformes, Podicipediformes, Pelecaniformes,
# Ciconiiformes, Sphenisciformes, Procellariiformes; and "Other" for the rest.
# Plot the tree with ggtree, coloring tip points by group using %<+%.
# Hint: use geom_tippoint(aes(color = Group)).

################################################################################
# SECTION 8: TIME-CALIBRATED TREES
################################################################################

# ---------- SLIDE: Joining tree with metadata using as.treedata ----------

# When metadata has many columns, join with full_join() + as.treedata().
# as_tibble(tree) produces a tibble with a "label" column for each tip.
# meta <- read_tsv("data/ncov_metadata.tsv")
# td <- full_join(as_tibble(tree), meta, by = c("label" = "strain")) |>
#         as.treedata()
# ggtree(td, mrsd = max(meta$decimal_date), as.Date = TRUE) +
#   geom_tippoint(aes(color = Nextstrain_clade))

# ---------- SLIDE: Time-calibrated trees with mrsd ----------

# mrsd = most-recent sampling date; as.Date = TRUE shows calendar x-axis;
# theme_tree2() adds the visible time axis.
beast_file <- system.file("examples/MCC_FluA_H3.tree", package = "ggtree")
beast_tree  <- read.beast(beast_file)

p_beast <- ggtree(beast_tree, mrsd = "2013-01-01") +
  geom_tiplab(size = 2, align = TRUE, linesize = 0.5) +
  theme_tree2()
p_beast

# ---------- SLIDE: Checking branch-length units ----------

# mrsd only works when branches are in years (time-calibrated tree).
# node.depth.edgelength() returns root-to-node distances;
# the maximum equals the temporal span of the tree.

# For a SARS-CoV-2 tree spanning ~6 years:
# max(node.depth.edgelength(tree))   # should be ~6, not ~0.001

################################################################################
# SECTION 9: ATTACHING HEATMAPS WITH GHEATMAP
################################################################################

# ---------- SLIDE: Attaching heatmap to the tree ----------

beast_file    <- system.file("examples/MCC_FluA_H3.tree", package = "ggtree")
beast_tree    <- read.beast(beast_file)
genotype_file <- system.file("examples/Genotype.txt", package = "ggtree")
genotype      <- read.table(genotype_file, sep = "\t", stringsAsFactors = FALSE)

p3 <- ggtree(beast_tree, mrsd = "2013-01-01") +
  geom_tiplab(size = 2, align = TRUE, linesize = 0.5) +
  theme_tree2()
head(genotype)

# ---------- SLIDE: Attaching heatmap to the tree ----------

gheatmap(p3, genotype, offset = 8, width = 0.6,
         colnames = FALSE, legend_title = "genotype") +
  scale_x_ggtree()

# ---------- SLIDE: pivot_wider for gheatmap ----------

# gheatmap requires a WIDE matrix: rows = tip names, columns = variables.
# If your data is long (one row per tip × variable), use pivot_wider first.
#
# Long format (one row per tip+mutation):
#   strain          mutation
#   Alpha/UK/001    S:N501Y
#   Alpha/UK/001    S:D614G
#   Delta/IND/002   S:L452R
#
# Convert to wide:
# mut_wide <- spike_long |>
#   filter(mutation %in% key_muts) |>
#   mutate(present = "present") |>
#   pivot_wider(names_from = mutation, values_from = present,
#               values_fill = "absent") |>
#   column_to_rownames("strain")

# ---------- SLIDE: gheatmap row-name alignment and pitfalls ----------

# CRITICAL: gheatmap aligns by ROW NAME, not row order.
# Two preparations are always needed before calling gheatmap:

# 1. Force every expected column to exist
#    (pivot_wider silently drops columns for mutations absent from all tips):
# for (m in setdiff(key_muts, colnames(mut_wide))) mut_wide[[m]] <- "absent"
# mut_wide <- mut_wide[, key_muts]

# 2. Force every tree tip to be a row
#    (tips with no mutations are dropped by the pivot):
# missing_tips <- setdiff(tree$tip.label, rownames(mut_wide))
# if (length(missing_tips) > 0) {
#   pad <- as.data.frame(
#     matrix("absent", nrow = length(missing_tips), ncol = ncol(mut_wide),
#            dimnames = list(missing_tips, colnames(mut_wide)))
#   )
#   mut_wide <- rbind(mut_wide, pad)
# }

# TASK 3:
# Using the genotype data bundled with ggtree (loaded above as `genotype`),
# inspect its structure with head() and class().
# Then call gheatmap() with offset = 5, width = 0.4, and colnames = TRUE.
# Try changing the fill colors with scale_fill_manual().
# Hint: the base plot p3 is already defined above.
