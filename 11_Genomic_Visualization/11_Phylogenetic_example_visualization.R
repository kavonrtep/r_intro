



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("biomformat", quietly = TRUE))
    BiocManager::install("bioformat")
if (!requireNamespace("phyloseq", quietly = TRUE))
    install.packages("ade4")
    BiocManager::install("phyloseq")


library(biomformat)  # for reading BIOM files
library(ape)         # for reading phylogenetic trees
library(tidyverse)

metadata <- read.table("data/qiime/sample-metadata.tsv", sep="\t", header=TRUE)
# Remove leading '#' from the first column name if present

# Set row names to sample IDs (first column)
rownames(metadata) <- metadata[[1]]
metadata <- metadata[ , -1]  # drop the explicit ID column now in rownames

phy_tree <- read.tree("data/qiime/tree.nwk")

biom_obj <- read_biom("data/qiime/feature-table.biom")
otu_mat <- as.matrix(biom_data(biom_obj))

tax_tab <- read.table("data/qiime/taxonomy.tsv", sep="\t", header=TRUE)

tax_df <- tax_tab %>%
  rename(FeatureID = `Feature.ID`) %>%
  separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";")

otu_counts <-  rowSums(otu_mat, na.rm = TRUE)

rel_abun <- otu_counts / sum(otu_counts) * 100

# 8. Build an OTU metadata table
otu_metadata <- tax_df %>%
  mutate(TotalCounts = otu_counts[FeatureID],
         RelAbundance = rel_abun[FeatureID]) %>%
  select(FeatureID, Kingdom:Species, TotalCounts, RelAbundance)


p <- ggtree(phy_tree) %<+% otu_metadata +                              # join metadata to tree :contentReference[oaicite:14]{index=14}
     geom_tippoint(aes(size = RelAbundance, color = Phylum)) +        # point size ~ relative abundance, color ~ Phylum
     scale_size_continuous(range = c(1, 6)) +                         # adjust point-size range
     theme_tree2() +                                                  # clean tree theme
     labs(color = "Phylum", size = "Rel. Abundance (%)")

p

p <- ggtree(phy_tree) %<+% otu_metadata +
     geom_tippoint(aes(size = RelAbundance, color = Phylum)) +
     geom_tiplab(aes(label = Genus), size = 2) +
     scale_size_continuous(range = c(1, 6)) +
     theme_tree2() +
     labs(color = "Phylum", size = "Rel. Abundance (%)")

p

p <- ggtree(phy_tree,branch.length = "none" ) %<+% otu_metadata +                              # join metadata to tree :contentReference[oaicite:14]{index=14}
     geom_tippoint(aes(size = RelAbundance, color = Phylum)) +        # point size ~ relative abundance, color ~ Phylum
     scale_size_continuous(range = c(1, 6)) +                         # adjust point-size range
     theme_tree2() +                                                  # clean tree theme
     labs(color = "Phylum", size = "Rel. Abundance (%)")

gheatmap(p, log(otu_mat))


p <- ggtree(phy_tree, branch.length = "none", layout = "circular") %<+% otu_metadata +                              # join metadata to tree :contentReference[oaicite:14]{index=14}
     geom_tippoint(aes(size = RelAbundance, color = Phylum)) +        # point size ~ relative abundance, color ~ Phylum
     scale_size_continuous(range = c(1, 6)) +                         # adjust point-size range
     theme_tree2() +                                                  # clean tree theme
     labs(color = "Phylum", size = "Rel. Abundance (%)")

p
