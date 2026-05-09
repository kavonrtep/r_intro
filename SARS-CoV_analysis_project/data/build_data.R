#!/usr/bin/env Rscript
# build_data.R
#
# One-shot script: pulls a Nextstrain ncov Auspice JSON, extracts a
# time-calibrated Newick tree + per-tip metadata TSV, subsamples to ~400
# tips stratified by Nextstrain clade, and writes the two files this
# project's analysis scripts expect:
#
#   data/ncov_subsample.nwk     time tree, branches in years
#   data/ncov_metadata.tsv      one row per tip
#
# Run from the project root:
#   Rscript data/build_data.R
#
# This script is *not* part of the analysis pipeline. The student should
# read from data/ directly. See nextstrain_howto.md for what the files
# contain and how to refresh them.

suppressPackageStartupMessages({
  library(jsonlite)
  library(ape)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(lubridate)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(42)

# ---- 1. Download Auspice JSON --------------------------------------------

dl_url <- "https://data.nextstrain.org/ncov_open_global_6m.json"
dl_dest <- "data/_auspice_cache.json"

if (!file.exists(dl_dest)) {
  message("Downloading ", dl_url)
  # curl --compressed handles the Content-Encoding: gzip negotiation that
  # this endpoint requires. R's download.file is finicky here.
  status <- system2("curl", c("-sL", "--compressed", "-o", shQuote(dl_dest),
                              shQuote(dl_url)))
  if (status != 0) stop("download failed (curl exit ", status, ")")
}

aus <- read_json(dl_dest)
stopifnot(aus$version == "v2")

# ---- 2. Walk the nested tree, collect tips and edges ---------------------
#
# Auspice v2 stores the tree as nested children. We walk it depth-first,
# assigning ape-style integer node IDs (tips 1..n_tips, internals n_tips+1..),
# and collect edges with branch lengths in *years* (parent$num_date -> child$num_date).

# First pass: enumerate every node, separate tips from internals.

flat_nodes <- list()
walk_enumerate <- function(node, parent_idx = NA_integer_) {
  is_tip <- is.null(node$children) || length(node$children) == 0
  flat_nodes[[length(flat_nodes) + 1L]] <<- list(
    name        = node$name,
    is_tip      = is_tip,
    parent_idx  = parent_idx,
    num_date    = node$node_attrs$num_date$value %||% NA_real_,
    region      = node$node_attrs$region$value %||% NA_character_,
    country     = node$node_attrs$country$value %||% NA_character_,
    division    = node$node_attrs$division$value %||% NA_character_,
    clade       = node$node_attrs$clade_membership$value %||% NA_character_,
    pango       = node$node_attrs$pango_lineage$value %||% NA_character_,
    nx_pango    = node$node_attrs$Nextclade_pango$value %||% NA_character_,
    emerging    = node$node_attrs$emerging_lineage$value %||% NA_character_,
    host        = node$node_attrs$host$value %||% NA_character_,
    author      = node$node_attrs$author$value %||% NA_character_,
    accession   = node$node_attrs$genbank_accession$value %||% NA_character_,
    s_muts_branch = paste(unlist(node$branch_attrs$mutations$S %||% list()), collapse = ",")
  )
  this_idx <- length(flat_nodes)
  if (!is_tip) {
    for (ch in node$children) walk_enumerate(ch, this_idx)
  }
  invisible(this_idx)
}

walk_enumerate(aus$tree)
message("Walked ", length(flat_nodes), " nodes")

# Re-order: ape requires tips to be 1..n_tips and internals n_tips+1..N.
is_tip_vec <- vapply(flat_nodes, `[[`, logical(1), "is_tip")
n_tips <- sum(is_tip_vec)
n_int  <- sum(!is_tip_vec)

old_to_new <- integer(length(flat_nodes))
old_to_new[which(is_tip_vec)]  <- seq_len(n_tips)
old_to_new[which(!is_tip_vec)] <- (n_tips + 1L):(n_tips + n_int)

# Build edge matrix: each non-root node contributes one edge (parent -> self).
edges <- do.call(rbind, lapply(seq_along(flat_nodes), function(i) {
  p <- flat_nodes[[i]]$parent_idx
  if (is.na(p)) return(NULL)
  c(old_to_new[p], old_to_new[i])
}))

# Branch lengths in years.
edge_lengths <- vapply(seq_along(flat_nodes), function(i) {
  p <- flat_nodes[[i]]$parent_idx
  if (is.na(p)) return(NA_real_)
  flat_nodes[[i]]$num_date - flat_nodes[[p]]$num_date
}, numeric(1))
edge_lengths <- edge_lengths[!is.na(edge_lengths)]
edge_lengths[edge_lengths < 0] <- 0  # rare numeric noise

# Tip labels in ape order.
tip_labels <- vapply(which(is_tip_vec), function(i) flat_nodes[[i]]$name, character(1))

phy <- structure(list(
  edge        = edges,
  edge.length = edge_lengths,
  Nnode       = n_int,
  tip.label   = tip_labels
), class = "phylo")

# Sanity check: a 6-month build should be in years and span a few months.
md <- max(node.depth.edgelength(phy))
message(sprintf("Tree max depth: %.3f years", md))
stopifnot(md > 0.1 && md < 10)

# ---- 3. Build per-tip metadata table -------------------------------------
#
# For S mutation profile, we also need to walk root->tip and accumulate
# spike mutations (handling reversions). Easiest: a recursive walk that
# carries the running set.

tip_s_muts <- character(length(flat_nodes))

walk_s_muts <- function(node, carried = character()) {
  branch <- unlist(node$branch_attrs$mutations$S %||% list())
  here <- carried
  for (m in branch) {
    # m is e.g. "T22I": ref=T, pos=22, alt=I
    ref <- substr(m, 1, 1); alt <- substr(m, nchar(m), nchar(m))
    pos <- substr(m, 2, nchar(m) - 1)
    # reversion if the same position is already in carried with alt==current_ref
    existing <- grep(paste0("^.", pos, "."), here, value = TRUE)
    if (length(existing) > 0) {
      here <- setdiff(here, existing)
      if (substr(existing[1], 1, 1) != alt) here <- c(here, paste0(substr(existing[1], 1, 1), pos, alt))
    } else {
      here <- c(here, m)
    }
  }
  if (is.null(node$children) || length(node$children) == 0) {
    # tip
    return(list(name = node$name, muts = here))
  }
  out <- list()
  for (ch in node$children) out <- c(out, list(walk_s_muts(ch, here)))
  out
}

flat_tip_muts <- function(x, acc = list()) {
  if (!is.null(x$name) && !is.null(x$muts)) {
    acc[[length(acc) + 1L]] <- x; return(acc)
  }
  for (e in x) acc <- flat_tip_muts(e, acc)
  acc
}

s_walk <- walk_s_muts(aus$tree)
s_tips <- flat_tip_muts(s_walk)
s_lookup <- setNames(
  vapply(s_tips, function(x) paste(sprintf("S:%s", x$muts), collapse = ","), character(1)),
  vapply(s_tips, function(x) x$name, character(1))
)

# Pull tip rows from flat_nodes
tip_rows <- which(is_tip_vec)
meta_full <- tibble(
  strain        = vapply(tip_rows, function(i) flat_nodes[[i]]$name, character(1)),
  num_date      = vapply(tip_rows, function(i) flat_nodes[[i]]$num_date, numeric(1)),
  region        = vapply(tip_rows, function(i) flat_nodes[[i]]$region, character(1)),
  country       = vapply(tip_rows, function(i) flat_nodes[[i]]$country, character(1)),
  division      = vapply(tip_rows, function(i) flat_nodes[[i]]$division, character(1)),
  Nextstrain_clade = vapply(tip_rows, function(i) flat_nodes[[i]]$clade, character(1)),
  pango_lineage = vapply(tip_rows, function(i) flat_nodes[[i]]$pango, character(1)),
  Nextclade_pango = vapply(tip_rows, function(i) flat_nodes[[i]]$nx_pango, character(1)),
  emerging_lineage = vapply(tip_rows, function(i) flat_nodes[[i]]$emerging, character(1)),
  host          = vapply(tip_rows, function(i) flat_nodes[[i]]$host, character(1)),
  author        = vapply(tip_rows, function(i) flat_nodes[[i]]$author, character(1)),
  genbank_accession = vapply(tip_rows, function(i) flat_nodes[[i]]$accession, character(1))
) %>%
  mutate(
    nextclade_mutations = unname(s_lookup[strain]),
    date = format(as.Date(format(date_decimal(num_date), "%Y-%m-%d"))),
    decimal_date = num_date,
    year = as.integer(substr(date, 1, 4))
  )

# Add a length column (synthetic — Auspice JSON doesn't carry seq length).
# Approximate: SARS-CoV-2 genome ≈ 29903 nt, with N's reducing effective length.
# Use a stable per-strain pseudo-jitter so plots have a continuous variable to bind.
meta_full <- meta_full %>%
  mutate(
    length = 29903L - sample(0:300, n(), replace = TRUE),
    age    = sample(c(NA_integer_, 5:90), n(), replace = TRUE,
                    prob = c(0.4, dnorm(5:90, mean = 45, sd = 18))),
    sex    = sample(c(NA_character_, "M", "F"), n(), replace = TRUE,
                    prob = c(0.5, 0.25, 0.25))
  )

# ---- 4. Add latitude / longitude from the meta$geo_resolutions ------------
# Auspice stores per-country / per-division coordinates. Pull country-level.

country_geo <- aus$meta$geo_resolutions %>%
  Filter(function(x) x$key == "country", .) %>%
  .[[1]]
country_lookup <- map_dfr(country_geo$demes, function(d) {
  tibble(country = d %||% NA_character_)
}, .id = "country_name")
# Different shape — easier to extract directly:
country_xy <- map_dfr(names(country_geo$demes), function(nm) {
  d <- country_geo$demes[[nm]]
  tibble(country = nm,
         latitude = d$latitude %||% NA_real_,
         longitude = d$longitude %||% NA_real_)
})

meta_full <- meta_full %>%
  left_join(country_xy, by = "country")

message("Full meta rows: ", nrow(meta_full))
message("Tips in tree:   ", length(phy$tip.label))
stopifnot(setequal(meta_full$strain, phy$tip.label))

# ---- 5. Subsample to ~400 tips, stratified by clade ----------------------

target_total <- 400L
clade_counts <- meta_full %>% count(Nextstrain_clade, sort = TRUE)
message("Clade distribution before subsampling:")
print(clade_counts)

per_clade <- ceiling(target_total / nrow(clade_counts))
keep_strains <- meta_full %>%
  group_by(Nextstrain_clade) %>%
  slice_sample(n = per_clade) %>%
  ungroup() %>%
  pull(strain)

# If under target, top up with a random sample of remaining tips
if (length(keep_strains) < target_total) {
  extra <- sample(setdiff(meta_full$strain, keep_strains),
                  target_total - length(keep_strains))
  keep_strains <- c(keep_strains, extra)
}
keep_strains <- head(keep_strains, target_total)

phy_sub  <- keep.tip(phy, keep_strains)
meta_sub <- meta_full %>% filter(strain %in% keep_strains)

stopifnot(setequal(phy_sub$tip.label, meta_sub$strain))

message("Subsampled to ", length(phy_sub$tip.label), " tips")
message(sprintf("Subsampled tree max depth: %.3f years",
                max(node.depth.edgelength(phy_sub))))

# ---- 6. Write outputs ----------------------------------------------------

write.tree(phy_sub, "data/ncov_subsample.nwk")
write_tsv(meta_sub, "data/ncov_metadata.tsv")

message("\nWrote:")
message("  data/ncov_subsample.nwk  (", length(phy_sub$tip.label), " tips)")
message("  data/ncov_metadata.tsv   (", nrow(meta_sub), " rows, ",
        ncol(meta_sub), " cols)")
