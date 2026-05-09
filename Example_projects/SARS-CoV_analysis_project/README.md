# SARS-CoV-2 Phylogeny Analysis

Small educational R project for visualizing a SARS-CoV-2 phylogenetic tree together with rich metadata (clade, lineage, spike mutations, region, host, sampling date). The goal is to teach the six principled visual encodings of phylogenetic metadata in `ggtree` on a real dataset where the biological story is well-known and the failure modes (sampling bias, naive ancestral state reconstruction, branch-unit confusion) are pedagogically valuable.

See [`project_description.md`](project_description.md) for the full background, hypothesis, and research questions.

## Structure

- `data/` - input files: time-calibrated Newick tree, Nextstrain-style metadata TSV, optional spike RBD alignment
- `scripts/` - R scripts for each analysis step (create as you go)
- `results/` - CSV/TSV intermediate tables (clade counts, mutation matrix, per-country sampling)
- `reports/` - Quarto report and rendered figures

`scripts/`, `results/`, and `reports/` do not exist yet — create them when the workflow needs them.

## Tools and style

- **tidyverse** for data handling
- **ggplot2** + **ggtree** + **ggtreeExtra** + **ggnewscale** for plotting
- **treeio**, **ape**, **phytools**, **Biostrings** for tree and alignment manipulation
- **Quarto** for the final report

Keep code short, clear, and functional. One script per logical step.

## Notes

- Data acquisition is documented in [`nextstrain_howto.md`](nextstrain_howto.md). The expected files are pre-cached in `data/` — do not refetch during analysis.
- Tooling and editor setup is in [`setup.md`](setup.md).
- Agent-specific guidance (project structure, R style, what to avoid) is in [`AGENTS.md`](AGENTS.md).
- A worked seed prompt for the coding agent is in [`example_prompt.md`](example_prompt.md).
